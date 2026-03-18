import os
import sys
import glob
import numpy as np
import pandas as pd
from pathlib import Path


def get_SNPs_in_LD(chr_sumstat_df, ld_dir, snp, r4_threshold, p_threshold=0.05):
    """
    For a given SNP, finds all other SNPs (dominance_pval < p_threshold) in LD (>= r4).
    """
    # Parse the input SNP string (Format: "chr:pos:ref:alt")
    try:
        chrom, pos, ref, alt = snp.split(':')
        chrom = int(chrom)
        pos = int(pos)
    except ValueError:
        raise ValueError(f"SNP {snp} is not in the expected 'chr:pos:ref:alt' format.")

    # Locate the correct LD chunk file on your hard drive
    chunk_files = glob.glob(os.path.join(ld_dir, f"chr{chrom}_*.npz"))
    matching_chunks = []
    
    for f in chunk_files:
        basename = os.path.basename(f).replace('.npz', '')
        parts = basename.split('_')
        start = int(parts[1])
        end = int(parts[2])
        
        if start <= pos <= end:
            target_gz = os.path.join(ld_dir, f"{basename}.gz")
            matching_chunks.append((f, target_gz))
            
    if not matching_chunks:
        print(f"Error: No LD chunks containing position {pos} in chromosome {chrom} found.")
        return pd.DataFrame(), set()
    
    master_pos_to_r = {}
    
    for target_npz, target_gz in matching_chunks:
        if not os.path.exists(target_gz):
            continue

        manifest = pd.read_csv(target_gz, compression="gzip", sep='\t')
        pos_col = "position"

        snp_matches = manifest.index[manifest[pos_col] == pos].tolist()
        if not snp_matches:
            continue

        target_idx = snp_matches[0]

        with np.load(target_npz) as ld_data:
            row_indices = ld_data['row']
            col_indices = ld_data['col']
            correlation_data = ld_data['data']
        
        row_matches = (row_indices == target_idx)
        col_matches = (col_indices == target_idx)

        linked_indices = np.concatenate([col_indices[row_matches], row_indices[col_matches]])
        r_values = np.concatenate([correlation_data[row_matches], correlation_data[col_matches]])

        r4_values = (r_values**2)**2 # Dominance LD is square of additive LD
    
        # Apply 2-times r-squared filter
        valid_mask = (r4_values >= r4_threshold)
        final_indices = linked_indices[valid_mask]
        final_r4_values = r4_values[valid_mask]
    
        # Get the base-pair positions for these linked SNPs
        ld_positions = manifest.iloc[final_indices][pos_col].values

        # Zip them together and update the master dictionary
        # If a position already exists, it harmlessly overwrites it with the same correlation
        chunk_dict = dict(zip(ld_positions, final_r4_values))
        master_pos_to_r.update(chunk_dict)
    
    master_pos_to_r[pos] = 1.0


    # Final filtering on your summary stats dataframe
    result_df = chr_sumstat_df[
        (chr_sumstat_df['pos'].isin(master_pos_to_r.keys())) & 
        (chr_sumstat_df['dominance_pval'] < p_threshold)
    ].copy()

    ld_snps = [v for v in result_df['variant'].tolist() if v != snp]

    lead_snp_row = chr_sumstat_df[chr_sumstat_df['variant'] == snp].copy()
    result_df = pd.concat([result_df, lead_snp_row]).drop_duplicates(subset=['variant'])

    result_df['indep_status'] = (result_df['variant'] == snp)
    result_df['indep_id'] = snp
    result_df['r4'] = result_df['pos'].map(master_pos_to_r)
    result_df = result_df.sort_values(by='pos').reset_index(drop=True)
    
    return result_df, set(ld_snps)


def merge_ld_blocks(indep_df, lead_df,  merge_window=250):
    """Merges LD blocks closer than merge_window(kb)"""

    merge_window_kb = merge_window*1000

    lead_df = lead_df[["variant", "indep_status", "indep_id", "r4"]]
    lead_df = lead_df.rename(columns={
        "indep_status": "lead_status",
        "indep_id": "lead_id", 
        "r4": "r4_lead"
    })
    
    merged_df = indep_df.merge(lead_df, on="variant", how="left").copy()
    merged_df["lead_status"] = merged_df["lead_status"].fillna(False).astype(bool)
    merged_df["lead_id"] = merged_df.groupby(["chr", "indep_id"])["lead_id"].transform("first")

    merged_df = merged_df.sort_values(by=["chr", "pos"])
    merged_df = merged_df.reset_index(drop=True)
    
    merged_df["indep_start"] = merged_df.groupby(["chr", "indep_id"])["pos"].transform("min")
    merged_df["indep_end"] = merged_df.groupby(["chr", "indep_id"])["pos"].transform("max")
    
    # Isolate unique LD blocks
    blocks = merged_df[['chr', 'indep_id', 'indep_start', 'indep_end']].drop_duplicates()
    blocks = blocks.sort_values(by=['chr', 'indep_start'])
    print(f"8. Merging starts with {len(blocks)} unique LD blocks.")
    
    ld_starts = {}
    ld_ends = {}
    
    # Iterate through each chromosome and merge close intervals
    for chrom, chrom_group in blocks.groupby('chr'):
        current_start = None
        current_end = None
        current_indep_ids = []
        
        for _, row in chrom_group.iterrows():
            if current_start is None:
                # Initialize the first block for this chromosome
                current_start = row['indep_start']
                current_end = row['indep_end']
                current_indep_ids.append(row['indep_id'])
            else:
                # Check if the gap is smaller than or equal to the threshold
                if (row['indep_start'] - current_end) <= merge_window_kb:
                    # Merge them: push the end boundary out to cover the new block
                    current_end = max(current_end, row['indep_end'])
                    current_indep_ids.append(row['indep_id'])
                else:
                    # The gap is too big. Save the finalized super-block for all its members
                    for i_id in current_indep_ids:
                        ld_starts[i_id] = current_start
                        ld_ends[i_id] = current_end
                    
                    # Reset the tracker for the next independent block
                    current_start = row['indep_start']
                    current_end = row['indep_end']
                    current_indep_ids = [row['indep_id']]
                    
        # Catch the final block at the end of the chromosome loop
        if current_start is not None:
            for i_id in current_indep_ids:
                ld_starts[i_id] = current_start
                ld_ends[i_id] = current_end
                
    # Map the newly calculated super-boundaries back to every SNP in the main dataframe
    merged_df['ld_start'] = merged_df['indep_id'].map(ld_starts)
    merged_df['ld_end'] = merged_df['indep_id'].map(ld_ends)

    cols_to_join = ["chr", "ld_start", "ld_end"]
    merged_df["ld_id"] = merged_df[cols_to_join].astype(str).agg(':'.join, axis=1)
    print(f"9. {merged_df['ld_id'].nunique()} independent LD blocks remained after merging.")

    return merged_df


def main(sumstat_path, ld_dir, phen_code, output_dir, p_threshold=(5e-8)/1060):

    # Read GWAS summary statistic file
    print("1. Reading GWAS summary statistic file...")
    sumstat_df = pd.read_csv(
        sumstat_path,
        sep="\t",
        compression="gzip",
        dtype={"variant":str, "rsid":str, "chr":int, "pos":int,
            "minor_AF":float, "low_confidence_variant":str, "N": int,
            "beta":float, "pval":float, "dominance_beta":float,
            "dominance_pval":float, "add_sig":int, "dom_sig":int
        }
    )

    # Get chromosomes
    chrom = sumstat_df["chr"].dropna().unique()
    print(f"2. {len(chrom)} chromosome(s) with at least one significant variant identified:", chrom)

    # Initialize output dataframe
    out_df = pd.DataFrame()
    indep_df = pd.DataFrame()
    lead_df = pd.DataFrame()


    # Loop through chromosomes
    for c in chrom:
        print(f"3. Processing chromosome {c}...")
        chr_df = sumstat_df[sumstat_df["chr"] == c].copy()

        # Get significant SNPs and sort them
        sorted_sig_df = chr_df[chr_df["dominance_pval"] < p_threshold].sort_values("dominance_pval", ascending=True).copy()
        print(f"{len(sorted_sig_df)} significant SNPs (p < {p_threshold}) are found.")
        
        # Initialize variables
        indep_sig_snps = []
        stage1_clumped = set()
        r4_threshold_1=(0.6)**2

        # Stage 1 clumping: loop through significant SNPs
        print(f"4. Identifying independent significant SNPs at r4 < {r4_threshold_1} and\n their LD blocks based on SNPs with p<0.05.")
        for _, row in sorted_sig_df.iterrows():
            snp = row['variant']
            if snp in stage1_clumped:
                continue
                
            indep_sig_snps.append(snp)
            
            # Find all other SNPs in LD to clump them
            snps_df, snps_set = get_SNPs_in_LD(chr_df, ld_dir, snp, r4_threshold_1, p_threshold=0.05)
            indep_df = pd.concat([indep_df, snps_df])
            stage1_clumped.update(snps_set)

        print(f"5. {len(indep_sig_snps)} independent significant SNPs are found and LD boundaries are identified.")

        # Stage 2 clumping : Define Lead SNPs
        stage2_df = chr_df[chr_df['variant'].isin(indep_sig_snps)].sort_values('dominance_pval')
       
        lead_snps = []
        stage2_clumped = set()
        r4_threshold_2=(0.1)**2

        print(f"6. Identifying lead SNPs at r4 < {r4_threshold_2}")
        for _, row in stage2_df.iterrows():
            snp = row['variant']
            if snp in stage2_clumped:
                continue
                
            lead_snps.append(snp)

            # Find all other SNPs in LD to clump them
            snps_df, snps_set = get_SNPs_in_LD(stage2_df, ld_dir, snp, r4_threshold_2, p_threshold=0.05)
            lead_df = pd.concat([lead_df, snps_df])
            stage2_clumped.update(snps_set)
        
        print(f"{len(lead_snps)} lead  SNPs are found.")

    # Merging the loci closer than 250 kb
    print("7. Merging the identified indepent LD blocks around the independent SNPs\n if the blocks are closer than 250 kb")
    out_df = merge_ld_blocks(indep_df, lead_df)
   
    # Save files
    print("10. Saving files...")
    out_path = f"{output_dir}/{phen_code}_sig_loci.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)
    #out_path = f"{output_dir}/{phen_name}_sig_loci.xlsx"
    #out_df.to_excel(out_path, index=False)
    print(f"Processing for {phen_code} is done!")


# ==========================================
# ==========================================

if __name__ == "__main__":
    
    # Check if a task ID was passed from the terminal
    # if len(sys.argv) < 2:
    #     print("Error: missing arguments")
    #     sys.exit(1)

    # task_id = int(sys.argv[1]) - 1
    # phen_codes = sys.argv[2]
    # ld_dir = sys.argv[3]  
    # sumstats_dir = sys.argv[4]
    # out_dir = sys.argv[5]
    
    # phen_list = pd.read_excel(phen_codes, usecols=["phenotype_code"])["phenotype_code"].sort_values(ascending=True).tolist()

    # if task_id < 0 or task_id >= len(phen_list):
    #     print(f"Stopping execution: task_id {task_id} is out of bounds (phenotype list size is {len(phen_list)}).")
    #     sys.exit(0)

    # sumstat_path =f"{sumstats_dir}/{phen_list[task_id]}_sig_SNPs.tsv.bgz"

    # print(f"Process starts for phenotype {phen_list[task_id]}")
    # main(sumstat_path, ld_dir, phen_list[task_id], out_dir)

    sumstat_path="/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sumstats_QCed/M72_sig_SNPs.tsv.bgz"
    phen_code="M72"
    ld_dir= "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/ld_files"
    out_dir= "/Users/sezgi/Documents/dominance_pleiotropy/loci_level"
    main(sumstat_path, ld_dir, phen_code, out_dir)

