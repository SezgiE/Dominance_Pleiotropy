import os
import glob
import numpy as np
import pandas as pd
import scipy.sparse

def get_SNPs_in_LD(chr_sumstat_df, ld_dir, snp, r2_threshold, p_threshold=0.05):
    """
    For a given SNP, finds all other SNPs (dominance_pval < p_threshold) in LD (>= r2).
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
    target_npz = None
    target_gz = None
    
    for f in chunk_files:
        basename = os.path.basename(f).replace('.npz', '')
        parts = basename.split('_')
        start = int(parts[1])
        end = int(parts[2])
        
        if start <= pos <= end:
            target_npz = f
            target_gz = os.path.join(ld_dir, f"{basename}.gz")
            break
            
    if not target_npz or not os.path.exists(target_gz):
        print(f"Error: LD chunk containing position {pos} not found in {ld_dir}.")
        return [],[]

    # Load the SNP manifest (.gz) to find our variant's exact row index
    manifest = pd.read_csv(target_gz, compression="gzip", sep='\t')
    pos_col = "position"
    
    snp_matches = manifest.index[manifest[pos_col] == pos].tolist()
    if not snp_matches:
        print(f"Warning: Position {pos} exists in the chunk window, but the SNP is missing from the LD reference panel.")
        return [],[]
        
    target_idx = snp_matches[0] # Grab the first match if multi-allelic

    with np.load(target_npz) as ld_data:
        row_indices = ld_data['row']
        col_indices = ld_data['col']
        correlation_data = ld_data['data']
        matrix_shape = ld_data['shape']
    
    # Extract only the correlations involving the target SNP
    row_matches = (row_indices == target_idx)
    col_matches = (col_indices == target_idx)

    linked_indices = np.concatenate([col_indices[row_matches], row_indices[col_matches]])
    r_values = np.concatenate([correlation_data[row_matches], correlation_data[col_matches]])
    r_values = (r_values**2)**2
    
    # Apply 2-times r-squared filter
    valid_mask = (r_values > r2_threshold)
    final_indices = linked_indices[valid_mask]
    final_r2_values = (r_values)[valid_mask]
    
    # Get the base-pair positions for these linked SNPs
    ld_positions = manifest.iloc[final_indices][pos_col].values
    pos_to_r2 = dict(zip(ld_positions, final_r2_values))
    pos_to_r2[pos] = 1.0

    # Final filtering on your summary stats dataframe
    result_df = chr_sumstat_df[
        (chr_sumstat_df['pos'].isin(ld_positions)) & 
        (chr_sumstat_df['dominance_pval'] < p_threshold)
    ]

    ld_snps = result_df['variant'].tolist()

    lead_snp_row = chr_sumstat_df[chr_sumstat_df['variant'] == snp].copy()
    result_df = pd.concat([result_df, lead_snp_row]).drop_duplicates(subset=['variant'])

    result_df['indep_status'] = (result_df['variant'] == snp)
    result_df['indep_id'] = snp
    result_df['r4'] = result_df['pos'].map(pos_to_r2)
    result_df = result_df.sort_values(by='pos').reset_index(drop=True)
    
    return result_df, set(ld_snps)


def merge_ld_blocks(ld_df, merge_window=250):
    """Merges LD blocks closer than merge_window(kb)"""
    print("Merging LD blocks")


def main(sumstat_path, ld_dir, output_path, p_threshold=(5e-8)/1060):

    # Read GWAS summary statistic file
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

    # Initialize output dataframe
    out_df = pd.DataFrame()
    indep_df = pd.DataFrame()
    lead_df = pd.DataFrame()


    # Loop through chromosomes
    for c in chrom:
        chr_df = sumstat_df[sumstat_df["chr"] == c].copy()

        # Get significant SNPs and sort them
        sorted_sig_df = chr_df[chr_df["dominance_pval"] < p_threshold].sort_values("dominance_pval", ascending=True).copy()
        
        # Initialize variables
        indep_sig_snps = []
        stage1_clumped = set()

        # Stage 1 clumping: loop through significant SNPs
        for _, row in sorted_sig_df.iterrows():
            snp = row['variant']
            if snp in stage1_clumped:
                continue
                
            indep_sig_snps.append(snp)
            
            # Find all other SNPs in LD (r2 >= 0.6) to clump them
            snps_df, snps_set = get_SNPs_in_LD(chr_df, ld_dir, snp, r2_threshold=(0.6)**2, p_threshold=0.05)
            indep_df = pd.concat([indep_df, snps_df])
            stage1_clumped.update(snps_set)


        # Stage 2 clumping : Define Lead SNPs (r2 < 0.1)
        stage2_df = chr_df[chr_df['variant'].isin(indep_sig_snps)].sort_values('dominance_pval')
       
        lead_snps = []
        stage2_clumped = set()

        for _, row in stage2_df.iterrows():
            snp = row['variant']
            if snp in stage2_clumped:
                continue
                
            lead_snps.append(snp)

            # Find all other SNPs in LD (r2 >= 0.1) to clump them
            snps_df, snps_set = get_SNPs_in_LD(stage2_df, ld_dir, snp, r2_threshold=(0.1)**2, p_threshold=0.05)
            lead_df = pd.concat([lead_df, snps_df])
            stage2_clumped.update(snps_set)
        
    # Save files
    out_name_indep = f"{output_path}/indep.xlsx"
    out_name_lead = f"{output_path}/lead.xlsx"
    indep_df.to_excel(out_name_indep, index=False)
    lead_df.to_excel(out_name_lead, index=False)


if __name__ == "__main__":
    
    sumstat_path = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sumstats_QCed/20002_1473_sig_SNPs.tsv.bgz"
    ld_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/ld_files"
    output_path = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level"
    
    main(sumstat_path, ld_dir, output_path)

