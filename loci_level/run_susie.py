import os
import sys
import glob
import tempfile
import subprocess
import numpy as np
import pandas as pd


def get_ld_matrix(snp_list, ld_dir):
    """
    Takes a list of SNPs and returns an NxN pairwise LD correlation matrix (r).
    """
    if not snp_list:
        return np.array([]).reshape(0, 0), None

    chroms = set([s[0] for s in snp_list])
    if len(chroms) > 1:
        raise ValueError("All SNPs in the locus must be on the same chromosome.")
    chrom = chroms.pop()

    positions_tmp = sorted([s[1] for s in snp_list])

    # Identify locus boundaries
    min_pos = positions_tmp[0]
    max_pos = positions_tmp[-1]

    # Locate the single encompassing LD chunk
    chunk_files = glob.glob(os.path.join(ld_dir, f"chr{chrom}_*.npz"))
    target_npz = None
    target_gz = None

    for npz_file in chunk_files:
        basename = os.path.basename(npz_file).replace(".npz", "")
        parts = basename.split("_")
        start = int(parts[1])
        end = int(parts[2])

        # Check if the chunk fully contains our locus range
        if start <= min_pos and end >= max_pos:
            target_npz = npz_file
            target_gz = os.path.join(ld_dir, f"{basename}.gz")
            break

    if not target_npz or not os.path.exists(target_gz):
        print(f"No single LD chunk found covering the range {min_pos} - {max_pos}.")
        return np.array([]).reshape(0, 0), None

    # Load manifest file for LD data
    manifest = pd.read_csv(target_gz, compression="gzip", sep="\t")
    pos_col = "position"

    manifest_subset = manifest[manifest[pos_col].isin(positions_tmp)]

    # Remove multiallelic variants
    manifest_subset = manifest_subset.drop_duplicates(subset=[pos_col], keep=False)

    if manifest_subset.empty:
        return np.array([]).reshape(0, 0), None

    # Map to chunk indices
    positions = sorted(manifest_subset[pos_col].tolist())
    pos_to_idx = {pos: i for i, pos in enumerate(positions)}
    chunk_idx_to_pos = dict(zip(manifest_subset.index, manifest_subset[pos_col]))
    chunk_indices = list(chunk_idx_to_pos.keys())

    # Initiate an identity matrix
    N = len(positions)
    R_matrix = np.eye(N)

    # Load the LD matrix and filter
    with np.load(target_npz) as ld_data:
        row_indices = ld_data["row"]
        col_indices = ld_data["col"]
        correlation_data = ld_data["data"]

    valid_mask = np.isin(row_indices, chunk_indices) & np.isin(
        col_indices, chunk_indices
    )

    valid_rows = row_indices[valid_mask]
    valid_cols = col_indices[valid_mask]
    valid_data = correlation_data[valid_mask]

    # Populate the NxN matrix symmetrically
    for r, c, val in zip(valid_rows, valid_cols, valid_data):
        p1 = chunk_idx_to_pos[r]
        p2 = chunk_idx_to_pos[c]
        i = pos_to_idx[p1]
        j = pos_to_idx[p2]

        R_matrix[i, j] = val
        R_matrix[j, i] = val

    # Sorting
    sort_order = np.argsort(positions)
    matrix_manifest = np.array(positions)[sort_order].tolist()
    ld_matrix = R_matrix[sort_order, :][:, sort_order]

    # Make sure its symmetrical
    np.fill_diagonal(ld_matrix, 1.0)
    ld_matrix = (ld_matrix + ld_matrix.T) / 2.0

    # Dominance LD
    ld_matrix_2 = ld_matrix ** 2

    # Make sure no NAs
    if np.isnan(ld_matrix_2).any():
        print("Warning: NaNs detected in LD matrix. Neutralizing to 0.")
        ld_matrix_2 = np.nan_to_num(ld_matrix_2, nan=0.0)

    return ld_matrix_2, matrix_manifest


def get_data(phen_code, sumstat_dir, loci_dir):
    """Gets  data and filtering out the MHC region."""

    MHC_chr = 6
    MHC_start = 25000000
    MHC_end = 34000000

    # Load Data
    raw_df = pd.read_csv(
        f"{sumstat_dir}/{phen_code}_sig_SNPs.tsv.bgz", sep="\t", compression="gzip"
    )
    loci_df = pd.read_csv(
        f"{loci_dir}/{phen_code}_sig_loci.tsv",
        sep="\t",
        usecols=["variant", "ld_start", "ld_end", "ld_id"],
    )

    if loci_df.empty:
        print("No loci to plot.")
        return None

    # Mask Raw SNPs and omit the one in the MCH regions
    mhc_snps = (
        (raw_df["chr"] == MHC_chr)
        & (raw_df["pos"] >= MHC_start)
        & (raw_df["pos"] <= MHC_end)
    )
    masked_df = raw_df[~mhc_snps].reset_index(drop=True)

    if masked_df.empty:
        return None

    # Process and Mask Loci Blocks
    unique_blocks = loci_df[["ld_id"]].dropna().drop_duplicates().reset_index(drop=True)
    unique_blocks[["chr", "start_bp", "end_bp"]] = (
        unique_blocks["ld_id"].str.split(":", expand=True).astype(int)
    )

    # Filtering out the LD block in the MHC region
    mhc_blocks = (
        (unique_blocks["chr"] == MHC_chr)
        & (unique_blocks["start_bp"] <= MHC_end)
        & (unique_blocks["end_bp"] >= MHC_start)
    )

    masked_blocks = unique_blocks[~mhc_blocks].reset_index(drop=True)

    if masked_blocks.empty:
        print(f"No locus to run SuSiE for '{phen_code}'. Skipping...")
        return None

    print(f"Found {len(masked_blocks)} unique loci to run SuSiE for.")

    return masked_df, masked_blocks


def run_SuSiE(
    sumstats_dir,
    loci_dir,
    sample_size,
    phen_code,
    ld_dir,
    r_script_path,
    output_dir,
    buffer_kb=50,
):
    """Executes SuSiE fine-mapping on independent loci."""

    # ------------------- Load data --------------------- #
    try:
        results = get_data(phen_code, sumstats_dir, loci_dir)

        if results is None:
            sys.exit(1)

        sumstat_data, unique_blocks = results

    except FileNotFoundError:
        print(f"No significant loci data for {phen_code} ")
        sys.exit(1)
    # ------------------- Load data --------------------- #
    

    all_susie_results = []

    # Create the output directory
    os.makedirs(output_dir, exist_ok=True)

    # Loop over the independen LD blocks
    for row in unique_blocks.itertuples(index=False):
        chrom = row.chr
        start_bp = row.start_bp
        end_bp = row.end_bp

        locus_string = f"{chrom}:{start_bp}:{end_bp}"

        # Define the LD  window (add the buffer)
        ld_start = start_bp - (buffer_kb * 1000)
        ld_end = end_bp + (buffer_kb * 1000)

        # Slice the LD block from the summary statistics
        block_df = (
            sumstat_data[
                (sumstat_data["chr"] == chrom)
                & (sumstat_data["pos"] >= ld_start)
                & (sumstat_data["pos"] <= ld_end)
            ]
            .copy()
            .sort_values(by="pos")
            .reset_index(drop=True)
        )

        # List of variant for which LD matrix will be obtained
        var_list = list(zip(block_df["chr"].astype(int), block_df["pos"].astype(int)))


        # Get the LD matrix
        print("Obtaining LD matrix...")
        try:
            ld_matrix, matrix_manifest = get_ld_matrix(var_list, ld_dir)
            if ld_matrix.size == 0 or matrix_manifest is None:
                continue
        except ValueError as e:
            print(f"Skipping locus {chrom}:{start_bp}-{end_bp} - {e}")
            continue


        # Remove the variants where  no LD information is available
        # Define the exact columns you want to keep
        cols_to_keep = ["variant", "chr", "pos", "minor_AF", "rsid", 
                        "dominance_beta", "dominance_se", "dom_z_score"]

        susie_df = (
            block_df.loc[block_df["pos"].isin(matrix_manifest), cols_to_keep]
            .copy()
            .sort_values(by="pos")
            .reset_index(drop=True)
        )


        # -------------------  Double Check Data Alignment ---------------------- #
        if (
            len(susie_df) != ld_matrix.shape[0]
            or susie_df["pos"].tolist() != matrix_manifest
        ):
            print(f"Alignment failed for locus {locus_string}. Skipping.")
            continue
        print(
            "The LD matrix is successfully obtained and aligned with summary statistics"
        )
        # -------------------  Double Check Data Alignment ---------------------- #


        #-------------------------- Run SuSie in R -------------------------------#
        print(
            f"Running SuSiE for locus {chrom}:{start_bp}-{end_bp} with {len(susie_df)} SNPs..."
        )

        # A temporary directory
        with tempfile.TemporaryDirectory() as tmpdir:
            susie_df_file = os.path.join(tmpdir, "susie_df.tsv")
            ld_file = os.path.join(tmpdir, "ld_matrix.csv")
            out_file = os.path.join(tmpdir, "susieR_results.tsv")

            # Write the data to disk for R to read
            susie_df.to_csv(susie_df_file, sep='\t', index=False)
            np.savetxt(ld_file, ld_matrix, delimiter=",")

            # Call the R script via the command line
            command = [
                "Rscript", 
                r_script_path, 
                susie_df_file, 
                ld_file, 
                str(sample_size), 
                out_file
            ]

            try:
                # Run the command and wait for it to finish
                subprocess.run(command, check=True, capture_output=True, text=True)
                
                # Read the results back from R
                r_results = pd.read_csv(out_file, sep='\t')
                
                # Merge all 7 fine-mapping columns back into your main dataframe
                new_columns = ["PIP", "CS", "CS_prob", "low_purity", "lead_r2", "post_mean", "post_sd", "lambda"]
                
                for col in new_columns:
                    susie_df[col] = r_results[col].values

                susie_df.insert(0, "phen_id", phen_code)
                susie_df.insert(1, "locus_id", locus_string)
                
                # Append to master list
                all_susie_results.append(susie_df)

            except subprocess.CalledProcessError as e:
                print(f"SuSiE R script failed for this locus.")
                print(f"R Error Output:\n{e.stderr}")
                continue

            #-------------------------- Run SuSie in R -------------------------------#


        # Save the final locus-specific dataframe and locus matrix
        os.makedirs(f"{output_dir}/susie_raw_files", exist_ok=True)
        np.savetxt(f"{output_dir}/susie_raw_files/{phen_code}_{chrom}:{start_bp}:{end_bp}_matrix.csv", ld_matrix, delimiter=",")
        susie_df.to_csv(f"{output_dir}/susie_raw_files/{phen_code}_{chrom}:{start_bp}:{end_bp}_data.csv", sep='\t', index=False)

    print("Merging all loci into master dataframe...")
    susie_res_df = pd.concat(all_susie_results, ignore_index=True)
    
    final_out = os.path.join(output_dir, f"{phen_code}_susie_res.tsv")
    susie_res_df.to_csv(final_out, sep='\t', index=False)
    print(f"Success! Saved results to {final_out}")

    
    return all_susie_results


if __name__ == "__main__":

    if len(sys.argv) < 4:
        print("Error: missing arguments")
        sys.exit(1)

    task_id = int(sys.argv[1])-1
    sumstats_dir = sys.argv[2]
    phen_dict_path = sys.argv[3]
    loci_dir = sys.argv[4]
    ld_dir = sys.argv[5]
    r_script_path = sys.argv[6]
    out_dir = sys.argv[7]

    # sumstats_dir = (
    #     "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sumstats_QCed"
    # )
    # phen_dict_path = (
    #     "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/phen_dict.xlsx"
    # )
    # loci_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sig_loci"
    # ld_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/ld_files"
    # r_script_path = "/Users/sezgi/Documents/dominance_pleiotropy/scripts/loci_level/run_SuSie_core.R"
    # out_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/susie_results"

    traits = pd.read_excel(phen_dict_path, 
                           usecols=["phenotype_code", "description", "n_non_missing"]
                           )

    traits = traits.sort_values(by="phenotype_code").reset_index(drop=True)
    row = traits.iloc[task_id]
    # traits= traits[traits["phenotype_code"] == "1747_2"]
    # row = traits.iloc[0]

    col_names = ["Phenotype", "Locus", "Total_sig", "Multi-allelic_sig", "Ratio"]
    df_dup = pd.DataFrame(columns=col_names)
    mult_SNPs = set()
    print(row)

    phen_code = row["phenotype_code"]
    phen_name = row["description"]
    sample_size = int(row["n_non_missing"])

    print(f"Processing: {phen_name} ({phen_code})")
    
    run_SuSiE(sumstats_dir, loci_dir, sample_size, phen_code, ld_dir, r_script_path,
              out_dir)
