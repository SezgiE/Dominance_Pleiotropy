import os
import sys
import glob
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter

# Import the required R packages
base = importr('base')
susieR = importr('susieR')


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

    return ld_matrix, matrix_manifest


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
    sumstat_data,
    unique_blocks,
    sample_size,
    phen_code,
    phen_name,
    ld_dir,
    output_dir,
    buffer_kb=50,
):
    """Executes SuSiE fine-mapping on independent loci."""
    
    all_susie_results = []

    # Create the output directory
    os.makedirs(output_dir, exist_ok=True)

    # Loop over the independen LD blocks
    for row in unique_blocks.itertuples(index=False):
        chrom = row.chr
        start_bp = row.start_bp
        end_bp = row.end_bp

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
        except ValueError as e:
            print(f"Skipping locus {chrom}:{start_bp}-{end_bp} - {e}")
            continue

        if ld_matrix.size == 0 or matrix_manifest is None:
            continue

        # Remove the variants where  no LD information is available
        susie_df = (
            block_df[block_df["pos"].isin(matrix_manifest)]
            .copy()
            .sort_values(by="pos")
            .reset_index(drop=True)
        )

        # --- THE DOUBLE CHECK STEP FOR DATA ALIGNMENT ---
        if (
            len(susie_df) != ld_matrix.shape[0]
            or susie_df["pos"].tolist() != matrix_manifest
        ):
            print(f"Alignment failed for locus {chrom}:{start_bp}-{end_bp}. Skipping.")
            continue
        # -----------------------------
        print(
            "The LD matrix is successfully obtained and aligned with summary statistics"
        )

        # Run SuSie     
        z_scores = susie_df["dom_z_score"].values

        print(
            f"Running SuSiE for locus {chrom}:{start_bp}-{end_bp} with {len(z_scores)} SNPs..."
        )

        if np.isnan(ld_matrix).any():
            print("Warning: NaNs detected in LD matrix. Neutralizing to 0.")
            ld_matrix = np.nan_to_num(ld_matrix, nan=0.0)
        
        ld_matrix_2 = ld_matrix ** 2
        np.savetxt(f"{output_dir}/{chrom}:{start_bp}:{end_bp}_matrix.csv", ld_matrix_2, delimiter=",")
        susie_df.to_csv(f"{output_dir}/{chrom}:{start_bp}:{end_bp}_data.csv", sep='\t', index=False)
        
        # Execute SuSiE
        try:
    
            with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
                
                # Run the model
                susie_fit = susieR.susie_rss(
                    z=z_scores, 
                    R=ld_matrix_2, 
                    n=sample_size, 
                    L=10,
                    estimate_residual_variance=False
                )
                
                # Extract PIPs (automatically converts to a safe Python array)
                pips = np.array(susie_fit.rx2("pip"))
                susie_df["PIP"] = pips
                
                # Extract Credible Sets
                susie_df["CS"] = 0 
                cs_res = susieR.susie_get_cs(susie_fit)
                
                if cs_res != ro.rinterface.NULL and "cs" in cs_res.names:
                    cs_list = cs_res.rx2("cs")
                    for i, cs_name in enumerate(cs_list.names):
                        r_indices = np.array(cs_list[i])
                        py_indices = [int(idx) - 1 for idx in r_indices] # R to Python index fix
                        cs_num = int(cs_name.replace('L', ''))
                        susie_df.loc[py_indices, "CS"] = cs_num

            # Append to master list
            all_susie_results.append(susie_df)
            
            # Save locus-specific results
            out_file = os.path.join(output_dir, f"{phen_code}_chr{chrom}_{start_bp}_{end_bp}_susie.tsv")
            susie_df.to_csv(out_file, sep='\t', index=False)

        except Exception as e:
            print(f"SuSiE R execution failed for this locus: {e}")
            continue
            
    return all_susie_results


if __name__ == "__main__":

    # if len(sys.argv) < 4:
    #     print("Error: missing arguments")
    #     sys.exit(1)

    # task_id = int(sys.argv[1])-1
    # sumstats_dir = sys.argv[2]
    # phen_dict_path = sys.argv[3]
    # loci_dir = sys.argv[4]
    # ld_dir = sys.argv[5]
    # out_dir = sys.argv[6]

    sumstats_dir = (
        "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sumstats_QCed"
    )
    phen_dict_path = (
        "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/phen_dict.xlsx"
    )
    loci_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sig_loci"
    ld_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/ld_files"
    sig_SNPs_path = "/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/significant_SNPs/all_sig_SNPs.tsv.gz"
    out_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/susie_results"

    traits = pd.read_excel(
        phen_dict_path, usecols=["phenotype_code", "description", "n_non_missing"]
    )

    traits = traits.sort_values(by="phenotype_code").reset_index(drop=True)
    # row = traits.iloc[task_id]
    traits= traits[traits["phenotype_code"] == "M72"]
    row = traits.iloc[0]

    col_names = ["Phenotype", "Locus", "Total_sig", "Multi-allelic_sig", "Ratio"]
    df_dup = pd.DataFrame(columns=col_names)
    mult_SNPs = set()
    print(row)

    phen_code = row["phenotype_code"]
    phen_name = row["description"]
    sample_size = int(row["n_non_missing"])

    print(f"Processing: {phen_code}")

    try:
        results = get_data(phen_code, sumstats_dir, loci_dir)

        if results is None:
            sys.exit(0)  # Skip to the next phenotype

        sumstat_data, unique_blocks = results

    except FileNotFoundError:
        # Skip to the next phenotype if the file/dir doesn't exist
        print(f"No significant loci data for {phen_code} ")
        sys.exit(0)
    
    run_SuSiE(sumstat_data, unique_blocks, sample_size, phen_code, phen_name, ld_dir,
              out_dir)
