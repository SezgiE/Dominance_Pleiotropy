import os
import sys
import glob
import numpy as np
import pandas as pd


def get_ld_matrix(snp_list, ld_dir):
    """
    Takes a list of SNPs and returns an NxN pairwise LD correlation matrix (r).
    """

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
        raise ValueError(
            f"No single LD chunk found covering the range {min_pos} - {max_pos}."
        )

    # Load manifest and map to chunk indices
    manifest = pd.read_csv(target_gz, compression="gzip", sep="\t")
    pos_col = "position"

    manifest_subset = manifest[manifest[pos_col].isin(positions_tmp)]

    if manifest_subset.empty:
        raise ValueError(f"The LD chunk does no include variants in the list.")

    positions = sorted(manifest_subset[pos_col].tolist())
    N = len(positions)
    pos_to_idx = {pos: i for i, pos in enumerate(positions)}
    R_matrix = np.eye(N)

    chunk_idx_to_pos = dict(zip(manifest_subset.index, manifest_subset[pos_col]))
    chunk_indices = list(chunk_idx_to_pos.keys())

    # Load sparse matrix and filter (Single read)
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


def run_SuSiE(sumstat_data, unique_blocks, sample_size, phen_code,
              phen_name, ld_dir, output_dir, buffer_kb=50):
    
    """Generates LocusZoom-style plots for each unique independent locus."""

    os.makedirs(output_dir, exist_ok=True)

    for row in unique_blocks.itertuples(index=False):
        chrom = row.chr
        start_bp = row.start_bp
        end_bp = row.end_bp

        # Define the LD  window (add the buffer)
        ld_start = start_bp - (buffer_kb * 1000)
        ld_end = end_bp + (buffer_kb * 1000)

        # Slice the specific window for the LD block
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

        var_list = list(zip(block_df["chr"].astype(int), block_df["pos"].astype(int)))

        ld_matrix, matrix_manifest = get_ld_matrix(var_list, ld_dir)

        susie_df = (
            block_df[block_df["pos"].isin(matrix_manifest)]
            .copy()
            .sort_values(by="pos")
            .reset_index(drop=True)
        )

        susie_df = susie_df[["dom_z_score"]]
        print(len(block_df))
        print(len(susie_df))
        print(ld_matrix.shape)
        print(susie_df.head())
        print(len(set(matrix_manifest)))

    return


if __name__ == "__main__":

    sumstat_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sumstats_QCed"
    phen_dict_path = "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/phen_dict_renamed.xlsx"
    loci_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sig_loci"
    ld_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/ld_files"
    out_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/susie_results"

    traits = pd.read_excel(
        phen_dict_path, usecols=["phenotype_code", "description", "n_non_missing"]
    )
    traits = traits[traits["phenotype_code"] == "1747_2"]

    for index, row in traits.iterrows():

        phen_code = row["phenotype_code"]
        phen_name = row["description"]
        sample_size = int(row["n_non_missing"])

        print(f"Processing: {phen_code}")

        try:
            results = get_data(phen_code, sumstat_dir, loci_dir)

            if results is None:
                continue  # Skip to the next phenotype

            sumstat_data, unique_blocks = results

        except FileNotFoundError:
            # Skip to the next phenotype if the file/dir doesn't exist
            print(f"No significant loci data for{phen_code} ")
            continue

        run_SuSiE(
            sumstat_data, unique_blocks, sample_size, phen_code, phen_name, ld_dir, out_dir
        )
