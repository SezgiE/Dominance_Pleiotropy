import os
import sys
import time
import shutil
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path


def load_wget_commands(add_file_path, dom_file_path):
    # Load the additive and dominance Excel files
    a_df = pd.read_excel(add_file_path)
    d_df = pd.read_excel(dom_file_path)

    # Create dictionaries mapping phenotype_code -> wget command
    add_wget_dict = dict(zip(a_df["phenotype_code"], a_df["wget"]))
    dom_wget_dict = dict(zip(d_df["phenotype_code"], d_df["wget"]))

    # Convert keys to sets for highly efficient comparison
    a_codes_set = set(add_wget_dict.keys())
    d_codes_set = set(dom_wget_dict.keys())

    # Check  both files have the exact same phenotype codes
    if a_codes_set != d_codes_set:
        missing_in_d = a_codes_set - d_codes_set
        missing_in_a = d_codes_set - a_codes_set
        raise ValueError(
            f"Error: Phenotype codes do not match between files!\n"
            f"Codes missing in Dominance file: {missing_in_d}\n"
            f"Codes missing in Additive file: {missing_in_a}"
        )

        # Sort the phenotype codes for consistent processing order
    phenotype_codes = sorted(list(a_codes_set))

    return add_wget_dict, dom_wget_dict, phenotype_codes


def execute_wget(wget_cmd, temp_dir):
    """Executes wget and returns the exact path to the downloaded file."""
    # Extract filename defined after '-O '
    filename = wget_cmd.split("-O ")[-1].strip()
    filepath = os.path.join(temp_dir, filename)

    # Execute download securely
    subprocess.run(wget_cmd, shell=True, cwd=temp_dir, check=True)
    return filepath


def get_maf_filter(phenotype_info, phenotype_code):
        
        phen_info = pd.read_csv( phenotype_info,sep="\t",compression="gzip",
        usecols=[
            "phenotype",
            "variable_type"
        ]        )
        
        phen_info = phen_info[phen_info["phenotype"] == phenotype_code]

        if phen_info.empty:
            raise ValueError(f"Phenotype code {phenotype_code} not found in phenotype information file.")
        
        variable_type = phen_info["variable_type"].values[0]
        
        if variable_type == "binary" or variable_type == "categorical" or variable_type == "ordinal":
            maf_threshold = 0.05
        else:
            maf_threshold = 0.01
        
        print(f"Phenotype {phenotype_code} is of type '{variable_type}'.")
        return maf_threshold


def preprocess_sumstats(file_add, file_dom, file_out, code, maf_threshold, p_threshold):
    """
    Merges additive and dominance summary statistics and filters additive only or both additive
    and dominance variants.
    """

    print("1. Reading Additive GWAS SumStats...")
    data_add = pd.read_csv(
        file_add,
        sep="\t",
        compression="gzip",
        usecols=[
            "variant",
            "minor_AF",
            "low_confidence_variant",
            "n_complete_samples",
            "beta",
            "pval",
        ],
    ).rename(columns={"n_complete_samples": "N"})

    print("3. Reading Dominance GWAS SumStats...")
    data_dom = pd.read_csv(
        file_dom,
        sep="\t",
        compression="gzip",
        usecols=["variant", "dominance_beta", "dominance_pval"],
    )

    print("4. Merging data...")
    data_merge_raw = pd.merge(data_add, data_dom, on="variant", how="inner").copy()

    # Filter for MAF > 0.01 if continious or 0.05 if binary trait
    print(f"Applying MAF filter: minor_AF > {maf_threshold}...")
    data_merged = data_merge_raw[data_merge_raw["minor_AF"] > maf_threshold].copy()
    print(f"After MAF filtering, {len(data_merged)} variants remain. {len(data_merge_raw) - len(data_merged)} variants removed.")

    print(f"5. Filtering for genome-wide significant SNPs based on  p-value threshold of {p_threshold}...")
    # Assess significance based on p-value thresholds for both additive and dominance
    data_merged["add_sig"] = np.where(data_merged["pval"] < p_threshold, 1, 0)
    data_merged["dom_sig"] = np.where(data_merged["dominance_pval"] < p_threshold, 1, 0)

    # Filter for SNPs that are significant in additive
    data_out = data_merged[(data_merged["add_sig"] == 1)].copy()
    data_out[str(code)] = data_out["add_sig"] + data_out["dom_sig"]

    # Write to a gzipped file
    data_out.to_csv(file_out, sep="\t", index=False, compression="gzip")

    print(f"Done! Saved {len(data_out)} SNPs to {file_out}")


def run_single_trait(task_index, add_excel, dom_excel, phenotype_info, temp_dir, output_dir, p_threshold):
    """Executes the pipeline for ONE specific trait based on the SLURM array index."""

    os.makedirs(output_dir, exist_ok=True)
    add_dict, dom_dict, phenotype_codes = load_wget_commands(add_excel, dom_excel)

    # Grab ONLY the phenotype assigned to this specific array task
    try:
        code = phenotype_codes[task_index]
    except IndexError:
        print(f"Task index {task_index} is out of range. Exiting.")
        return

    # Create a UNIQUE temporary directory
    file_path = os.path.join(temp_dir, f"sumstats_{code}")
    os.makedirs(file_path, exist_ok=True)

    start_time = time.time()
    print(f"--- Starting workflow for phenotype: {code} ---")

    add_filepath, dom_filepath = None, None

    try:
        print("Downloading summary statistics...")
        start_time_download = time.time()
        add_filepath = execute_wget(add_dict[code], temp_dir)
        dom_filepath = execute_wget(dom_dict[code], temp_dir)
        time_download = time.time() - start_time_download
        print(f"Download completed in {time_download:.2f} seconds")

        out_filepath = os.path.join(output_dir, f"{code}_sig_SNPs.tsv.bgz")

        # Determine MAF threshold based on phenotype type
        maf_threshold = get_maf_filter(phenotype_info, code)
        
        print("Preprocessing sumstats...")
        time_preprocess_start = time.time()
        preprocess_sumstats(add_filepath, dom_filepath, out_filepath, code, maf_threshold, p_threshold)
        time_preprocess = time.time() - time_preprocess_start
        print(f"Preprocessing completed in {time_preprocess:.2f} seconds")

    except Exception as e:
        print(f"ERROR processing phenotype {code}: {e}")

    finally:
        # Delete the unique temp directory and everything inside it
        print("Cleaning up temporary sumstats...")
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

        elapsed_time = time.time() - start_time
        print(f"Finished cycle for {code} in {elapsed_time:.2f} seconds\n")


# ==========================================
# ==========================================

if __name__ == "__main__":
    # Check if a task ID was passed from the terminal
    if len(sys.argv) < 2:
        print("Usage: python script.py <task_index>")
        sys.exit(1)

    # SLURM arrays usually start at 1, but Python lists start at 0.
    # We subtract 1 to align them perfectly.
    task_index = int(sys.argv[1]) - 1
    tmp_dir = sys.argv[2]  # Grab the scratch path from bash
    out_dir = sys.argv[3]  # Grab the final save path from bash
    p_threshold = float(sys.argv[4])  # Grab the p-value threshold from bash
    
    if p_threshold is None:
        p_threshold = 5e-8  # Default to genome-wide significance if not provided

    # Get the directory where the script is located
    base_dir = Path(__file__).resolve().parent.parent

    # Define paths cleanly
    additive_excel = base_dir / "UKB_sumstats_Neale" / "a_sumStats.xlsx"
    dominance_excel = base_dir / "UKB_sumstats_Neale" / "d_sumStats.xlsx"
    phenotype_info = base_dir / "UKB_sumstats_Neale" / "all_phenotypes_info.tsv.gz"

    run_single_trait(task_index, additive_excel, dominance_excel, phenotype_info, tmp_dir, out_dir, p_threshold)

