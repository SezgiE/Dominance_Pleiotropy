import os
import sys
import glob
import tempfile
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict


def get_data_dict(data_dir):
    
    path = Path(data_dir)

    pheno_to_loci = defaultdict(list)

    for f in path.glob('*__*matrix.csv'):

        parts = f.stem.split('__')
        
        if len(parts) == 2:
            phen_id = parts[0]
            locus_id = parts[1].replace('_matrix', '')
    
            pheno_to_loci[phen_id].append(locus_id)
    
    return dict(pheno_to_loci)


def get_data(data_dir, phen_info, phen_code, locus_id):

    data_df = pd.read_csv(f"{data_dir}/{phen_code}__{locus_id}_data.csv", sep='\t')
    ld_matrix = pd.read_csv(f"{data_dir}/{phen_code}__{locus_id}_matrix.csv", header=None)

    N = phen_info[phen_info["phenotype_code"] == phen_code]["n_non_missing"].values[0]
    k_covariates = 25
    n_cases = phen_info[phen_info["phenotype_code"] == phen_code]["n_cases"].fillna(0).values[0]
    t_type = phen_info[phen_info["phenotype_code"] == phen_code]["variable_type"].values[0]
    
    prop_cases = n_cases / N
    trait_type = "quant" if t_type == "continuous_irnt" else "cc"
    
    return data_df, ld_matrix, N, k_covariates, prop_cases, trait_type


def check_overlap(locus1, locus2):
    """
    Returns True if two locus strings (chr_start_end) overlap.
    """
    chr1, start1, end1 = locus1.split('_')
    chr2, start2, end2 = locus2.split('_')
    
    # Must be on the same chromosome
    if chr1 != chr2:
        return False
        
    # Check for coordinate overlap (start1 <= end2 AND start2 <= end1)
    return int(start1) <= int(end2) and int(start2) <= int(end1)


def run_coloc(r_script_path, trait_1_df, ld1, n1, k1, prop_cases1, t1_type, trait_2_df, ld2, n2, k2, prop_cases2, t2_type):

    #-------------------------- Run Coloc in R -------------------------------#

    # A temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:

        # trait 1
        df1_file = os.path.join(tmpdir, "df1.tsv")
        ld1_file = os.path.join(tmpdir, "ld1.csv")
        trait_1_df.to_csv(df1_file, sep='\t', index=False)
        np.savetxt(ld1_file, ld1, delimiter=",")
        
        # trait 2
        df2_file = os.path.join(tmpdir, "df2.tsv")
        ld2_file = os.path.join(tmpdir, "ld2.csv")
        trait_2_df.to_csv(df2_file, sep='\t', index=False)
        np.savetxt(ld2_file, ld2, delimiter=",")
        
        # tmp output for coloc in R
        out_file = os.path.join(tmpdir, "coloc_results.tsv")


        # Call the R script via the command line
        command = [
            "Rscript",
            r_script_path,
            df1_file,
            ld1_file,
            str(n1),
            str(k1),
            str(prop_cases1),
            t1_type,
            df2_file,
            ld2_file,
            str(n2),
            str(k2),
            str(prop_cases2),
            t2_type,
            out_file   
        ]

        try:
            # Run the command
            subprocess.run(command, check=True, capture_output=True, text=True)

            try:
                coloc_results = pd.read_csv(out_file, sep='\t')

                return coloc_results

            except FileNotFoundError:
                print(f"Colocolization is not succesfull for this locus pair. No results to read. Skipping...")
                return None

            except pd.errors.EmptyDataError:
                print(f"Colocolization is not succesfull for this locus pair. No results to read. Skipping...")
                return None
        except subprocess.CalledProcessError as e:
            print(f"Coloc R script failed for this locus pair.")
            print(f"R Error Output:\n{e.stderr}")
            return None


def main(r_script_path, phen_info, data_dir, phen_code, out_dir):

    locus_ids = phen_to_loci.get(phen_code, [])
    print(f"Starting coloc for phenotype: {phen_code}")
    print(f"Number of loci for phenotype {phen_code}: {len(locus_ids)}")

    result_dfs = []

    for locus in locus_ids:

        print(f"Processing phenotype {phen_code} locus: {locus}")

        # Get data for locus 1
        trait_1_df, ld1, n1, k1, prop_cases1, trait1_type = get_data(data_dir, phen_info, phen_code, locus)


        # Scan the entire dictionary for overlaps
        for phen2, loci2_list in phen_to_loci.items():
            for locus2 in loci2_list:
                
                # Skip checking the exact same phenotype and locus against itself
                if phen2 == phen_code and locus2 == locus:
                    continue
                
                if check_overlap(locus, locus2):

                    print(f"Overlap found between '{phen_code} | locus: {locus}' and '{phen2} | locus: {locus2}'. Proceeding with coloc.")

                    # Get data for locus 2
                    trait_2_df, ld2, n2, k2, prop_cases2, trait2_type = get_data(data_dir, phen_info, phen2, locus2)

                    # Run coloc
                    print(f"Colocalizing:\nphen1: {phen_code} locus: {locus}\nphen2: {phen2} locus: {locus2}")

                    coloc_results = run_coloc(r_script_path, trait_1_df, ld1, n1, k1, prop_cases1, trait1_type, 
                                              trait_2_df, ld2, n2, k2, prop_cases2, trait2_type)
                    
                    if coloc_results is None:
                        print(f"Colocalization failed. Skipping...")
                        continue

                    # Save coloc results
                    coloc_results["phen1"] = phen_code
                    coloc_results["locus1"] = locus
                    coloc_results["phen2"] = phen2
                    coloc_results["locus2"] = locus2

                    # Adjust column order
                    front_cols = ["phen1", "locus1", "phen2", "locus2"]
                    other_cols = [col for col in coloc_results.columns if col not in front_cols]
                    
                    coloc_results = coloc_results[front_cols + other_cols]
                    
                    result_dfs.append(coloc_results)
                    print(f"Coloc completed for {locus} and {locus2}. Results added to the list.")
                
                else:
                    print(f"No overlap between '{phen_code} | locus: {locus}' and '{phen2} | locus: {locus2}'.\nSkipping coloc for this pair.")

    if result_dfs:
        # Concatenate all coloc results
        final_coloc_df = pd.concat(result_dfs, ignore_index=True)

        output_path = os.path.join(out_dir, f"{phen_code}_coloc_results.tsv")
        final_coloc_df.to_csv(output_path, index=False, sep='\t')


if __name__ == "__main__":

    # Define paths
    r_script_path = "/Users/sezgi/Documents/dominance_pleiotropy/scripts/loci_level/run_coloc_core.R"
    data_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/susie_results/susie_raw_files/"
    phen_info_path = "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/phen_dict_renamed.xlsx"
    out_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/coloc_results/results_by_phenotype"

    phen_info = pd.read_excel(phen_info_path, usecols=["phenotype_code", 
                                                       "description", "n_non_missing", 
                                                       "variable_type", "n_cases", "category"])

    phen_codes = phen_info["phenotype_code"].tolist()
    category_map = dict(zip(phen_info['phenotype_code'], phen_info['category']))
    phen_code_to_name = dict(zip(phen_info['phenotype_code'], phen_info['description']))
    
    #phen_to_loci = get_data_dict(data_dir)

    # for phen_code in phen_to_loci.keys():
        
    #     main(r_script_path, phen_info, data_dir, phen_code, out_dir)

    merged = pd.concat(pd.read_csv(f, sep='\t') for f in glob.glob(f'{out_dir}/*coloc_results.tsv'))
    merged = merged[(merged['cs_H4'] >= 0.8) & (merged['PIP'] >= 0.8)]

    merged["cat1"] = merged['phen1'].map(category_map)
    merged["phen_name1"] = merged['phen1'].map(phen_code_to_name)


    snp_info = merged.groupby('variant').agg(
        n_traits=('phen1', 'nunique'),
        n_categories=('cat1', lambda x: len(set(x))),
        phenotypes=('phen_name1', lambda x: ', '.join(x.unique())),
        phen_codes=('phen1', lambda x: ', '.join(x.unique()))
    ).reset_index()

    snp_info.to_csv(f'{out_dir}/snp_info.tsv', sep='\t', index=False)
    merged.to_csv(f'{out_dir}/merged_coloc.tsv', sep='\t', index=False)