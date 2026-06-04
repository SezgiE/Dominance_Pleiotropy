import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def fdr_correction(df, alpha_thresh=0.05):
    
    raw_pvals = df["p_value"]
    reject, pvals_corrected, _, _ = multipletests(raw_pvals, alpha=alpha_thresh, method='fdr_bh')
    df["fdr_corrected_p"] = pvals_corrected
    df["is_significant"] = reject

    return df


def fisher_test(df):

    cols_to_fill = ["mapped_in_group", "unmapped_in_group"]
    df[cols_to_fill] = df[cols_to_fill].fillna(0)

    total_mapped = df["mapped_in_group"].sum()
    total_unmapped = df["unmapped_in_group"].sum()

    # Mapped genes but NOT in this specific group
    df["mapped_out_group"] = total_mapped - df["mapped_in_group"]

    # Unmapped genes and NOT in this specific group            
    df["unmapped_out_group"] = total_unmapped - df["unmapped_in_group"] 
    
    def apply_fisher(row):
        table = [
            [row["mapped_in_group"], row["unmapped_in_group"]],
            [row["mapped_out_group"], row["unmapped_out_group"]]
        ]
        return pd.Series(fisher_exact(table, alternative='greater'))
    
    df[["odds_ratio", "p_value"]] = df.apply(apply_fisher, axis=1)
    
    return df


def process_data_enrich(background_phenotypes_path, pleio_phenotypes_path):
    
    background_ph = pd.read_csv(background_phenotypes_path, sep="\t")
    coloc_ph = pd.read_csv(pleio_phenotypes_path, sep="\t")

    pleio_ph_codes = set(coloc_ph['phen_code'].unique())
    background_ph['is_mapped'] = background_ph['phenotype_code'].isin(pleio_ph_codes)

    fisher_data = []
    for category in background_ph['category'].unique():
        subset = background_ph[background_ph['category'] == category]
        
        fisher_data.append({
            "group_name": category,
            "mapped_in_group": subset['is_mapped'].sum(),
            "unmapped_in_group": (~subset['is_mapped']).sum()
        })

    fisher_df = pd.DataFrame(fisher_data)

    return fisher_df


if __name__ == "__main__":

    background_phenotypes_path = "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/supp_table1.tsv"
    pleio_phenotypes_path = '/Users/sezgi/Documents/dominance_pleiotropy/loci_level/coloc_results/coloc_snp_info.tsv'
    output_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/loci_results"
    
    fisher_df = process_data_enrich(background_phenotypes_path, pleio_phenotypes_path)
    print(fisher_df)
    fisher_df = fisher_test(fisher_df)
    summary_df = fdr_correction(fisher_df)
    summary_df.to_csv(f"{output_dir}/category_enrich_summary.tsv", sep="\t", index=False)
