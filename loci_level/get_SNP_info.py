import numpy as np
import pandas as pd
from scipy.stats import norm


def std_beta(dom_z, maf, n):
    
    """Standardizes the z value based on MAF and sample size"""

    std_b = dom_z / np.sqrt(2 * maf * (1 - maf) * (n + dom_z**2))
    se = 1 / np.sqrt(2 * maf * (1 - maf) * (n + dom_z**2))
    
    return std_b, se


def get_snp_info(snp, phen_code, sumstats_dir):

    cols = [
        "variant", "minor_AF", "N", "chr", "rsid", "pos",
        "add_z_score", "dom_z_score", "add_log10_pval", "dom_log10_pval"
    ]
    
    sumstats = pd.read_csv(
        f"{sumstats_dir}/{phen_code}_sig_SNPs.tsv.bgz", 
        sep='\t',
        compression='gzip',
        usecols=cols
    )
    
    snp_row = sumstats[sumstats["variant"] == snp].iloc[0]
    
    # Split the snp string once
    snp_parts = snp.split(":")

    return {
        "variant": snp,
        "rsID": snp_row["rsid"],
        "chr": snp_row["chr"],
        "pos": snp_row["pos"],
        "A1": snp_parts[3],
        "A2": snp_parts[2],
        "maf": snp_row["minor_AF"],
        "sample_size": snp_row["N"],
        "add_z": snp_row["add_z_score"],
        "dom_z": snp_row["dom_z_score"],
        "add_log10_pval": snp_row["add_log10_pval"],
        "dom_log10_pval": snp_row["dom_log10_pval"]
    }


if __name__ == "__main__":

    sumstats_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sumstats_QCed"
    trait_info = pd.read_excel("/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/phen_dict_renamed.xlsx", 
                        usecols=["phenotype_code", "description", "category"]
                        )

    coloc_snps = pd.read_csv("/Users/sezgi/Documents/dominance_pleiotropy/loci_level/coloc_results/coloc_snps.tsv", sep='\t')
    snps = coloc_snps["variant"].unique()


    results = []
    
    for snp in snps:

        # Get the traits associated with the SNP
        phen_codes = coloc_snps[coloc_snps["variant"] == snp]["phen_codes"].str.split(", ").iloc[0]


        # Loop through each trait and get the SNP information
        for phen_code in phen_codes:

            snp_info = get_snp_info(snp, phen_code, sumstats_dir)
            std_add_b, std_add_se = std_beta(snp_info["add_z"], snp_info["maf"], snp_info["sample_size"])
            std_dom_b, std_dom_se = std_beta(snp_info["dom_z"], snp_info["maf"], snp_info["sample_size"])

            # Append the results to the DataFrame
            results.append({
                "phen_code": phen_code,
                "phen_name": trait_info[trait_info["phenotype_code"] == phen_code]["description"].iloc[0],
                "category": trait_info[trait_info["phenotype_code"] == phen_code]["category"].iloc[0],
                "variant": snp_info["variant"],
                "rsID": snp_info["rsID"],
                "CHR": snp_info["chr"],
                "BP": snp_info["pos"],
                "A1": snp_info["A1"],
                "A2": snp_info["A2"],
                "maf": snp_info["maf"],
                "sample_size": snp_info["sample_size"],
                "add_z": snp_info["add_z"],
                "dom_z": snp_info["dom_z"],
                "add_log10_pval": snp_info["add_log10_pval"],
                "dom_log10_pval": snp_info["dom_log10_pval"],
                "std_add_b": std_add_b,
                "std_add_se": std_add_se,
                "std_dom_b": std_dom_b,
                "std_dom_se": std_dom_se
            })


    # Save the results to a TSV file
    result_df = pd.DataFrame(results)        
    result_df.to_csv("/Users/sezgi/Documents/dominance_pleiotropy/loci_level/coloc_results/coloc_snp_info.tsv", sep='\t', index=False)


    # FUMA input file
    result_df = result_df.sort_values('dom_log10_pval', ascending=False).drop_duplicates(subset=['variant'], keep='first').copy()
    result_df["P"] = 2 * norm.sf(np.abs(result_df['std_dom_b'] / result_df['std_dom_se']))

    fuma_df = result_df[["rsID", "CHR", "BP", "A1", "A2", "P", "std_dom_b", "std_dom_se"]].copy()
    fuma_df = fuma_df.rename(columns={'std_dom_b': 'Beta', 'std_dom_se': 'SE'})

    fuma_independent_snps = fuma_df[["rsID", "CHR", "BP"]].copy()

    # Save the files
    fuma_path = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/fuma_input"
    fuma_df.to_csv(f'{fuma_path}/fuma_input.txt', sep='\t', index=False, encoding='ascii')
    fuma_independent_snps.to_csv(f'{fuma_path}/fuma_indp_snps.txt', sep='\t', index=False, encoding='ascii')

