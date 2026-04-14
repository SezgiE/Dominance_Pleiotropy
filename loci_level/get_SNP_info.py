import pandas as pd
import numpy as np


def std_beta(dom_z, maf, n):
    
    """Standardizes the z value based on MAF and sample size"""

    std_b = dom_z / np.sqrt(2 * maf * (1 - maf) * (n + dom_z**2))
    
    return std_b


def get_snp_info(snp, phen_code, sumstats_dir):

    sumstats = pd.read_csv(f"{sumstats_dir}/{phen_code}_sig_SNPs.tsv.bgz", sep='\t',
                           compression='gzip',
                           usecols=["variant", "chr", "pos", "rsid", 
                                    "add_z_score", "dom_z_score", "minor_AF", "N",
                                    "add_log10_pval", "dom_log10_pval"]
                                   )
    
    add_log10_pval = sumstats[sumstats["variant"] == snp]["add_log10_pval"].values[0]
    dom_log10_pval = sumstats[sumstats["variant"] == snp]["dom_log10_pval"].values[0]
    maf = sumstats[sumstats["variant"] == snp]["minor_AF"].values[0]
    add_z = sumstats[sumstats["variant"] == snp]["add_z_score"].values[0]
    dom_z = sumstats[sumstats["variant"] == snp]["dom_z_score"].values[0]
    sample_size = sumstats[sumstats["variant"] == snp]["N"].values[0]

    return add_z, dom_z, maf, sample_size, add_log10_pval, dom_log10_pval


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

            add_z, dom_z, maf, sample_size, add_log10_pval, dom_log10_pval = get_snp_info(snp, phen_code, sumstats_dir)
            std_add_b = std_beta(add_z, maf, sample_size)
            std_dom_b = std_beta(dom_z, maf, sample_size)
            
            # Append the results to the DataFrame
            results.append({
                "variant": snp,
                "phen_code": phen_code,
                "phen_name": trait_info[trait_info["phenotype_code"] == phen_code]["description"].values[0],
                "category": trait_info[trait_info["phenotype_code"] == phen_code]["category"].values[0],
                "maf": maf,
                "add_z": add_z,
                "add_log10_pval": add_log10_pval,
                "dom_z": dom_z,
                "dom_log10_pval": dom_log10_pval,
                "sample_size": sample_size,
                "std_add_b": std_add_b,
                "std_dom_b": std_dom_b
            })


    result_df = pd.DataFrame(results)        
    result_df.to_csv("/Users/sezgi/Documents/dominance_pleiotropy/loci_level/loci_results/snp_info.tsv", sep='\t', index=False)
