import os
import pandas as pd

def filter_gtex(tissue_name, snps_eqtl_path, gene_eqtl_path, susie_eqtl_path, snps_pleio):

    snps_eqtl = pd.read_parquet(snps_eqtl_path)
    genes_eqtl = pd.read_csv(gene_eqtl_path, sep="\t", compression="gzip")
    susie_eqtl = pd.read_parquet(susie_eqtl_path)
    

    # Adjust SNP_b38 list
    snps_pleio.columns = ["chr", "start", "pos", "variant_id_b37"]
    snps_pleio = snps_pleio[["chr", "pos", "variant_id_b37"]]


    # Adjust SNPs and Susie eqtl
    susie_eqtl[["chr", "pos", "ref", "alt", "build"]] = susie_eqtl["variant_id"].str.split("_", expand=True)
    susie_eqtl = susie_eqtl.drop(columns=["variant_id"])

    snps_eqtl[["chr", "pos", "ref", "alt", "build"]] = snps_eqtl["variant_id"].str.split("_", expand=True)
    snps_eqtl = snps_eqtl.drop(columns=["variant_id"])


    # Check case-sensitivity and whitespace
    snps_pleio["chr"] = snps_pleio["chr"].astype(str).str.strip().str.lower()
    snps_pleio["pos"] = snps_pleio["pos"].astype(str).str.strip()

    susie_eqtl["chr"] = susie_eqtl["chr"].astype(str).str.strip().str.lower()
    susie_eqtl["pos"] = susie_eqtl["pos"].astype(str).str.strip()

    snps_eqtl["chr"] = snps_eqtl["chr"].astype(str).str.strip().str.lower()
    snps_eqtl["pos"] = snps_eqtl["pos"].astype(str).str.strip()
    snps_eqtl = snps_eqtl[["phenotype_id", "chr", "pos", "start_distance", 
                           "slope", "slope_se", "pval_nominal"]]


    # Overlapping causal variants
    snps_overlapped = pd.merge(susie_eqtl, snps_pleio, on=["chr", "pos"], how="inner")
    snps_overlapped = pd.merge(snps_overlapped, snps_eqtl, on=["phenotype_id", "chr", "pos"], how="inner")


    # Rename the first column to "gene_id"
    snps_overlapped.rename(columns={snps_overlapped.columns[0]: "gene_id"}, inplace=True)


    # Obtain information from gene data
    genes_eqtl_filt = genes_eqtl[["gene_id", "gene_start", "gene_end", "strand", 
                                  "beta_shape1", "beta_shape2", "pval_beta","qval"]]
    
    snps_overlapped = pd.merge(snps_overlapped, genes_eqtl_filt, on="gene_id", how="inner")


    # Filter based on Gene qvals
    clean_df = snps_overlapped[snps_overlapped["qval"] <= 0.05]
    clean_df.insert(0, 'tissue_name', tissue_name)

    clean_df = clean_df[[
        "tissue_name", "gene_id", "gene_name", "gene_start", "gene_end", "biotype", "strand", 
        "beta_shape1", "beta_shape2", "pval_beta", "qval", "variant_id_b37", "chr", "pos", "ref", "alt", 
        "build", "start_distance", "af", "cs_id", "cs_size", "pip", "slope", "slope_se", "afc", "afc_se", "pval_nominal"
    ]]

    return clean_df


if __name__ == "__main__":

    output_file = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/merged_gtex_susie_snps.tsv"

    gtex_susie_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_v11_Susie"
    gtex_snps_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_Analysis_v11_eQTL"
    snps_pleio = pd.read_csv("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/b37_to_b38/pleio_snps_b38.bed", sep="\t", header=None)
    

    # Get tissue names from the GTEx SuSiE directory
    gtex_susie_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_v11_Susie"
    suffix = ".v11.eQTLs.SuSiE_summary.parquet"
    files = [f for f in os.listdir(gtex_susie_dir) if f.endswith(suffix)]
    tissue_names = [f.replace(suffix, "") for f in files]


    result_dfs = pd.DataFrame()

    for tissue_name in tissue_names:

        print(f"Processing {tissue_name}...")

        # Construct file paths for the current tissue
        snps_eqtl_path = os.path.join(gtex_snps_dir, f"{tissue_name}.v11.eQTLs.signif_pairs.parquet")
        gene_eqtl_path = os.path.join(gtex_snps_dir, f"{tissue_name}.v11.eGenes.txt.gz")
        susie_eqtl_path = os.path.join(gtex_susie_dir, f"{tissue_name}.v11.eQTLs.SuSiE_summary.parquet")
       
        # Filter and merge data for the current tissue
        result_df = filter_gtex(tissue_name, snps_eqtl_path, gene_eqtl_path, susie_eqtl_path, snps_pleio)
        result_dfs = pd.concat([result_dfs, result_df], ignore_index=True)


    result_dfs.to_csv(output_file, sep="\t", index=False)