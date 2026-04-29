import os
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

    cols_to_fill = ["pleio_signal_itt", "non_pleio_signal_itt"]
    df[cols_to_fill] = df[cols_to_fill].fillna(0)

    total_pleio = df["pleio_signal_itt"].sum()
    total_non_pleio = df["non_pleio_signal_itt"].sum()

    # Pleio signal but not in this tissue (ntt)
    df["pleio_signal_ntt"] = total_pleio - df["pleio_signal_itt"]

    # Non-Pleio and not in this tissue (ntt)              
    df["non_pleio_signal_ntt"] = total_non_pleio - df["non_pleio_signal_itt"] 
    

    def apply_fisher(row):
        table = [
            [row["pleio_signal_itt"], row["non_pleio_signal_itt"]],
            [row["pleio_signal_ntt"], row["non_pleio_signal_ntt"]]
        ]
        return pd.Series(fisher_exact(table, alternative='greater'))
    
    df[["odds_ratio", "p_value"]] = df.apply(apply_fisher, axis=1)
    

    return df
   

def get_total_snps_and_pips(gtex_susie_dir, files, tissue_names, snps_pleio, output_dir):

    print("Obtaning total PIPs...")
    cols_to_read = ["phenotype_id", "cs_id", "variant_id", "pip"]  
    df_list = [
        pd.read_parquet(os.path.join(gtex_susie_dir, f), columns=cols_to_read).assign(tissue=t) 
        for f, t in zip(files, tissue_names)
    ]
    
    combined_df = pd.concat(df_list, ignore_index=True)


    combined_df[["ENSG", "version"]] = combined_df["phenotype_id"].str.split(".", n=2, expand=True).iloc[:, :2]
    combined_df[["chr", "pos"]] = combined_df["variant_id"].str.split("_", n=2, expand=True).iloc[:, :2]
    
    combined_df["var_id"] = combined_df["chr"].astype(str) + ":" + combined_df["pos"].astype(str)
    snps_pleio["var_id"] = snps_pleio["chr"].astype(str) + ":" + snps_pleio["pos"].astype(str)

    
    print("Total independent causal signal: ", combined_df.groupby(["tissue", 'ENSG'])['cs_id'].nunique().sum())
    
    
    # Fuma Gene2func background genes
    background_genes = combined_df.drop_duplicates(subset=["ENSG"])
    print("Total Unique Causal Genes: ", len(background_genes))
    background_genes = background_genes[["ENSG"]]  
    background_genes.to_csv(f'{output_dir}/fuma_background_genes.txt', index=False)
    

    # Combine all the PIP values
    all_pip_vals = combined_df["pip"].dropna().tolist()
    pleio_pip_vals= pd.merge(combined_df, snps_pleio[["var_id"]].drop_duplicates(), on="var_id", how="inner")["pip"].dropna().tolist()

    df_all = pd.DataFrame({"PIP": all_pip_vals,
                           "source": "all"
                           })
    
    df_pleio = pd.DataFrame({"PIP": pleio_pip_vals,
                             "source": "pleio"
                             })

    pip_df = pd.concat([df_all, df_pleio], ignore_index=True)

    pip_df["source"] = pip_df["source"].astype("category")
    pip_df["PIP"] = pip_df["PIP"].astype("float32")
    

    return pip_df.to_parquet(f"{output_dir}/merged_pip_values.parquet", index=False)


def gtex_analyze_tissues(tissue_name, snps_eqtl_path, gene_eqtl_path, 
                susie_eqtl_path, snps_pleio):

    snps_eqtl = pd.read_parquet(snps_eqtl_path)
    genes_eqtl = pd.read_csv(gene_eqtl_path, sep="\t", compression="gzip")
    susie_eqtl = pd.read_parquet(susie_eqtl_path)


    # Adjust SNPs and Susie eqtl
    susie_eqtl[["ENSG", "version"]] = susie_eqtl["phenotype_id"].str.split(".", n=2, expand=True).iloc[:, :2]
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
    snps_overlapped = pd.merge(snps_overlapped, snps_eqtl, on=["phenotype_id", "chr", "pos"], how="left")


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
        "tissue_name", "gene_id", "ENSG", "version", "gene_name", "gene_start", "gene_end", "biotype", "strand", 
        "beta_shape1", "beta_shape2", "pval_beta", "qval", "variant_id_b37", "chr", "pos", "ref", "alt", 
        "build", "start_distance", "af", "cs_id", "cs_size", "pip", "slope", "slope_se", "afc", "afc_se", "pval_nominal"
    ]]


    # Summarize tissue for Fisher's test
    pleio_tissue = clean_df.groupby('ENSG')['cs_id'].nunique().sum()

    tissue_all = susie_eqtl.groupby('ENSG')['cs_id'].nunique().sum()
    
    # pleio_tissue = (clean_df["chr"].astype(str) + ":" + clean_df["pos"].astype(str)).tolist(
    # tissue_all = (susie_eqtl["chr"].astype(str) + ":" + susie_eqtl["pos"].astype(str)).tolist()

    A = pleio_tissue                      # Pleio signal in this tissue (itt)
    B = tissue_all  - pleio_tissue        # Non-Pleio signal in this tissue (itt)
        
    tissue_summary = {
        "name": tissue_name,
        "pleio_signal_itt": A,
        "non_pleio_signal_itt": B
    }

    # Summarize Gene biotypes for Fisher's test
    pleio_biotype = clean_df.groupby(["biotype", "ENSG"])["cs_id"].nunique().groupby("biotype").sum()
    tissue_all_biotype = susie_eqtl.groupby(["biotype", "ENSG"])["cs_id"].nunique().groupby("biotype").sum()

    biotype_df = pd.DataFrame({
        "A_b": pleio_biotype, 
        "tissue_total": tissue_all_biotype
    }).fillna(0)
    
    biotype_df["B_b"] = biotype_df["tissue_total"] - biotype_df["A_b"]
    
    # list of dictionaries for each biotype
    biotype_summary_list = []
    for biotype, row in biotype_df.iterrows():
        biotype_summary_list.append({
            "name": tissue_name,
            "biotype": biotype,
            "pleio_signal_itt": int(row["A_b"]),
            "non_pleio_signal_itt": int(row["B_b"])
        })
    
    biotype_summary_df = pd.DataFrame(biotype_summary_list) if biotype_summary_list else pd.DataFrame()

    return clean_df, tissue_summary, biotype_summary_df


if __name__ == "__main__":

    # Define output directory
    output_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/gtex_res"

    # GTEx files directories
    gtex_susie_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_v11_Susie"
    gtex_snps_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_Analysis_v11_eQTL"
    
    # Read pleiotropic SNPs
    snps_pleio = pd.read_csv("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/b37_to_b38/pleio_snps_b38.bed", sep="\t", header=None)
    snps_pleio.columns = ["chr", "start", "pos", "variant_id_b37"]
    snps_pleio = snps_pleio[["chr", "pos", "variant_id_b37"]]


    #------------------------------------ Analysis by Tissue Type --------------------------------#
    # Get tissue names from the GTEx SuSiE directory
    gtex_susie_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_v11_Susie"
    suffix = ".v11.eQTLs.SuSiE_summary.parquet"
    files = [f for f in os.listdir(gtex_susie_dir) if f.endswith(suffix)]
    tissue_names = [f.replace(suffix, "") for f in files]

    get_total_snps_and_pips(gtex_susie_dir, files, tissue_names, snps_pleio, output_dir)

    result_list = []
    tissue_summaries = []
    biotype_summaries = []

    for tissue_name in tissue_names:

        print(f"Processing {tissue_name}...")

        # Construct file paths for the current tissue
        snps_eqtl_path = os.path.join(gtex_snps_dir, f"{tissue_name}.v11.eQTLs.signif_pairs.parquet")
        gene_eqtl_path = os.path.join(gtex_snps_dir, f"{tissue_name}.v11.eGenes.txt.gz")
        susie_eqtl_path = os.path.join(gtex_susie_dir, f"{tissue_name}.v11.eQTLs.SuSiE_summary.parquet")
       
        # Filter and merge data for the current tissue
        result_df, tissue_summary, biotype_summary = gtex_analyze_tissues(tissue_name, snps_eqtl_path, gene_eqtl_path, 
                                            susie_eqtl_path, snps_pleio)
        result_list.append(result_df)
        tissue_summaries.append(tissue_summary)
        biotype_summaries.append(biotype_summary)


    result_df = pd.concat(result_list, ignore_index=True) if result_list else pd.DataFrame()
    result_df.to_csv(f"{output_dir}/gtex_susie_pleio_snps.tsv", sep="\t", index=False)

    # Fuma Gene2func pleiotropic genes (gene of interest)
    pleio_genes = result_df[['ENSG']].drop_duplicates()
    print("Total Pleiotropic Causal Genes: ", len(pleio_genes))    
    pleio_genes.to_csv(f'{output_dir}/fuma_pleio_genes.txt', index=False)

    # Tissue Enrichment Test
    tissue_summary_df = pd.DataFrame(tissue_summaries) if tissue_summaries else pd.DataFrame()
       
    tissue_summary_df = fisher_test(tissue_summary_df)
    tissue_summary_df = fdr_correction(tissue_summary_df)

    # Biotype Enrichment Test
    merged_biotype_summary_df = pd.concat(biotype_summaries, ignore_index=True) if biotype_summaries else pd.DataFrame()
    biotype_results = merged_biotype_summary_df.groupby(["biotype"], as_index=False)[["pleio_signal_itt", "non_pleio_signal_itt"]].sum()
    biotype_results.rename(columns={"biotype": "name"}, inplace=True)

    biotype_results = fisher_test(biotype_results)
    biotype_results = fdr_correction(biotype_results)

    # Merge Enrichment results
    tissue_summary_df.insert(0, 'type', "tissue")
    biotype_results.insert(0, 'type', "biotype")

    eqtl_summary_df = pd.concat([tissue_summary_df, biotype_results], ignore_index=True)
    eqtl_summary_df.to_csv(f"{output_dir}/eqtl_summary.tsv", sep="\t", index=False)

