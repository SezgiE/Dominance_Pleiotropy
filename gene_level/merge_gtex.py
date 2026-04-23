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


def get_total_snps_and_pips(gtex_susie_dir, files, tissue_names, snps_pleio, output_dir):

    print("Obtaning total PIPs...")
    cols_to_read = ["variant_id", "biotype", "pip"]  # Specify your required columns
    df_list = [pd.read_parquet(os.path.join(gtex_susie_dir, f), columns=cols_to_read).assign(tissue=t) for f, t in zip(files, tissue_names)]
    combined_df = pd.concat(df_list, ignore_index=True)
    
    combined_df[["chr", "pos"]] = combined_df["variant_id"].str.split("_", n=2, expand=True).iloc[:, :2]
    combined_df["var_id"] = combined_df[["chr", "pos"]].astype(str).agg(':'.join, axis=1)
    snps_pleio["var_id"] = snps_pleio[["chr", "pos"]].astype(str).agg(':'.join, axis=1)

    set1 = set(combined_df["var_id"])
    set2 = set(snps_pleio["var_id"])
    total_snps = set1.union(set2)

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
    
    pip_df.to_parquet(f"{output_dir}/merged_pip_values.parquet", index=False)

    gtex_susie_snps_biotype = combined_df[["var_id", "biotype"]]

    return total_snps, gtex_susie_snps_biotype
    

def gtex_tissue(tissue_name, snps_eqtl_path, gene_eqtl_path, 
                susie_eqtl_path, snps_pleio, total_snps):

    snps_eqtl = pd.read_parquet(snps_eqtl_path)
    genes_eqtl = pd.read_csv(gene_eqtl_path, sep="\t", compression="gzip")
    susie_eqtl = pd.read_parquet(susie_eqtl_path)


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
        "tissue_name", "gene_id", "gene_name", "gene_start", "gene_end", "biotype", "strand", 
        "beta_shape1", "beta_shape2", "pval_beta", "qval", "variant_id_b37", "chr", "pos", "ref", "alt", 
        "build", "start_distance", "af", "cs_id", "cs_size", "pip", "slope", "slope_se", "afc", "afc_se", "pval_nominal"
    ]]


    # Fisher's Exact Test
    pleio_tissue = set(clean_df["chr"].astype(str) + ":" + clean_df["pos"].astype(str))
    pleio_all = set(snps_pleio["chr"].astype(str) + ":" + snps_pleio["pos"].astype(str))
    tissue_all = set(susie_eqtl["chr"].astype(str) + ":" + susie_eqtl["pos"].astype(str))

    A = len(pleio_tissue)                                      # Pleio and in this tissue
    B = len(tissue_all - pleio_tissue)                         # Non-Pleio in but in this tissue
    C = len(pleio_all - tissue_all)                            # Pleio in but not in this tissue
    D = len(total_snps.difference(pleio_all, tissue_all))     # Non-Pleio and not in this tissue 

    odds_ratio, p_val = fisher_exact([[A, B], [C, D]], alternative='greater')
        
    results_fisher = {
        "name": tissue_name,
        "pleio_variants": A,
        "total_sig_variants": A + B,
        "odds_ratio": odds_ratio,
        "p_value": p_val
    }


    return clean_df, results_fisher


def gtex_biotype(pleio_bio_results, pleio_all, gtex_susie_snps_biotype, total_snps):

    print("Running Fisher's Exact test for biotype enrichment")

    biotypes = pleio_bio_results["biotype"].unique().tolist()

    fisher_res_list = []
    
    for biotype in biotypes:

        print(f"Testing biotype: {biotype}")
        # Bio and pleio
        pleio_bio_df = pleio_bio_results[pleio_bio_results["biotype"] == biotype]
        pleio_bio = set(pleio_bio_df["chr"].astype(str) + ":" + pleio_bio_df["pos"].astype(str))
        
        # Pleio
        pleio_all = set(snps_pleio["chr"].astype(str) + ":" + snps_pleio["pos"].astype(str))
        
        # Bio
        susie_biotype = gtex_susie_snps_biotype[gtex_susie_snps_biotype["biotype"] == biotype]
        bio_all = set(susie_biotype["var_id"])
        
        # Fisher's Exact Test
        A = len(pleio_bio)                                         # Pleio and in this bio
        B = len(bio_all - pleio_bio)                               # Non-Pleio in but in this bio
        C = len(pleio_all - bio_all)                               # Pleio in but not in this bio
        D = len(total_snps.difference(pleio_all, bio_all))         # Non-Pleio and not in this tissue 

        odds_ratio, p_val = fisher_exact([[A, B], [C, D]], alternative='greater')
            
        results_fisher = {
            "name": biotype,
            "pleio_variants": A,
            "total_sig_variants": A + B,
            "odds_ratio": odds_ratio,
            "p_value": p_val
        }
        
        fisher_res_list.append(results_fisher)

    return fisher_res_list


if __name__ == "__main__":

    # Define output directory
    output_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex"

    # GTEx files directories
    gtex_susie_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_v11_Susie"
    gtex_snps_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_Analysis_v11_eQTL"
    
    # Read pleiotropic SNPs
    snps_pleio = pd.read_csv("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/b37_to_b38/pleio_snps_mult_b38.bed", sep="\t", header=None)
    snps_pleio.columns = ["chr", "start", "pos", "variant_id_b37"]
    snps_pleio = snps_pleio[["chr", "pos", "variant_id_b37"]]


    #------------------------------------ Analysis by Tissue Type --------------------------------#
    # Get tissue names from the GTEx SuSiE directory
    gtex_susie_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_v11_Susie"
    suffix = ".v11.eQTLs.SuSiE_summary.parquet"
    files = [f for f in os.listdir(gtex_susie_dir) if f.endswith(suffix)]
    tissue_names = [f.replace(suffix, "") for f in files]

    # Total unique SNPs across all tissues + pleiotropic SNPs for Fisher's test
    total_snps, gtex_susie_snps_biotype = get_total_snps_and_pips(gtex_susie_dir, files, tissue_names, snps_pleio, output_dir)

    result_list = []
    fisher_tissue_list = []

    for tissue_name in tissue_names:

        print(f"Processing {tissue_name}...")

        # Construct file paths for the current tissue
        snps_eqtl_path = os.path.join(gtex_snps_dir, f"{tissue_name}.v11.eQTLs.signif_pairs.parquet")
        gene_eqtl_path = os.path.join(gtex_snps_dir, f"{tissue_name}.v11.eGenes.txt.gz")
        susie_eqtl_path = os.path.join(gtex_susie_dir, f"{tissue_name}.v11.eQTLs.SuSiE_summary.parquet")
       
        # Filter and merge data for the current tissue
        result_df, fisher_tissue_res = gtex_tissue(tissue_name, snps_eqtl_path, gene_eqtl_path, 
                                            susie_eqtl_path, snps_pleio, total_snps)
        result_list.append(result_df)
        fisher_tissue_list.append(fisher_tissue_res)


    result_df = pd.concat(result_list, ignore_index=True) if result_list else pd.DataFrame()
    result_df.to_csv(f"{output_dir}/gtex_susie_pleio_snps.tsv", sep="\t", index=False)


    fisher_tissue_df = pd.DataFrame(fisher_tissue_list) if fisher_tissue_list else pd.DataFrame()
    fisher_tissue_df = fdr_correction(fisher_tissue_df)


    #------------------------------------ Analysis by Biotype --------------------------------#
    fisher_bio_res = gtex_biotype(result_df, snps_pleio, gtex_susie_snps_biotype, total_snps)
    fisher_bio_df = pd.DataFrame(fisher_bio_res) if fisher_bio_res else pd.DataFrame()

    fisher_bio_df = fdr_correction(fisher_bio_df)

    # Merged
    fisher_tissue_df.insert(0, 'type', "tissue")
    fisher_bio_df.insert(0, 'type', "biotype")

    fisher_all_df = pd.concat([fisher_tissue_df, fisher_bio_df], ignore_index=True)
    fisher_all_df.to_csv(f"{output_dir}/fisher_results.tsv", sep="\t", index=False)