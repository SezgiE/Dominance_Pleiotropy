import numpy as np
import pandas as pd
import gseapy as gp
from gseapy import Biomart
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
bm = Biomart()


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


def process_data_enrich(background_genes_path, magma_res_path):
    
    # Whole genes in chr1-22 based on Ensembl b38
    bg_genes = pd.read_csv(background_genes_path, sep='\t')

    # Mapped genes from MAGMA
    mapped_genes = set()
    with open(magma_res_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) > 2:
                mapped_genes.add(parts[0])

    bg_genes['is_mapped'] = bg_genes['ensembl_gene_id'].isin(mapped_genes)

    # Loop through chromosomes to prepare input dataframe
    chrom_data = []
    for chrom in bg_genes['chromosome_name'].unique():
        subset = bg_genes[bg_genes['chromosome_name'] == chrom]
        chrom_data.append({
            "group_name": chrom,
            "mapped_in_group": subset['is_mapped'].sum(),
            "unmapped_in_group": (~subset['is_mapped']).sum()
        })

    chr_input = pd.DataFrame(chrom_data)


    biotype_data = []
    for biotype in bg_genes['gene_biotype'].unique():
        subset = bg_genes[bg_genes['gene_biotype'] == biotype]
        biotype_data.append({
            "group_name": biotype,
            "mapped_in_group": subset['is_mapped'].sum(),
            "unmapped_in_group": (~subset['is_mapped']).sum()
        })

    biotype_input = pd.DataFrame(biotype_data)


    return chr_input, biotype_input


def enrichment_getsets(magma_res_path, output_dir):
    
    # Mapped genes from MAGMA
    mapped_genes = set()
    with open(magma_res_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) > 2:
                mapped_genes.add(parts[0])
    
    # ENSG to gene symbol Map
    gene_id_dict ={'ensembl_gene_id': mapped_genes}
    results = bm.query(dataset='hsapiens_gene_ensembl',
                    attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'go_id'],
                    filters=gene_id_dict)
    gene_symbols = results["external_gene_name"].dropna().unique().tolist()

    # Run enrichment
    enr = gp.enrichr(
        gene_list=gene_symbols,
        gene_sets=['MSigDB_Hallmark_2020','KEGG_2026', 
                   "Reactome_Pathways_2024", "GO_Biological_Process_2025",
                   "GO_Molecular_Function_2025", "GO_Cellular_Component_2025",
                   "GTEx_Tissues_V8_2023"],
        organism='human',
        outdir=None
    )
    
    # Extract results and filter for significance (FDR < 0.05)
    results = enr.results
    sig_results = results[results['Adjusted P-value'] < 0.05].copy()
    
    # Calculate -log10(FDR) for plotting
    sig_results['log10(FDR)'] = -np.log10(sig_results['Adjusted P-value'])
    
    # Sort by significance
    sig_results = sig_results.sort_values('log10(FDR)', ascending=True)
    sig_results[["GoIs", "Total"]] = sig_results["Overlap"].str.split("/", n=2, expand=True).iloc[:, :2].astype(int)
    sig_results["Overlap"] = sig_results["GoIs"] / sig_results["Total"]

    sig_results.to_csv(f'{output_dir}/magma_gene_enrich_res.tsv', index=False, header = True, sep="\t")
    
    return sig_results


if __name__ == "__main__":

    background_genes_path = '/Users/sezgi/Documents/dominance_pleiotropy/gene_level/magma/magma_v1/ensembl_genes_raw38.tsv'
    magma_res_path = '/Users/sezgi/Documents/dominance_pleiotropy/gene_level/magma/magma_v1/magma_pleio_mapping.genes.annot'

    output_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/magma"

    
    fisher_input_chr, fisher_input_bio = process_data_enrich(background_genes_path, magma_res_path)
    
    # CHR enrichment
    chr_enrich_res = fisher_test(fisher_input_chr)
    chr_enrich_res = fdr_correction(chr_enrich_res)

    # Biotype enrichment
    bio_enrich_res = fisher_test(fisher_input_bio)
    bio_enrich_res = fdr_correction(bio_enrich_res)

    # Merge CHR and biotype Enrichment results
    chr_enrich_res.insert(0, 'type', "Chromosome")
    bio_enrich_res.insert(0, 'type', "Gene biotype")

    summary_df = pd.concat([chr_enrich_res, bio_enrich_res], ignore_index=True)
    summary_df.to_csv(f"{output_dir}/magma_enrich_summary.tsv", sep="\t", index=False)

    # Geneset enrichment
    geneset_enrich_res = enrichment_getsets(magma_res_path, output_dir)