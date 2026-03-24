import os
import sys
import numpy as np
import pandas as pd


def get_data(phen_code, sumstat_dir, loci_dir):
    """Gets and merges locus data."""

    raw_df = pd.read_csv(f"{sumstat_dir}/{phen_code}_sig_SNPs.tsv.bgz", sep='\t', compression='gzip')
    loci_df = pd.read_csv(f"{loci_dir}/{phen_code}_sig_loci.tsv", sep='\t', 
                          usecols=["variant", "indep_status", "indep_id", "r4", 
                                   "lead_status", "lead_id","r4_lead",  "indep_start", 
                                   "indep_end", "ld_start", "ld_end", "ld_id"])
    
    if loci_df.empty:
        print("No loci to plot.")
        sys.exit(0)
    
    merged_df = raw_df.merge(loci_df, on="variant", how="left")

    print(merged_df[merged_df["dominance_pval"] == 0]["dominance_pval"].head())
    return merged_df


def plot_regional_association(merged_df, genes_df, phen_code, phen_name, output_dir, buffer_kb=50):
    """Generates LocusZoom-style plots for each unique independent locus."""
    
    os.makedirs(output_dir, exist_ok=True)

    # Identify the unique ld blocks
    unique_blocks = merged_df['ld_id'].dropna().unique()
    n_blocks = len(unique_blocks)
    print(f"Found {len(unique_blocks)} unique loci to run SuSiE for.")

    merged_df["-log10P"] = -np.log10(merged_df['dominance_pval'])


    for idx, block in enumerate(unique_blocks):
        print(f"Plotting locus: {block}")
        
        locus_data = merged_df[merged_df['ld_id'] == block].copy()
        
        chrom = locus_data['chr'].iloc[0]
        ld_start = locus_data['ld_start'].min()
        ld_end = locus_data['ld_end'].max()
        
        # Define the plot window (add the buffer)
        plot_start = ld_start - (buffer_kb * 1000)
        plot_end = ld_end + (buffer_kb * 1000)
        
        # Slice the specific window for plotting
        window_df = merged_df[
            (merged_df['chr'] == chrom) & 
            (merged_df['pos'] >= plot_start) & 
            (merged_df['pos'] <= plot_end)
        ].copy()
        

        # SNPs not in locus_data will get NaN, which we fill with -1 (for background grey)
        window_df['r4'] = window_df['r4'].fillna(-1)

        window_genes = genes_df[
            (genes_df['chr'] == chrom) & 
            (genes_df['end'] >= plot_start) & 
            (genes_df['start'] <= plot_end)
        ].copy()

        

if __name__ == "__main__":
    
    sumstat_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sumstats_QCed"
    phen_dict_path = "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/phen_dict_renamed.xlsx"
    loci_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sig_loci"
    out_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/susie_results"

    traits = pd.read_excel(phen_dict_path, usecols=["phenotype_code", "description"])



    for index, row in traits.iterrows():
        
        phen_code = row['phenotype_code']
        phen_name = row['description']

        try:
            merged_df = get_data(phen_code, sumstat_dir, loci_dir)
        except FileNotFoundError:
            # Skip to the next phenotype if the file/dir doesn't exist
            print(f"No significant loci data for{phen_code} ")
            continue
        
        print(f"Processing: {phen_code}")
        #plot_regional_association(merged_df, genes_df, phen_code, phen_name, out_dir)