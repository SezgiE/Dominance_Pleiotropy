import pandas as pd
import glob
import os

def sig_phen_summary(trait_df, phen_info_df, add_pleiotropic_snps, dom_pleiotropic_snps, trait_code):
    """
    Summarizes significant SNPs per chromosome for a given trait.
    """
    # Filter phen_info for the specific trait
    phen_info = phen_info_df[phen_info_df['phenotype_code'] == trait_code]
    
    # Exclude mhc region (chr6: 25Mb-34Mb)
    mhc_region = (trait_df['chr'] == 6) & (trait_df['pos'] >= 25000000) & (trait_df['pos'] <= 34000000)
    trait_df['is_mhc'] = mhc_region
    
    # boolean conditions
    trait_df['is_add'] = trait_df['add_sig'] == 1
    trait_df['is_dom'] = trait_df['dom_sig'] == 2
    trait_df['is_both'] = trait_df[trait_code] == 3

    # Is pleiotropic for add/dom or both
    trait_df['is_add_pleiotropic'] = trait_df['is_add'] & trait_df['variant'].isin(add_pleiotropic_snps)
    trait_df['is_dom_pleiotropic'] = trait_df['is_dom'] & trait_df['variant'].isin(dom_pleiotropic_snps)
    trait_df['is_both_pleiotropic'] = trait_df['is_add_pleiotropic'] & trait_df['is_dom_pleiotropic']

    # Group by chr and sum the boolean columns
    summary = trait_df.groupby(['chr', 'is_mhc']).agg(
        add_sig=('is_add', 'sum'),
        dom_sig=('is_dom', 'sum'),
        both_sig=('is_both', 'sum'),
        add_pleio =('is_add_pleiotropic', 'sum'),
        dom_pleio =('is_dom_pleiotropic', 'sum'),
        both_pleio =('is_both_pleiotropic', 'sum')
    ).reset_index()
    
    # Add phenotype_code for the merge
    summary['phenotype_code'] = trait_code
    
    # Filter out chromosomes with zero add/dom signals
    summary = summary[(summary['dom_sig'] > 0)]
    
    # Merge with phen_info (similar to left_join/inner_join)
    summary = pd.merge(phen_info, summary, on='phenotype_code', how='inner')
    
    # Reorder columns: phen_info columns first, then chr, then the signal counts
    final_cols = list(phen_info.columns) + ['chr', 'is_mhc', 'add_sig', 'dom_sig', 'both_sig', 'add_pleio', 'dom_pleio', 'both_pleio']
    
    return summary[final_cols]


if __name__ == "__main__":
    
    # Directories
    trait_df_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sumstats_QCed"
    phen_info_path = "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/supp_table1.tsv"
    all_sig_snps_path = "/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/significant_SNPs/all_sig_SNPs.tsv.gz"
    output_path = "/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/results/phen_summary.tsv"

    # Load phenotype metadata
    cols_to_use = ["phenotype_code", "phenotype_description", "category", "variable_type", 
                   "source", "n_non_missing", "n_missing", "n_controls", "n_cases"]
    
    phen_info_df = pd.read_csv(phen_info_path, sep='\t', usecols=cols_to_use)

    # Load all significant SNPs (if needed for further filtering, not currently used in the summary)
    all_sig_snps_df = pd.read_csv(all_sig_snps_path, sep='\t', compression='gzip', 
                                  usecols=["variant", "add_sig_total", "dom_sig_total"])
    
    add_pleiotropic_snps = set(all_sig_snps_df[all_sig_snps_df['add_sig_total'] > 1]['variant'])
    dom_pleiotropic_snps = set(all_sig_snps_df[all_sig_snps_df['dom_sig_total'] > 1]['variant'])
    
    # Find all "_sig_SNPs.tsv.bgz" files
    file_pattern = os.path.join(trait_df_dir, "*_sig_SNPs.tsv.bgz")
    trait_files = glob.glob(file_pattern)
    
    print(f"Found {len(trait_files)} files to process.")
    
    phen_summaries = []
    
    for file_path in trait_files:
        # Extract trait_code from filename
        filename = os.path.basename(file_path)
        trait_code = filename.replace("_sig_SNPs.tsv.bgz", "")
        
        # Read the trait dataframe
        trait_df = pd.read_csv(file_path, sep='\t', compression='gzip')
        
        # Get summary
        phen_summary = sig_phen_summary(trait_df, phen_info_df, add_pleiotropic_snps, dom_pleiotropic_snps, trait_code)
        phen_summaries.append(phen_summary)
        
    # Combine all summaries and save
    if phen_summaries:

        final_df = pd.concat(phen_summaries, ignore_index=True)
        # this should return True if all rows for a given phenotype_code have is_mhc = True, otherwise False
        final_df['MHC_filtered'] = final_df.groupby('phenotype_code')['is_mhc'].transform('all')
        
        final_df = final_df.sort_values(by=['MHC_filtered', 'phenotype_code', 'chr'],
                                        ascending=[True, True, True]).reset_index(drop=True)

        final_df.insert(0, 'trait_number', pd.factorize(final_df['phenotype_code'])[0] + 1)
        
        final_df.to_csv(output_path, sep='\t', index=False)
        print(f"Saved final summary to: {output_path}")
    else:
        print("No summaries were generated.")