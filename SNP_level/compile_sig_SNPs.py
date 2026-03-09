import pandas as pd
import collections
import json
import glob
import os


def compile_significant_snps(input_dir, variant_info, output_dir, hwe_sig=1e-6, info_threshold=0.9):

    # 1. Find all files in the input directory
    all_files = glob.glob(os.path.join(input_dir, "*.tsv.bgz"))
    all_files.sort()

    dfs = []

    var_dict = collections.defaultdict(lambda: {
        "sig_add_traits": [],
        "sig_dom_traits": [],
        "add_betas": {},
        "dom_betas": {},
        "add_pvals": {},
        "dom_pvals": {}
    })

    print(f"Found {len(all_files)} files. Extracting 'variant' and the last column...")

    # 2. Read each file, extract 'variant' and the last column (0 if not significant, 1 if only additive, 2 if both dominance & additive)  
    for f in all_files:

        filename = os.path.basename(f)
        phen_code = filename.replace("_sig_SNPs.tsv.bgz", "")

        try:
            
            # Read data
            df_raw = pd.read_csv(f, sep='\t', compression='gzip', 
                             usecols=['variant', "low_confidence_variant", phen_code, "pval", 
                                      "beta", "dominance_beta",  "dominance_pval"], 
                             dtype={'variant': str, phen_code: int, "pval": float, "beta": float, 
                                    "dominance_beta": float,  "dominance_pval": float})
            
            df = df_raw[df_raw["low_confidence_variant"].astype(str).str.lower() == 'false'].copy()
            
            removed_var = len(df_raw)-len(df)

            if removed_var > 0:
                print(f"{phen_code}: Removed {removed_var} low-confidence variants.")

            df.drop_duplicates(subset=['variant'], keep='first', inplace=True)
            

            records = df.to_dict('records')
            
            for row in records:
                vid = row['variant']
                status = row[phen_code]
                
                # Only store metrics if the variant is significant in this trait (1 or 2)
                if status > 0:
                    var_dict[vid]["sig_add_traits"].append(phen_code)
                    var_dict[vid]["add_pvals"][phen_code] = row['pval']
                    var_dict[vid]["add_betas"][phen_code] = row['beta']
                    
                    if status == 2: # Both dominance & additive
                        var_dict[vid]["sig_dom_traits"].append(phen_code)
                        var_dict[vid]["dom_pvals"][phen_code] = row['dominance_pval']
                        var_dict[vid]["dom_betas"][phen_code] = row['dominance_beta']

            
            # Set variant as the index for fast multi-file merging
            df = df[['variant', phen_code]].copy()
            df.set_index('variant', inplace=True)
            dfs.append(df)
            
        except Exception as e:
            print(f"Skipping {f} - error reading file: {e}")

    print(f"Successfully loaded {len(dfs)} files. Building the final matrix...")


    # 3. Concatenate everything side-by-side based on the 'variant' index in a significance matrix format
    merged_matrix = pd.concat(dfs, axis=1)
    merged_matrix.fillna(0, inplace=True)

    # Bring 'variant' back out of the index into a normal column
    merged_matrix.reset_index(inplace=True)
    

    # 4. Calculate total significant counts for additive and dominance across all phenotypes
    pheno_cols = [c for c in merged_matrix.columns if c != 'variant']
    merged_matrix['add_sig_total'] = merged_matrix[pheno_cols].isin([1, 2]).sum(axis=1)
    merged_matrix['dom_sig_total'] = (merged_matrix[pheno_cols] == 2).sum(axis=1)
    print(f"Final matrix shape after merging: {merged_matrix.shape}.")

    # Convert the dictionary to a DataFrame and append to merged_matrix
    print("Appending dictionary columns (JSONs and trait lists)...")
    dict_rows = []
    for vid, data in var_dict.items():
        dict_rows.append({
            "variant": vid,
            "sig_add_traits": ", ".join(data["sig_add_traits"]),
            "sig_dom_traits": ", ".join(data["sig_dom_traits"]),
            "add_betas": json.dumps(data["add_betas"]),
            "dom_betas": json.dumps(data["dom_betas"]),
            "add_pvals": json.dumps(data["add_pvals"]),
            "dom_pvals": json.dumps(data["dom_pvals"])
        })

    dict_df = pd.DataFrame(dict_rows)
    
    # Merge the new columns onto your existing wide matrix
    merged_matrix = merged_matrix.merge(dict_df, on='variant', how='left')
    print(merged_matrix.head())

    # 5. Read the variant info to get chr, rsID, info, and minor_AF for each variant
    print("Reading variant information for filtering...")
    var_info = pd.read_csv(variant_info, sep='\t', compression='gzip',
                            usecols=["variant", "chr", "rsid", "info", "minor_AF", "p_hwe"],
                            dtype={"variant": str, "chr": str, "rsid": str, "info": float, "minor_AF": float, "p_hwe": float})
    
    print(f"Total variants in variant information file: {len(var_info)}.")
    var_info = var_info[var_info['variant'].isin(merged_matrix['variant'])]
    print(f"Variants in the significance matrix that have variant info: {len(var_info)}.")


    # 6. Apply HWE, INFO, Chromosome filters and remove INDELs and multi-allelic variants
    print("Applying filters: HWE p-value > 1e-6, info > 0.9, and removing sex chromosomes...")
    var_info_HWE = var_info[var_info['p_hwe'] >= hwe_sig].copy()
    print(f"removed {len(var_info) - len(var_info_HWE)} variants based on HWE p-value threshold of {hwe_sig}. Remaining: {len(var_info_HWE)}.")
    var_info_HEW_INFO = var_info_HWE[var_info_HWE['info'] > info_threshold].copy()
    print(f"removed {len(var_info_HWE) - len(var_info_HEW_INFO)} variants based on INFO score threshold of {info_threshold}. Remaining: {len(var_info_HEW_INFO)}.")

    # Removing X (and adding Y, XY, and MT just to be perfectly safe for strict autosome-only analysis)
    var_info_HWE_INFO_sex = var_info_HEW_INFO[~var_info_HEW_INFO['chr'].isin(['X', 'Y', 'XY', 'MT'])].copy()
    print(f"removed {len(var_info_HEW_INFO) - len(var_info_HWE_INFO_sex)} variants based on chromosome filter.")

    # Filter out INDELS
    var_filtered_IND = var_info_HWE_INFO_sex[var_info_HWE_INFO_sex['variant'].apply(lambda x: 
                                                                                len(x.split(':')[2]) == 1 and len(x.split(':')[3]) == 1)].copy()
    print(f"removed {len(var_info_HWE_INFO_sex) - len(var_filtered_IND)} variants based on INDEL filter. Remaining: {len(var_filtered_IND)}.")

    # Filter out diallelic variants (those with more than one alternate allele)
    var_filtered_IND['pos'] = var_filtered_IND['variant'].apply(lambda x: int(x.split(':')[1]))
    var_filtered_IND_DIAL = var_filtered_IND.drop_duplicates(subset=['chr', 'pos'], keep=False).copy()
    print(f"removed {len(var_filtered_IND) - len(var_filtered_IND_DIAL)} variants based on diallelic filter. Remaining: {len(var_filtered_IND_DIAL)}.")


    # 7. Merge the filtered variant info back to the  significance matrix to keep only high-quality variants
    print(f"After filtering steps, {len(var_filtered_IND_DIAL)} variants remain.")
    print("Merging variant annotations back to the significance matrix...")
    merged_matrix = merged_matrix.merge(var_filtered_IND_DIAL, on='variant', how='inner')


    # 8. Save the final output
    output_all = f"{output_dir}/all_sig_SNPs.tsv.gz"
    print(f"Saving {len(merged_matrix)} high-quality unique SNPs to {output_all}...")
    merged_matrix.to_csv(output_all, sep='\t', index=False, compression='gzip')
    print("Done!")


if __name__ == "__main__":

# Path setup
    input_dir = "/Users/sezgi/Documents/overlapped_SNPs/significant_SNPs/sig_sumstats_merged"
    variant_info = "/Users/sezgi/Documents/overlapped_SNPs/UKB_sumstats_Neale/variants.tsv.bgz"
    output_dir = "/Users/sezgi/Documents/overlapped_SNPs/significant_SNPs"

    compile_significant_snps(input_dir, variant_info, output_dir)