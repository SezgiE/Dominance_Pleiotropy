import glob
import pandas as pd

info_df = pd.concat((pd.read_csv(f, sep='\t') for f in glob.glob("/Users/sezgi/Documents/dominance_pleiotropy/infos/*.info")), ignore_index=True)

sig_var = pd.read_csv("/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/significant_SNPs/all_sig_SNPs.tsv.gz", 
                       sep='\t', compression='gzip', usecols=["variant", "chr", "pos", "rsid", "dom_sig_total"],
                        dtype={"variant": str, "chr": str, "rsid": str, "dom_sig_total":int})
sig_var["key"] = sig_var['chr'].astype(str) + ":" + sig_var['pos'].astype(str)

var_info = pd.read_csv("/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/variants.tsv.bgz", 
                       sep='\t', compression='gzip', usecols=["variant", "chr", "rsid", "info", "minor_AF", "p_hwe"],
                        dtype={"variant": str, "chr": str, "rsid": str, "info": float, "minor_AF": float, "p_hwe": float})


var_info = var_info[var_info['minor_AF'] >= 1e-6].copy()
var_info_HWE = var_info[var_info['p_hwe'] >= 1e-6].copy()
var_info_HEW_INFO = var_info_HWE[var_info_HWE['info'] > 0.9].copy()
var_info_HWE_INFO_sex = var_info_HEW_INFO[~var_info_HEW_INFO['chr'].isin(['X', 'Y', 'XY', 'MT'])].copy()
var_filtered_IND = var_info_HWE_INFO_sex[var_info_HWE_INFO_sex['variant'].apply(lambda x: 
                                                                            len(x.split(':')[2]) == 1 and len(x.split(':')[3]) == 1)].copy()

var_filtered_IND['pos'] = var_filtered_IND['variant'].apply(lambda x: int(x.split(':')[1]))
var_filtered_IND_DIAL = var_filtered_IND.drop_duplicates(subset=['chr', 'pos'], keep=False).copy()

var_filtered_IND_DIAL["key"] = var_filtered_IND_DIAL['chr'].astype(str) + ":" + var_filtered_IND_DIAL['pos'].astype(str)
info_df["key"] = info_df["CHR"].astype(str) + ":" + info_df["POS"].astype(str)


merged_matrix = var_filtered_IND_DIAL.merge(info_df, left_on='rsid', right_on='SNP', how="inner")
merged_matrix_key = var_filtered_IND_DIAL.merge(info_df, on="key", how="inner")

print(len(info_df))
print(len(merged_matrix))
print(len(merged_matrix_key))


sig_var_dom = sig_var[sig_var["dom_sig_total"] > 0]
sig_var_merged = sig_var_dom.merge(info_df, left_on='rsid', right_on='SNP', how="inner")
sig_var_merged_key = sig_var_dom.merge(info_df, on="key", how="inner")

print(len(sig_var))
print(len(sig_var_dom))
print(len(sig_var_merged))
print(len(sig_var_merged_key))