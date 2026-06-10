import pandas as pd
import numpy as np
import json
from scipy import stats

def parse_betas(beta_str):
    if pd.isna(beta_str):
        return {}
    
    cleaned_str = str(beta_str).replace('""', '"')
    
    try:
        return json.loads(cleaned_str)
    except json.JSONDecodeError:
        return {}


def calculate_corr(data):
    # Correlation Tests
    n_obs = len(data)
    df = n_obs - 2

    # Calculate Pearson correlations
    add_res = stats.pearsonr(data['maf'], data['add_beta'])
    dom_res = stats.pearsonr(data['maf'], data['dom_beta'])

   
    add_ci = add_res.confidence_interval(confidence_level=0.95)
    dom_ci = dom_res.confidence_interval(confidence_level=0.95)

    
    def apa_fmt(val, decimals):
        formatted = f"{val:.{decimals}f}"
        if formatted.startswith("0."): return formatted[1:]
        if formatted.startswith("-0."): return "-" + formatted[2:]
        return formatted

    # Formatted prints
    print("Additive Betas vs MAF:")
    print(f"r({df}) = {apa_fmt(add_res.statistic, 2)}, p = {add_res.pvalue}, 95% CI [{apa_fmt(add_ci.low, 2)}, {apa_fmt(add_ci.high, 2)}]\n")

    print("Dominance Betas vs MAF:")
    print(f"r({df}) = {apa_fmt(dom_res.statistic, 2)}, p = {dom_res.pvalue}, 95% CI [{apa_fmt(dom_ci.low, 2)}, {apa_fmt(dom_ci.high, 2)}]")

# For all SNPs
# Load the data
file_path = "/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/significant_SNPs/all_sig_SNPs.tsv.gz"
cols_to_use = ['variant', 'minor_AF', 'add_betas', 'dom_betas']

df = pd.read_csv(file_path, sep='\t', usecols=cols_to_use)

# Expand to long format and align traits
records = []
for row in df.itertuples(index=False):
    variant = row.variant
    maf = row.minor_AF
    
    add_dict = parse_betas(row.add_betas)
    dom_dict = parse_betas(row.dom_betas)
    
    # Get all unique traits for this variant across both additive and dominance dicts
    traits = set(add_dict.keys()).union(set(dom_dict.keys()))
    
    # Append a flat record for each trait
    for t in traits:
        records.append({
            'variant': variant,
            'maf': maf,
            'trait': t,
            'add_beta': add_dict.get(t, np.nan),
            'dom_beta': dom_dict.get(t, np.nan)
        })

# Create the expanded DataFrame
long_df = pd.DataFrame(records)

# Drop rows where either beta or MAF is missing to ensure clean correlation tests
clean_df = long_df.dropna(subset=['maf', 'add_beta', 'dom_beta'])
print(f"For all SNPs")
calculate_corr(clean_df)


# For d-pleiotropic SNPs
# Load the data
df_pleio = pd.read_csv("/Users/sezgi/Documents/dominance_pleiotropy/loci_level/coloc_results/coloc_snp_info.tsv", sep='\t')
df_pleio = df_pleio.rename(columns={'std_add_b': 'add_beta', 'std_dom_b': 'dom_beta'})

print(f"\nFor d-pleiotropic SNPs")
calculate_corr(df_pleio)