import pandas as pd
from pathlib import Path


def in_mhc(item, is_loci=False):
    parts = item.split(':')
    if parts[0] == '6':
        if is_loci:
            start = int(parts[1])
            end = int(parts[2])
            return (start <= 34000000) and (end >= 25000000)
        else:
            pos = int(parts[1])
            return 25000000 <= pos <= 34000000
    return False


all_snps_path = Path("/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/significant_SNPs")
sig_loci_path = Path("/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sig_loci")

all_snps = pd.read_csv(all_snps_path / "all_sig_snps.tsv.gz", sep='\t', 
                       usecols=["variant", "dom_sig_total"])
pleio_snps = set(all_snps[all_snps['dom_sig_total'] > 1]['variant'])


dfs = {}
idp_snps = set()
idp_lead_snps = set()
idp_loci = set()
snps_to_remove = {"1:107603330:C:G",  "1:25736299:A:T" , "17:79628408:C:G",
                  "3:121715432:T:C", "3:121712980:C:T", "3:121712051:A:C"}

for file_path in sig_loci_path.glob("*_sig_loci.tsv"):

    phen_code = file_path.name.replace("_sig_loci.tsv", "")
    
    dfs[phen_code] = pd.read_csv(file_path, sep='\t', usecols=["indep_id", "lead_id", "ld_id"])

    
    idp_snps.update(dfs[phen_code]["indep_id"].tolist())
    idp_lead_snps.update(dfs[phen_code]["lead_id"].tolist())
    idp_loci.update(dfs[phen_code]["ld_id"].tolist())


idp_snps = {x for x in idp_snps if not in_mhc(x, is_loci=False)}
idp_lead_snps = {x for x in idp_lead_snps if not in_mhc(x, is_loci=False)}
idp_loci = {x for x in idp_loci if not in_mhc(x, is_loci=True)}


idp_snps.difference_update(snps_to_remove)
idp_lead_snps.difference_update(snps_to_remove)

pleio_idp_snps = idp_snps.intersection(pleio_snps)
pleio_idp_lead_snps = idp_lead_snps.intersection(pleio_snps)

summary_dict = {
    'idp_snps': pd.Series([item.split(':')[0] for item in idp_snps]).value_counts(),
    'idp_lead_snps': pd.Series([item.split(':')[0] for item in idp_lead_snps]).value_counts(),
    'idp_loci': pd.Series([item.split(':')[0] for item in idp_loci]).value_counts(),
    'pleio_idp_snps': pd.Series([item.split(':')[0] for item in pleio_idp_snps]).value_counts(),
    'pleio_idp_lead_snps': pd.Series([item.split(':')[0] for item in pleio_idp_lead_snps]).value_counts()
}

df_summary = pd.DataFrame(summary_dict).fillna(0).astype(int)

df_summary = df_summary.reset_index().rename(columns={'index': 'Chromosome'})

df_summary = df_summary.assign(
    sort_key=pd.to_numeric(df_summary['Chromosome'], errors='coerce').fillna(999)
).sort_values('sort_key').drop(columns='sort_key').reset_index(drop=True)

df_summary.loc[len(df_summary)] = df_summary.sum(numeric_only=True)
df_summary.at[len(df_summary)-1, 'Chromosome'] = 'Total'

df_summary.to_csv("/Users/sezgi/Documents/dominance_pleiotropy/loci_level/loci_results/loci_description_summary.tsv", sep='\t', index=False)