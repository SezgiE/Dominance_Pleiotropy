import pandas as pd


def format_qc_table(df):
    return df.rename(columns={
        'chr': 'Chromosome', 'non-filtered': 'Initial Variants (N)', 
        'filter_MAF': 'MAF > 0.01 (N)', 'filter_HWE': 'HWE P ≥ 1e-6 (N)', 
        'filter_INFO': 'INFO > 0.9 (N)', 'filter_Autosomes': 'Autosomes Only (N)', 
        'filter_INDELs': 'Biallelic SNPs (N)', 'QCed variants': 'Passed QC (N)', 
        'add_sig': 'Add. Sig. (N)', 'add_pleiotropy': 'Add. Pleiotropy (N)', 
        'dom_sig': 'Dom. Sig. (N)', 'dom_pleiotropy': 'Dom. Pleiotropy (N)', 
        'dom_pleiotropy_category': 'Dom. Pleio. Categories (N)'
    })


def snp_desc(snp_info_path, sig_SNPs_path, phen_dict_path):

    df = pd.read_csv(snp_info_path, sep="\t", compression="gzip",
                     usecols=["variant", "chr", "pos", "rsid", "info", "minor_AF", "p_hwe"],
                     dtype={"variant":str, "chr":str, "pos":int, "rsid":str, "info":float, 
                            "minor_AF":float, "p_hwe":float})

    out = pd.DataFrame(index=df['chr'].unique())
    out['raw_total'] = df.groupby('chr').size()

    df = df[df['minor_AF'] > 0.01]
    out['filter_MAF'] = df.groupby('chr').size()

    df = df[df['p_hwe'] >= 1e-6]
    out['filter_HWE'] = df.groupby('chr').size()

    df = df[df['info'] > 0.9]
    out['filter_INFO'] = df.groupby('chr').size()

    df = df[~df['chr'].isin(['X', 'Y', 'XY', 'MT'])]
    out['filter_Autosomes'] = df.groupby('chr').size()

    df = df[df['variant'].str.match(r'^[^:]+:[^:]+:[a-zA-Z]:[a-zA-Z]$', na=False)]
    out['filter_INDELs'] = df.groupby('chr').size()

    df = df.drop_duplicates(subset=['chr', 'pos'], keep=False)
    out['QCed variants'] = df.groupby('chr').size()

    output_df = out.fillna(0).astype(int).reset_index().rename(columns={'index': 'chr'}).sort_values('chr')
    
    sig_df = pd.read_csv(sig_SNPs_path, sep="\t", compression="gzip", 
                        usecols=["variant", "chr", "pos", "add_sig_total", "dom_sig_total",
                                  "sig_add_traits", "sig_dom_traits"],
                        dtype={"variant":str, "chr":str, "pos":int, "add_sig_total":int, "dom_sig_total":int,
                                  "sig_add_traits":str, "sig_dom_traits":str})
    
    names_list = pd.read_excel(phen_dict_path, usecols=["phenotype_code", "category"])

    sig_counts = sig_df.groupby('chr').agg(
        add_sig=('add_sig_total', lambda x: (x > 0).sum()),
        add_pleiotropy=('add_sig_total', lambda x: (x > 1).sum()),
        dom_sig=('dom_sig_total', lambda x: (x > 0).sum()),
        dom_pleiotropy=('dom_sig_total', lambda x: (x > 1).sum())
    )

    pleio_snps = sig_df[sig_df['dom_sig_total'] > 1][['variant', 'chr', 'sig_dom_traits']].copy()
    pleio_snps['sig_dom_traits'] = pleio_snps['sig_dom_traits'].str.replace(r"[()\s]", "", regex=True).str.split(",")

    cat_counts = pleio_snps.explode('sig_dom_traits').merge(names_list, left_on='sig_dom_traits', right_on='phenotype_code', how='inner') \
        .groupby(['chr', 'variant'])['category'].nunique()

    dom_pleio_cat = cat_counts[cat_counts > 1].groupby('chr').size().rename('dom_pleiotropy_category')

    final_sig_counts = sig_counts.join(dom_pleio_cat).fillna(0).astype(int).reset_index()
    output_df = output_df.merge(final_sig_counts, on='chr', how='left').fillna(0)

    output_df['chr'] = output_df['chr'].astype(str)

    output_df['sort_key'] = pd.to_numeric(output_df['chr'], errors='coerce').fillna(float('inf'))
    output_df = output_df.sort_values(['sort_key', 'chr']).drop(columns=['sort_key'])

    numeric_cols = output_df.columns[1:]
    totals = pd.DataFrame([output_df[numeric_cols].sum()], columns=numeric_cols)
    totals.insert(0, 'chr', 'Total')

    output_df = pd.concat([output_df, totals], ignore_index=True)
    output_df[numeric_cols] = output_df[numeric_cols].astype(int)

    final_table = format_qc_table(output_df)

    # save for Excel:
    output_1 = f"{output_path}/SNP_descriptives.xlsx"
    output_2 = f"{output_path}/SNP_descriptives_plotting.xlsx"

    final_table.to_excel(output_1, index=False)
    output_df.to_excel(output_2, index=False)


if __name__ == "__main__":
   
    snp_info_path = "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/variants.tsv.bgz"
    sig_SNPs_path = "/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/significant_SNPs/all_sig_SNPs.tsv.gz"
    phen_dict_path = "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/phen_dict_renamed.xlsx"
    output_path = "/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/results/"

    snp_desc(snp_info_path, sig_SNPs_path, phen_dict_path)