import pandas as pd
import glob
import os


directory_path = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/mult_results"
file_pattern = os.path.join(directory_path, "*_trait.xlsx")

all_files = glob.glob(file_pattern)

df_list = []

for file in all_files:

    df = pd.read_excel(file)
    df_list.append(df)

merged_df = pd.concat(df_list, ignore_index=True)

output_path = os.path.join(directory_path, "mult_by_traits.xlsx")
merged_df.to_excel(output_path, index=False)


# By Chromosome
directory_path = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/mult_results"
file_pattern = os.path.join(directory_path, "*_mult_pleio.xlsx")

all_pleio_files = glob.glob(file_pattern)
df_pleios = []

for file in all_pleio_files:

    df = pd.read_excel(file)
    df_pleios.append(df)

merged_df = pd.concat(df_pleios, ignore_index=True)

mult_SNPs = set(merged_df["variant"].tolist())

df_mult = pd.DataFrame(mult_SNPs, columns=["variant"])
df_mult["mult"] = 1

sig_df = pd.read_csv(
    "/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/significant_SNPs/all_sig_SNPs.tsv.gz",
    sep="\t",
    compression="gzip",
    usecols=[
        "variant",
        "rsid",
        "chr",
        "pos",
        "add_sig_total",
        "dom_sig_total",
        "sig_add_traits",
        "sig_dom_traits",
    ],
    dtype={
        "variant": str,
        "rsid": str,
        "chr": str,
        "pos": int,
        "add_sig_total": int,
        "dom_sig_total": int,
        "sig_add_traits": str,
        "sig_dom_traits": str,
    },
)

merged_df = sig_df.merge(df_mult, on="variant", how="left")
merged_df["mult"] = merged_df["mult"].fillna(0).astype(int)

merged_df["is_multi"] = (merged_df["dom_sig_total"] > 1) & (merged_df["mult"] > 0)

sig_counts = (
    merged_df.groupby("chr")
    .agg(
        dom_sig=("dom_sig_total", lambda x: (x > 0).sum()),
        dom_pleiotropy=("dom_sig_total", lambda x: (x > 1).sum()),
        multi_pleio=("is_multi", "sum"),
        mult_var_ids=(("variant"), lambda x: ", ".join(x[merged_df["is_multi"] == True].tolist())),
        mult_var_rsids=(("rsid"), lambda x: ", ".join(x[merged_df["is_multi"] == True].tolist()))
    )
    .reset_index()
)

sig_counts["Ratio"] = sig_counts["multi_pleio"] / sig_counts["dom_pleiotropy"]

sig_counts["chr"] = sig_counts["chr"].astype(int)
sig_counts = sig_counts.sort_values("chr").reset_index(drop=True)

output_path = os.path.join(directory_path, "mult_pleio.xlsx")
sig_counts.to_excel(output_path, index=False)
