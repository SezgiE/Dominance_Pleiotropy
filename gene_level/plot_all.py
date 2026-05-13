import os
import textwrap
import numpy as np
import pandas as pd
import gseapy as gp
from gseapy import Biomart
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
from matplotlib_venn import venn2, venn2_circles

bm = Biomart()


def set_style():
    """Hardcodes matplotlib parameters to format the plot."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 12,
        'axes.labelsize': 12,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'axes.linewidth': 1.0,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'pdf.fonttype': 42,
        'ps.fonttype': 42
    })


def enrichment_getsets(mapped_genes):
    
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
    
    return sig_results


def plot_intersection(all_genes_df, output_dir):

    eqtl_genes = set(all_genes_df["gene_id_eqtl"].dropna())
    pos_genes = set(all_genes_df["gene_id_pos"].dropna())

    intersected_genes = eqtl_genes.intersection(pos_genes)
    eqtl_only_genes = eqtl_genes - pos_genes
    pos_only_genes = pos_genes - eqtl_genes

    gene_to_cat = {}
    for col in ["gene_id_eqtl", "gene_id_pos"]:
        subset = all_genes_df[[col, "category"]].dropna()
        gene_to_cat.update(dict(zip(subset[col], subset["category"])))


    int_cats = pd.Series([gene_to_cat[g] for g in intersected_genes if g in gene_to_cat]).value_counts()
    eqtl_only_cats = pd.Series([gene_to_cat[g] for g in eqtl_only_genes if g in gene_to_cat]).value_counts()
    pos_only_cats = pd.Series([gene_to_cat[g] for g in pos_only_genes if g in gene_to_cat]).value_counts()

   
    # Combine into a dataframe and convert to percentages for a clean comparison
    cat_df = pd.DataFrame({'Intersected': int_cats, 'eQTL-only': eqtl_only_cats, 'Positional-only': pos_only_cats}).fillna(0)
    cat_df.index = [textwrap.fill(str(cat), width=25) for cat in cat_df.index]
    cat_df_pct = cat_df.div(cat_df.sum(axis=0), axis=1) * 100

    set_style()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), gridspec_kw={'width_ratios': [2, 2.5], 'wspace': 0.3})
    
    # Panel A: Venn Diagram
    v = venn2([eqtl_genes, pos_genes], 
              set_labels=('eQTL Mapped Genes', 'Positionally Mapped Genes'),
              set_colors=('#E64B35', '#4DBBD5'), 
              alpha=0.8, 
              ax=axes[0])
    current_x, current_y = v.set_labels[1].get_position()
    v.set_labels[1].set_position((current_x - 0.3, current_y - 0.05))
    
    c = venn2_circles([eqtl_genes, pos_genes], 
                      linestyle='solid', linewidth=1.0, color='black', ax=axes[0])
    axes[0].set_title('A.', loc='left', fontweight='bold', fontsize=14, x=-0.1)

    # Panel B: Category Distribution (Stacked Bar)
    npg_palette = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', 
                   '#8491B4', '#91D1C2', '#DC0000', '#7E6148', '#B09C85']
    
    # Transpose so x-axis represents the Groups, and stacks represent Categories
    cat_df_pct.T.plot(kind='bar', stacked=True, ax=axes[1], color=npg_palette[:len(cat_df_pct)], edgecolor='black', linewidth=0.5)
    
    axes[1].set_title('B.', loc='left', fontweight='bold', fontsize=14, x=-0.1)
    axes[1].set_ylabel('Percentage of genes (%)')
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=0)
    
    # Clean up axes and move legend outside the plot area
    axes[1].legend(title='Gene Category', bbox_to_anchor=(0.98, 1), loc='upper left', frameon=False)

    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, 'gene_intersection_panels.pdf')
    plt.savefig(out_path, format='pdf', bbox_inches='tight', transparent=True)
    plt.close()


def plot_geneset(df, category_label, output_dir):

    # Data Prep
    df["Term"] = df["Term"].str.replace(r"\s*\([^)]*\)", "", regex=True)
    mask = df["Gene_set"] == "KEGG_2026"
    df.loc[mask, "Term"] = df.loc[mask, "Term"].str.title()
    df["Gene_set"] = df["Gene_set"].str.replace("_", " ")
    df_sorted = df.sort_values(
        ["Gene_set", "Term"], ascending=[False, False]
    ).reset_index(drop=True)

    colors = ["#FF9F1C", "#00A087", "#E64B35", "#4DBBD5", "#7E818D", "#3C5488"]
    unique_dbs = df_sorted["Gene_set"].unique()
    db_color_map = {db: colors[i % len(colors)] for i, db in enumerate(unique_dbs)}
    row_colors = df_sorted["Gene_set"].map(db_color_map)

    all_genes = set()

    for genes_str in df_sorted["Genes"]:
        all_genes.update(str(genes_str).split(";"))
    unique_genes = sorted(list(all_genes))

    # Initialize the plot
    set_style()
    fig, (ax_l, ax_r, ax_m) = plt.subplots(
        1,
        3,
        figsize=(16, len(df_sorted) * 0.18),
        sharey=True,
        # 2. Lowered the matrix ratio from 3 to 1.5 (Tweak this specific number!)
        gridspec_kw={"width_ratios": [1.5, 1.5, 2.8], "wspace": 0},
    )

    # Horizontal background shading
    for i in range(len(df_sorted)):
        row_alpha = 0.15 if i % 2 == 0 else 0.05
        for ax in [ax_l, ax_r, ax_m]:
            ax.axhspan(i - 0.4, i + 0.4, color="gray", alpha=row_alpha, zorder=0, lw=0)

    # Overlap histogram
    ax_l.barh(
        range(len(df_sorted)),
        -df_sorted["Overlap"],
        color=row_colors,
        edgecolor="black",
        linewidth=1,
        height=0.7,
        zorder=3,
    )
    ax_l.set_xlim(right=0)
    ax_l.set_xlabel("Overlap Proportion", fontsize=12, labelpad=10)
    ax_l.xaxis.set_major_formatter(
        ticker.FuncFormatter(lambda x, pos: "" if x == 0 else f"{abs(x):.2f}")
    )

    # -log10(FDR) plot
    ax_r.scatter(
        df_sorted["log10(FDR)"],
        range(len(df_sorted)),
        color=row_colors,
        edgecolor="black",
        linewidth=1,
        s=30,
        zorder=3,
    )
    ax_r.set_xlim(left=0, right=df_sorted["log10(FDR)"].max() + 1)

    ax_r.set_xlabel(r"$-\log_{10}(\text{FDR-adj. p < 0.05})$", fontsize=12, labelpad=10)

    # Y-labels
    ax_l.set_yticks(range(len(df_sorted)))
    ax_l.set_yticklabels(df_sorted["Term"], ha="right", fontsize=9)

    # The shared spine in the middle
    ax_l.spines["right"].set_visible(True)
    ax_r.spines["left"].set_visible(True)

    # Ensure a vertical line exists at x=0 for both
    ax_l.axvline(0, color="black", linewidth=1, zorder=4)
    ax_r.axvline(0, color="black", linewidth=1, zorder=4)

    # Matrix Panel (Far Right)
    for i, row in df_sorted.iterrows():
        row_genes = str(row["Genes"]).split(";")
        for j, gene in enumerate(unique_genes):
            if gene in row_genes:
                ax_m.scatter(
                    j,
                    i,
                    color=row_colors[i],
                    marker="s",
                    s=50,
                    edgecolor="black",
                    linewidth=1,
                    zorder=3,
                )

    # Alternating vertical background shading
    for j in range(len(unique_genes)):
        if j % 2 == 0:  # Only shades every other column
            ax_m.axvspan(j - 0.5, j + 0.5, color="gray", alpha=0.08, zorder=0, lw=0)

    # Configure the matrix
    ax_m.set_xticks(range(len(unique_genes)))
    ax_m.set_xticklabels(unique_genes, rotation=90, ha="center", fontsize=8)
    ax_m.xaxis.tick_bottom()
    ax_m.tick_params(axis="x", direction="inout", length=8, width=1)

    # Lock the limits
    ax_m.set_aspect("equal")
    ax_m.set_xlim(-0.5, len(unique_genes) - 0.5)
    ax_l.set_ylim(-1, len(df_sorted) - 0.5)

    # Final cleanup
    for ax in [ax_l, ax_r, ax_m]:
        for spine in ["top", "right", "left", "bottom"]:
            if ax == ax_m and spine == "top":
                continue
            if spine != "bottom":
                ax.spines[spine].set_visible(False)

        ax.tick_params(axis="y", left=False, right=False)

    # Create legend
    legend_handles = [
        mpatches.Patch(
            facecolor=db_color_map[db], edgecolor="black", linewidth=1, label=db
        )
        for db in sorted(db_color_map.keys())
    ]

    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.6, 0.99),
        ncol=len(db_color_map),
        frameon=False,
        fontsize=10,
        handlelength=1.2,
        columnspacing=2.0,
    )

    
    fig.suptitle(category_label, fontsize=16, fontweight='bold', y=1.02)
    category_label = category_label.replace(" ", "_").replace("/", "_").lower()

    plt.tight_layout()
    plt.savefig(
        f"{output_dir}/enrich_by_category/{category_label}_geneset.pdf", dpi=600, bbox_inches="tight"
    )
    

if __name__ == "__main__":

    all_genes_path = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/genes_all/genes_all.tsv"
    all_genes_df = pd.read_csv(all_genes_path, sep='\t')
    
    output_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/plots"

    # Intersection between MAGMA and eQTL
    #plot_intersection(all_genes_df, output_dir)


    # Gene set enrichment by category
    categories = all_genes_df["category"].dropna().unique()
    
    enrich_results = {}
    enrich_output_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/genes_all"
    
    for category in categories:

        category_df = all_genes_df[all_genes_df["category"] == category]
        category_genes = category_df["gene_id_eqtl"].dropna().unique().tolist() + category_df["gene_id_pos"].dropna().unique().tolist()
        category_genes = list(set(category_genes))
        
        enrich_res = enrichment_getsets(category_genes)
        enrich_results[category] = enrich_res
        
        plot_geneset(enrich_res, category, enrich_output_dir)

    enrich_results_df = pd.concat(enrich_results, names=['Category', 'Index']).reset_index(level=0)
    enrich_results_df.to_csv(os.path.join(enrich_output_dir, "geneset_category_results.tsv"), sep='\t', index=False)