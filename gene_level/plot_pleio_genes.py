import os
import textwrap
import itertools
import numpy as np
import pandas as pd
import gseapy as gp
from gseapy import Biomart
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from collections import defaultdict
import matplotlib.patches as mpatches
from matplotlib_venn import venn2, venn2_circles

#bm = Biomart()


def set_style():
    """Hardcodes matplotlib parameters to format the plot."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 12,
            "axes.labelsize": 12,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.linewidth": 1.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def std_expression(gtex_med_TPM_path, mapped_genes):

    gtex_med_tpm = pd.read_csv(gtex_med_TPM_path, sep="\t", skiprows=2)
    gtex_med_tpm["ENSG"] = gtex_med_tpm["Name"].str.split(".").str[0]

    gtex_tpm_pleio = gtex_med_tpm[gtex_med_tpm["ENSG"].isin(mapped_genes)].copy()

    # Standardize row-wise for numeric columns
    numeric_cols = gtex_tpm_pleio.select_dtypes(include=[np.number]).columns

    gtex_tpm_pleio_std = gtex_tpm_pleio.loc[
        (gtex_tpm_pleio[numeric_cols] != 0).any(axis=1)
    ].copy()
    gtex_tpm_pleio_std[numeric_cols] = gtex_tpm_pleio[numeric_cols].apply(
        lambda row: (row - row.mean()) / row.std(), axis=1
    )

    return gtex_tpm_pleio_std


def enrichment_getsets(mapped_genes):

    # ENSG to gene symbol Map
    gene_id_dict = {"ensembl_gene_id": mapped_genes}
    results = bm.query(
        dataset="hsapiens_gene_ensembl",
        attributes=["ensembl_gene_id", "external_gene_name", "entrezgene_id", "go_id"],
        filters=gene_id_dict,
    )
    gene_symbols = results["external_gene_name"].dropna().unique().tolist()

    # Run enrichment
    enr = gp.enrichr(
        gene_list=gene_symbols,
        gene_sets=[
            "MSigDB_Hallmark_2020",
            "KEGG_2026",
            "Reactome_Pathways_2024",
            "GO_Biological_Process_2025",
            "GO_Molecular_Function_2025",
            "GO_Cellular_Component_2025",
            "GTEx_Tissues_V8_2023",
        ],
        organism="human",
        outdir=None,
    )

    # Extract results and filter for significance (FDR < 0.05)
    results = enr.results
    sig_results = results[results["Adjusted P-value"] < 0.05].copy()

    # Calculate -log10(FDR) for plotting
    sig_results["log10(FDR)"] = -np.log10(sig_results["Adjusted P-value"])

    # Sort by significance
    sig_results = sig_results.sort_values("log10(FDR)", ascending=True)
    sig_results[["GoIs", "Total"]] = (
        sig_results["Overlap"].str.split("/", n=2, expand=True).iloc[:, :2].astype(int)
    )
    sig_results["Overlap"] = sig_results["GoIs"] / sig_results["Total"]

    return sig_results


def plot_intersection(all_genes_df, output_dir):

    eqtl_genes = set(all_genes_df["gene_id_eqtl"].dropna())
    pos_genes = set(all_genes_df["gene_id_pos"].dropna())

    intersected_genes = eqtl_genes.intersection(pos_genes)
    eqtl_only_genes = eqtl_genes - pos_genes
    pos_only_genes = pos_genes - eqtl_genes

    gene_to_cat = defaultdict(set)
    for col in ["gene_id_eqtl", "gene_id_pos"]:
        subset = all_genes_df[[col, "category"]].dropna()
        for gene, cat in zip(subset[col], subset["category"]):
            
            if isinstance(cat, str) and ',' in cat:
                for c in cat.split(','):
                    gene_to_cat[gene].add(c.strip())
            else:
                gene_to_cat[gene].add(cat)

   
    int_cats = pd.Series([cat for g in intersected_genes if g in gene_to_cat for cat in gene_to_cat[g]]).value_counts()
    eqtl_only_cats = pd.Series([cat for g in eqtl_only_genes if g in gene_to_cat for cat in gene_to_cat[g]]).value_counts()
    pos_only_cats = pd.Series([cat for g in pos_only_genes if g in gene_to_cat for cat in gene_to_cat[g]]).value_counts()


    # Combine into a dataframe and convert to percentages for a clean comparison
    cat_df = pd.DataFrame(
        {
            "Intersected": int_cats,
            "eQTL-only": eqtl_only_cats,
            "Positional-only": pos_only_cats,
        }
    ).fillna(0)
    cat_df.index = [textwrap.fill(str(cat), width=25) for cat in cat_df.index]
    cat_df_pct = cat_df.div(cat_df.sum(axis=0), axis=1) * 100

    set_style()
    fig, axes = plt.subplots(2, 1, figsize=(8, 10), gridspec_kw={"hspace": 0.3})

    # Panel A: Venn Diagram
    v = venn2(
        [eqtl_genes, pos_genes],
        set_labels=("eQTL Mapped Genes", "Positionally Mapped Genes"),
        set_colors=("#E64B35", "#4DBBD5"),
        alpha=0.8,
        ax=axes[0],
    )
    current_x, current_y = v.set_labels[1].get_position()
    v.set_labels[1].set_position((current_x - 0.3, current_y - 0.05))

    c = venn2_circles(
        [eqtl_genes, pos_genes],
        linestyle="solid",
        linewidth=1.0,
        color="black",
        ax=axes[0],
    )
    axes[0].annotate("A.", xy=(-0.372, 1.05), xycoords="axes fraction", 
                    fontweight="bold", fontsize=14, va="top", ha="right")

    # Panel B: Category Distribution (Stacked Bar)
    custom_palette = ["#3C5488", '#DC0000','#00A087','#91D1C2',
                       "#B09C85","#8491B4","#E64B35","#C59316", '#4DBBD5']

    # Transpose so x-axis represents the Groups, and stacks represent Categories
    cat_df_pct.T.plot(
        kind="bar",
        stacked=True,
        ax=axes[1],
        color=custom_palette[: len(cat_df_pct)],
        edgecolor="black",
        linewidth=0.5,
    )

    axes[1].annotate("B.", xy=(-0.1, 1.05), xycoords="axes fraction", 
                    fontweight="bold", fontsize=14, va="top", ha="right")
    
    axes[1].set_ylabel("Percentage of category annotations (%)")
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=0)

    # Clean up axes and move legend outside the plot area
    axes[1].legend(
        title="Gene Category", bbox_to_anchor=(0.98, 1), loc="upper left", frameon=False,
        title_fontproperties={"weight": "bold"})

    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "gene_intersection_panels.pdf")
    plt.savefig(out_path, format="pdf", bbox_inches="tight", transparent=True)
    plt.close()


def plot_pleio_genes(pleio_genes, all_genes, std_exp_data, output_dir):

    # --------------------------- Data preparation for Panel A ---------------------------
    records = []
    for gene, (source, phens) in pleio_genes.items():
        for p in phens:
            records.append({"Gene": gene, "Phenotype": p, "Source": source})
    df_plot = pd.DataFrame(records)

    gene_name_map = all_genes[["gene_id_pos", "gene_name_pos"]].drop_duplicates(
        subset=["gene_id_pos"]
    )

    gene_name_map_eqtl = all_genes[["gene_id_eqtl", "gene_name_eqtl"]].drop_duplicates(
        subset=["gene_id_eqtl"]
    )

    df_plot = df_plot.merge(
        gene_name_map, left_on="Gene", right_on="gene_id_pos", how="left"
    )

    df_plot = df_plot.merge(
        gene_name_map_eqtl, left_on="Gene", right_on="gene_id_eqtl", how="left"
    )

    df_plot["gene_name_pos"] = df_plot["gene_name_pos"].combine_first(
        df_plot["gene_name_eqtl"]
    )

    df_plot["gene_name_pos"] = df_plot["gene_name_pos"].combine_first(df_plot["Gene"])

    df_plot = df_plot.drop(columns=["gene_id_pos"])

    df_plot["Source"] = df_plot["Source"].astype("category")
    df_plot["Source"] = df_plot["Source"].cat.rename_categories(
        {"both": "Both", "positional": "Positional"}
    )

    phen_cat_map = all_genes[["phen_name", "category"]].drop_duplicates(
        subset=["phen_name"]
    )
    df_plot = df_plot.merge(
        phen_cat_map, left_on="Phenotype", right_on="phen_name", how="left"
    )

    # Sort axes for consistent display
    custom_source_order = ["eQTL", "Both", "Positional"]
    df_plot["Source"] = pd.Categorical(
        df_plot["Source"], categories=custom_source_order, ordered=True
    )
    genes = (
        df_plot[["gene_name_pos", "Source"]]
        .drop_duplicates()
        .sort_values(["Source", "gene_name_pos"])["gene_name_pos"]
        .tolist()
    )
    phenotypes = (
        df_plot[["Phenotype", "category"]]
        .drop_duplicates()
        .sort_values(["category", "Phenotype"])["Phenotype"]
        .tolist()
    )

    color_map = {
        "Blood biochemistry": "#3C5488",
        "Blood count": "#DC0000",
        "Medical conditions": "#8491B4",
        "Eyesight": "#91D1C2",
        "Sun exposure": "#C59316",
        "Visual acuity": "#4DBBD5",
        "Spirometry": "#E64B35",
    }

    gene_to_source = df_plot.set_index("gene_name_pos")["Source"].to_dict()
    source_color_map = {"eQTL": "#E0E0E0", "Both": "#878787", "Positional": "#1A1A1A"}

    # --------------------------- Heatmap Data Prep: Panel B ---------------------------
    heatmap_data = (
        std_exp_data.set_index("Description").drop(columns=["ENSG", "Name"]).T
    )
    heatmap_data.index = heatmap_data.index.str.replace("_", " ")
    aligned_genes = [g for g in genes if g in heatmap_data.columns]
    heatmap_data = heatmap_data[aligned_genes]
    data_matrix = heatmap_data.values

    # Plotting
    set_style()
    fig, (ax, ax_b) = plt.subplots(
        1, 2, figsize=(11, 10), gridspec_kw={"wspace": 0.3, "width_ratios": [1, 1.25]}
    )

    # ---------------------------------------------------------
    # PANEL A: Matrix of Pleiotropic Genes and Phenotypes
    # ---------------------------------------------------------
    ax.text(
        -0.5,
        1.05,
        "A.",
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
    )
    all_combinations = list(itertools.product(genes, phenotypes))
    bg_genes = [combo[0] for combo in all_combinations]
    bg_phens = [combo[1] for combo in all_combinations]

    # Plot empty squares
    ax.scatter(
        bg_genes,
        bg_phens,
        facecolors="none",
        edgecolors="black",
        s=150,
        marker="s",
        linewidths=0.5,
        zorder=1,
        alpha=0.2,
    )

    for i in range(len(phenotypes)):
        if i % 2 == 0:
            ax.axhspan(
                i - 0.5, i + 0.5, facecolor="gray", zorder=0, edgecolor="none", alpha=0.1
            )
    
    for j in range(len(genes)):
        if j % 2 == 0:
            ax.axvspan(
                j - 0.5, j + 0.5, facecolor="gray", zorder=0, edgecolor="none", alpha=0.12 
            )

    edge_alpha = 0.2
    for category, color in color_map.items():
        subset = df_plot[df_plot["category"] == category]
        ax.scatter(
            subset["gene_name_pos"],
            subset["Phenotype"],
            c=color,
            s=150,
            marker="s",
            edgecolors=(0, 0, 0, edge_alpha),
            linewidths=0.5,
            zorder=3,
        )

    # Format axes and aesthetics
    ax.set_xticks(range(len(genes)))
    ax.xaxis.tick_top()
    ax.set_xticklabels(genes, rotation=45, ha="left", va="bottom")

    track_height = 0.4
    for i, gene in enumerate(genes):
        source = gene_to_source.get(gene)
        color = source_color_map.get(source, "lightgray")

        rect = mpatches.Rectangle(
            (i - 0.5, -0.5 - track_height),
            width=1.0,
            height=track_height,
            facecolor=color,
            edgecolor="white",
            linewidth=0.5,
            clip_on=False,
        )
        ax.add_patch(rect)

    ax.set_yticks(range(len(phenotypes)))
    ax.set_yticklabels(phenotypes)

    # Replaces the previous ax.set_ylim line to invert the y-axis:

    ax.set_xlim(-0.5, len(genes) - 0.5)
    ax.set_ylim(len(phenotypes) - 0.5, -0.5)

    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(False)

    # Ticks
    ax.tick_params(axis="y", which="both", length=2)
    ax.tick_params(
        axis="x",
        which="major",
        length=8,
        width=0.6,
        direction="out",
        color="black",
        pad=5,
    )

    handles_cat = [
        mpatches.Patch(color=color, label=str(label).capitalize())
        for label, color in color_map.items()
    ]
    legend_cat = ax.legend(
        handles=handles_cat,
        title="Category",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        fontsize=8,
        title_fontproperties={"weight": "bold", "size": 10},
    )
    ax.add_artist(legend_cat)

    handles_source = [
        mpatches.Patch(color=color, label=label)
        for label, color in source_color_map.items()
    ]
    legend_source = ax.legend(
        handles=handles_source,
        title="Mapping",
        bbox_to_anchor=(1.25, -0.05),
        loc="center",
        frameon=False,
        fontsize=8,
        ncol=3,
        title_fontproperties={"weight": "bold", "size": 10},
    )

    # ---------------------------------------------------------
    # PANEL B: Expression Heatmap
    # ---------------------------------------------------------
    ax_b.text(
        -0.15,
        1.05,
        "B.",
        transform=ax_b.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
    )

    v_max = np.nanmax(np.abs(data_matrix))

    c = ax_b.pcolormesh(
        data_matrix,
        cmap="RdYlBu_r",
        vmin=-v_max,
        vmax=v_max,
        edgecolors="white",
        linewidth=0.5,
    )

    for i in range(6, data_matrix.shape[0], 6):
        ax_b.axhline(
            i, color="black", linewidth=0.5, linestyle="--", alpha=0.2, zorder=2
        )

    for j in range(3, data_matrix.shape[1], 3):
        ax_b.axvline(
            j, color="black", linewidth=0.5, linestyle="--", alpha=0.2, zorder=2
        )

    ax_b.invert_yaxis()

    ax_b.set_xticks(np.arange(data_matrix.shape[1]) + 0.5)
    ax_b.set_yticks(np.arange(data_matrix.shape[0]) + 0.5)
    ax_b.yaxis.tick_right()

    ax_b.xaxis.tick_top()

    ax_b.set_xticklabels(aligned_genes, rotation=45, ha="left", va="bottom")
    ax_b.set_yticklabels(heatmap_data.index)

    ax_b.tick_params(axis="x", top=True, bottom=False, length=8, width=0.6)
    ax_b.tick_params(axis="y", left=False, right=True, length=3, width=0.5)

    for i, gene in enumerate(aligned_genes):
        source = gene_to_source.get(gene)
        color = source_color_map.get(source, "lightgray")

        rect = mpatches.Rectangle(
            (i, 0),
            width=1.0,
            height=-track_height - 0.4,
            facecolor=color,
            edgecolor="white",
            linewidth=0.5,
            clip_on=False,
        )
        ax_b.add_patch(rect)

    ax_b.spines["top"].set_visible(True)
    ax_b.spines["right"].set_visible(True)
    ax_b.spines["left"].set_visible(False)
    ax_b.spines["bottom"].set_visible(False)

    # Add colorbar
    cbar = fig.colorbar(c, ax=ax_b, location="left", shrink=0.3, aspect=30, pad=0.05)
    cbar.outline.set_visible(False)

    cbar.set_label("Median Expression (Z-score)", rotation=90, labelpad=5, fontsize=8)
    cbar.outline.set_linewidth(0.2)

    ax_b.set_ylabel("")
    ax_b.set_xlabel("")

    plt.savefig(
        f"{output_dir}/pleiotropy_matrix.pdf",
        dpi=600,
        bbox_inches="tight",
        bbox_extra_artists=(legend_cat, legend_source),
    )


def plot_geneset(enrich_path, output_dir):

    df = pd.read_csv(enrich_path, sep="\t")

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

    plt.tight_layout()
    plt.savefig(
        f"{output_dir}/pleio__genes_enrichment_plot.pdf", dpi=600, bbox_inches="tight"
    )


if __name__ == "__main__":

    gtex_med_TPM_path = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_expression/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct"
    all_genes_path = (
        "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/genes_all/genes_all.tsv"
    )
    all_genes_df = pd.read_csv(all_genes_path, sep="\t")

    output_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/plots"

    # Intersection between MAGMA and eQTL
    plot_intersection(all_genes_df, output_dir)

    # Gene set enrichment for pleiotropic genes
    pleiotropic_pos_df = all_genes_df.groupby("gene_id_pos").filter(
        lambda group: group["category"].nunique() > 1
    )
    pleiotropic_eqtl_df = all_genes_df.groupby("gene_id_eqtl").filter(
        lambda group: group["category"].nunique() > 1
    )

    pos_pleio_genes = set(pleiotropic_pos_df["gene_id_pos"].dropna().unique())
    eqtl_pleio_genes = set(pleiotropic_eqtl_df["gene_id_eqtl"].dropna().unique())
    all_pleio_genes = pos_pleio_genes.union(eqtl_pleio_genes)

    pos_phenotypes = (
        pleiotropic_pos_df.groupby("gene_id_pos")["phen_name"]
        .unique()
        .apply(list)
        .to_dict()
    )
    eqtl_phenotypes = (
        pleiotropic_eqtl_df.groupby("gene_id_eqtl")["phen_name"]
        .unique()
        .apply(list)
        .to_dict()
    )

    pleio_genes = {}
    for gene_id in all_pleio_genes:
        is_pos = gene_id in pos_pleio_genes
        is_eqtl = gene_id in eqtl_pleio_genes

        if is_pos and is_eqtl:
            source = "both"

            phen_list_pos = pos_phenotypes.get(gene_id, [])
            phen_list_eqtl = eqtl_phenotypes.get(gene_id, [])

            phenotypes = sorted(list(set(phen_list_pos + phen_list_eqtl)))
        elif is_pos:
            source = "positional"
            phenotypes = sorted(pos_phenotypes.get(gene_id, []))
        else:  # is_eqtl must be true
            source = "eQTL"
            phenotypes = sorted(eqtl_phenotypes.get(gene_id, []))

        pleio_genes[gene_id] = (source, phenotypes)

    std_exp_data = std_expression(gtex_med_TPM_path, list(pleio_genes.keys()))

    plot_pleio_genes(pleio_genes, all_genes_df, std_exp_data, output_dir)
    
    # # enrich_res = enrichment_getsets(list(pleio_genes.keys()))
    # # enrich_path = os.path.join(output_dir, "pleio_genes_enrichment_results.tsv")
    # # enrich_res.to_csv(enrich_path, sep="\t", index=False)

    # # plot_geneset(enrich_path, output_dir)
    
