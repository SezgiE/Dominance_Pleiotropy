import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches


def merge_magma_genes(magma_res_path, manuel_res_path):

    # Mapped genes from MAGMA
    mapped_genes = set()

    with open(magma_res_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) > 2:
                mapped_genes.add(parts[0])

    with open(manuel_res_path, "r") as f:
        for line_m in f:
            if line_m.startswith("#"):
                continue
            parts_m = line_m.strip().split()
            if len(parts_m) > 2:
                mapped_genes.add(parts_m[0])

    return mapped_genes


def set_style():
    """Hardcodes matplotlib parameters to format the plot."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 12,
            "axes.labelsize": 12,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
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


def plot_geneset(magma_gene_enrich_path, output_dir):

    df = pd.read_csv(magma_gene_enrich_path, sep="\t")

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
        f"{output_dir}/magma_geneset_enrichment_plot.pdf", dpi=600, bbox_inches="tight"
    )


def plot_heat_and_enrich(df, magma_summary_enrich_path, output_dir):

    # Load Enrichment Data 
    enrich_df = pd.read_csv(
        magma_summary_enrich_path,
        sep="\t",
        usecols=[
            "type",
            "group_name",
            "mapped_in_group",
            "p_value",
            "fdr_corrected_p",
            "is_significant",
        ],
    )

    # Filter for Chromosomes and drop 0 counts
    chr_enrich = enrich_df[
        (enrich_df["type"] == "Chromosome") & (enrich_df["mapped_in_group"] > 0)
    ].copy()
    chr_enrich = enrich_df[
        (enrich_df["type"] == "Chromosome") & (enrich_df["mapped_in_group"] > 0)
    ].copy()
    chr_enrich["group_name"] = chr_enrich["group_name"].astype(int)
    chr_enrich = chr_enrich.sort_values("group_name", ascending=True)

    # Filter for Biotypes and drop 0 counts
    bio_enrich = enrich_df[
        (enrich_df["type"] == "Gene biotype") & (enrich_df["mapped_in_group"] > 0)
    ].copy()
    bio_enrich["group_name"] = bio_enrich["group_name"].str.replace("_", " ")
    bio_enrich = bio_enrich.sort_values("mapped_in_group", ascending=False)

    # Heatmap Data Prep
    heatmap_data = df.set_index("Description").drop(columns=["ENSG", "Name"]).T
    heatmap_data.index = heatmap_data.index.str.replace("_", " ")
    heatmap_data = heatmap_data.sort_index(axis=0).sort_index(axis=1)
    data_matrix = heatmap_data.values
    print("Data matrix shape:", data_matrix.shape)
    
    # Setup Layout
    set_style()
    fig = plt.figure(
        figsize=(len(heatmap_data.columns) * 0.21, len(heatmap_data) * 0.15 + 10)
    )
    gs_top = fig.add_gridspec(
        1, 2, top=0.95, bottom=0.80, left=-0.05, right=0.95, wspace=0.12
    )

    # gs_bottom
    gs_bottom = fig.add_gridspec(1, 1, top=0.73, bottom=0.05, left=0.1, right=0.95)

    ax_a = fig.add_subplot(gs_top[0, 0])  # Top Left
    ax_b = fig.add_subplot(gs_top[0, 1])  # Top Right
    ax_c = fig.add_subplot(gs_bottom[0, 0])  # Bottom

    # ---------------------------------------------------------
    # PANEL A: Chromosomes
    # ---------------------------------------------------------
    x_pos_a = np.arange(len(chr_enrich))
    
    # Plot the bars using the numerical x positions
    ax_a.bar(x_pos_a, chr_enrich['mapped_in_group'], 
             color='#3C5488', edgecolor='black', linewidth=1)
             
    # Add significance markers centered above vertical bars
    for i in range(len(chr_enrich)):
        if chr_enrich['is_significant'].iloc[i]:
            ax_a.text(x_pos_a[i], chr_enrich['mapped_in_group'].iloc[i] + (chr_enrich['mapped_in_group'].max() * 0.02), 
                      '**', ha='center', va='bottom', fontsize=12, color='black')
        elif chr_enrich['p_value'].iloc[i] < 0.05:
            ax_a.text(x_pos_a[i], chr_enrich['mapped_in_group'].iloc[i] + (chr_enrich['mapped_in_group'].max() * 0.02), 
                      '*', ha='center', va='bottom', fontsize=12, color='black')
            
    ax_b.text(
        -0.15,
        0.95,
        "** FDR-adjusted p < 0.05\n* Nominal p < 0.05",
        transform=ax_b.transAxes,
        ha="right",
        va="top",
        fontsize=10,
        )

    # Apply the chromosome names as x-tick labels
    ax_a.set_xticks(x_pos_a)
    ax_a.set_xticklabels(chr_enrich['group_name'].astype(str), rotation=0)
    ax_a.set_ylabel('Mapped Genes (N)')
    ax_a.set_xlabel('Chromosomes')
   
    ax_c.text(-14, -29, "A.", fontsize=16, fontweight="bold", va="top")


    # ---------------------------------------------------------
    # PANEL B: Biotypes
    # ---------------------------------------------------------
    bio_enrich["group_name"] = bio_enrich["group_name"].str.title()

    x_pos_b = np.arange(len(bio_enrich))

    colors = ["#FF9F1C", "#00A087", "#E64B35", "#4DBBD5", "#7E818D", "#3C5488"]
    biotype_colors = {
        "Protein Coding": "#00A087",
        "Lncrna": "#FF9F1C",
        "Transcribed Unitary Pseudogene": "#E64B35",
        "Tec": "#4DBBD5",
        "Processed Pseudogene": "#8B1C62",
        "Transcribed Processed Pseudogene": "#B484AA",
        "Scarna": "#3C5488",
    }

    for i in range(len(bio_enrich)):
        biotype_name = bio_enrich["group_name"].iloc[i]
        bar_color = biotype_colors.get(biotype_name, "#7E818D")

        ax_b.bar(
            x_pos_b[i],
            bio_enrich["mapped_in_group"].iloc[i],
            color=bar_color,
            width=0.7,
            edgecolor="black",
            linewidth=1.0,
            label=biotype_name,
        )

        # Add significance markers centered above vertical bars
        if bio_enrich["is_significant"].iloc[i]:
            ax_b.text(
                x_pos_b[i],
                bio_enrich["mapped_in_group"].iloc[i]
                + (bio_enrich["mapped_in_group"].max() * 0.02),
                "**",
                ha="center",
                va="bottom",
                fontsize=12,
                color="black",
            )
        elif bio_enrich["p_value"].iloc[i] < 0.05:
            ax_b.text(
                x_pos_b[i],
                bio_enrich["mapped_in_group"].iloc[i]
                + (bio_enrich["mapped_in_group"].max() * 0.02),
                "*",
                ha="center",
                va="bottom",
                fontsize=12,
                color="black",
            )

    # Add the Legend
    main_legend = ax_b.legend(
        frameon=True,
        edgecolor="black",
        fancybox=False,
        framealpha=1.0,
        loc="upper right",
        ncol=1,
        fontsize=10,
    )
    ax_b.add_artist(main_legend)

    # Add the p-value key
    ax_b.text(
        0.98,
        0.45,
        "** FDR-adjusted p < 0.05\n* Nominal p < 0.05",
        transform=ax_b.transAxes,
        ha="right",
        va="top",
        fontsize=10,
    )

    ax_b.set_xticks(x_pos_b)
    ax_b.set_xticklabels([])  # Hides the x-axis text

    # Applied your custom labels
    ax_b.set_ylabel("Mapped Genes (N)")
    ax_b.set_xlabel("Gene Biotypes")

    # Title
    ax_c.text(25, -29, "B.", fontsize=16, fontweight="bold", va="top")

    # ---------------------------------------------------------
    # PANEL C: Heatmap
    # ---------------------------------------------------------
    v_max = np.nanmax(np.abs(data_matrix))

    c = ax_c.pcolormesh(
        data_matrix,
        cmap="RdYlBu_r",
        vmin=-v_max,
        vmax=v_max,
        edgecolors="white",
        linewidth=0.5,
    )

    for i in range(5, data_matrix.shape[0], 5):
        ax_c.axhline(
            i, color="black", linewidth=0.5, linestyle="--", alpha=0.2, zorder=2
        )

    for j in range(5, data_matrix.shape[1], 5):
        ax_c.axvline(
            j, color="black", linewidth=0.5, linestyle="--", alpha=0.2, zorder=2
        )

    ax_c.invert_yaxis()

    ax_c.set_xticks(np.arange(data_matrix.shape[1]) + 0.5)
    ax_c.set_yticks(np.arange(data_matrix.shape[0]) + 0.5)

    ax_c.xaxis.tick_top()

    ax_c.set_xticklabels(heatmap_data.columns, rotation=45, ha="left")
    ax_c.set_yticklabels(heatmap_data.index)

    ax_c.tick_params(axis="x", top=True, bottom=False, length=5, width=0.6)
    ax_c.tick_params(axis="y", left=True, right=False, length=3, width=0.5)

    for spine in ax_c.spines.values():
        spine.set_visible(False)

    # Add colorbar at the bottom
    cbar = fig.colorbar(
        c, ax=ax_c, orientation="horizontal", shrink=0.4, aspect=30, pad=0.02
    )
    cbar.outline.set_visible(False)
    cbar.set_label("Median Expression (Z-score)")
    cbar.outline.set_linewidth(0.2)

    ax_c.set_ylabel("")
    ax_c.set_xlabel("")
    ax_c.text(-14, -6, "C.", fontsize=16, fontweight="bold", va="top")

    # ---------------------------------------------------------
    # Final Output
    # ---------------------------------------------------------
    # tight_layout handles gridspecs very well, but we don't need it to mess up our hspace
    plt.savefig(
        f"{output_dir}/magma_expression_heatmap_combined.pdf",
        dpi=600,
        bbox_inches="tight",
    )


if __name__ == "__main__":

    magma_gene_enrich_path = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/magma/magma_gene_enrich_res.tsv"
    magma_summary_enrich_path = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/magma/magma_enrich_summary.tsv"
    magma_res_path = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/magma/magma_v1/magma_pleio_mapping.genes.annot"
    manuel_res_path = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/magma/magma_v1/manuel_pleio_mapping.genes.annot"
    gtex_med_TPM_path = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_expression/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct"

    output_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/plots"

    plot_geneset(magma_gene_enrich_path, output_dir)

    mapped_genes = merge_magma_genes(magma_res_path, manuel_res_path)
    std_exp_data = std_expression(gtex_med_TPM_path, mapped_genes)
    plot_heat_and_enrich(std_exp_data, magma_summary_enrich_path, output_dir)
