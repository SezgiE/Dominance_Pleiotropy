import textwrap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.patches as mpatches
import matplotlib.transforms as mtransforms


def get_data(data_path, phen_info_path):

    # Load the data
    coloc_merged = pd.read_csv(data_path, sep="\t")
    phen_info = pd.read_excel(phen_info_path)

    active_traits = coloc_merged["phen1"].unique()
    phen_info = phen_info[phen_info["phenotype_code"].isin(active_traits)]

    # Group and summarize
    snp_info = (
        coloc_merged.groupby("variant")["phen1"]
        .agg(n="nunique", phenotypes=lambda x: ", ".join(x.unique()))
        .reset_index()
    )

    # Parse CHR and BP from variant string (e.g., 10:71094504:T:C)
    snp_info[["chr", "pos"]] = (
        snp_info["variant"].str.split(":", n=2, expand=True).iloc[:, 0:2]
    )
    snp_info["chr"] = snp_info["chr"].astype(int)
    snp_info["pos"] = snp_info["pos"].astype(int)
    # Calculate EQUAL WIDTH cumulative positions for the X-axis
    snp_info = snp_info.sort_values(["chr", "pos"]).reset_index(drop=True)
    snp_info["bp_cum"] = 0.0

    # Get sorted unique chromosomes to ensure they plot in 1, 2, 3 order
    chromosomes = sorted(snp_info["chr"].unique())

    for i, chr_val in enumerate(chromosomes):
        mask = snp_info["chr"] == chr_val
        chr_max_bp = snp_info.loc[mask, "pos"].max()

        # Safety check to prevent division by zero
        if chr_max_bp == 0:
            chr_max_bp = 1

        # 1. Normalize the BP to a 0.0 - 1.0 scale
        # 2. Add 'i' to shift it into its dedicated, equally-sized slot
        snp_info.loc[mask, "bp_cum"] = i + (snp_info.loc[mask, "pos"] / chr_max_bp)

    return snp_info, phen_info


def set_style():
    """Hardcodes matplotlib parameters to format the plot."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8,
            "axes.labelsize": 9,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "axes.linewidth": 1.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def plot_coloc_results(snp_info, phen_info, out_dir):

    # Sort by category, then alphabetically
    phen_info = phen_info.sort_values(["category", "phenotype_code"]).reset_index(
        drop=True
    )

    chromosomes = sorted(snp_info['chr'].unique())
    num_chrs = len(chromosomes)

    # Distribute traits evenly across the exact same X-span as the SNPs
    phen_info["trait_x"] = np.linspace(0.05, num_chrs - 0.05, len(phen_info))

    # Colors
    unique_categories = phen_info["category"].unique()
    custom_palette = [
        "#00A05B",
        "#DC0000",
        "#631879",
        "#E64B35",
        "#4DBBD5",
        "#C59316",
        "#3C5488",
    ]

    new_palette = [mcolors.to_rgba(c, alpha=1) for c in custom_palette]
    top_to_bottom_cats = list(reversed(unique_categories))

    category_color_dict = {
        cat: new_palette[i % len(new_palette)]
        for i, cat in enumerate(top_to_bottom_cats)
    }

    # Map the color to the phenotype dataframe so we can easily look it up later
    phen_info["color"] = phen_info["category"].map(category_color_dict)

    # 2. Setup Figure
    trait_y_level = snp_info["n"].max() * 1.15
    fig, ax = plt.subplots(figsize=(10, 5))
    phen_dict = phen_info.set_index("phenotype_code").to_dict("index")

    # 3. Draw Connecting Lines (Network)
    max_traits = snp_info['n'].max()

    for _, row in snp_info.iterrows():
        snp_x = row['bp_cum']
        snp_y = row['n']
        traits = [t.strip() for t in row['phenotypes'].split(',')]
        
        # INCREASE alpha for higher pleiotropy (more traits)
        # Scale from a base transparency (0.05) up to a solid highlight (0.6)
        # Using (n / max_n) ensures the 'hubs' are the most opaque things on the plot
        hub_alpha = 0.05 + (0.55 * (len(traits) / max_traits))
        
        # Also make the lines slightly thicker for hubs to give them more 'weight'
        hub_linewidth = 0.5 + (0.7 * (len(traits) / max_traits))

        for t in traits:
            if t in phen_dict:
                t_info = phen_dict[t]
                ax.plot(
                    [snp_x, t_info["trait_x"]],
                    [snp_y, trait_y_level],
                    color=t_info["color"],
                    alpha=hub_alpha, 
                    linewidth=hub_linewidth,
                    zorder=1,
                    solid_capstyle='round' # Makes line endings cleaner
                )
                # --------------------------------------------------------------------

    # 4. Draw Bottom Axis (SNPs / Manhattan)
    chr_colors = ['#404040', '#A0A0A0'] # Dark grey, Light grey
    
    for i, chr_val in enumerate(chromosomes):
        subset = snp_info[snp_info['chr'] == chr_val]
        color = chr_colors[i % 2]
        
        # 1. Scatter points for this chromosome
        ax.scatter(subset['bp_cum'], subset['n'], color=color, s=15, zorder=3, edgecolors='none')

        # 2. ADD THIS: A colored bar at the base (y=0) to indicate chromosome boundaries
        # We draw a thin rectangle from i to i+1 at the very bottom
        rect = patches.Rectangle((i, 0.0), width=1, height=0.02, 
                                 facecolor=color, edgecolor='none', 
                                 transform=ax.get_xaxis_transform(), clip_on=False, zorder=2)
        ax.add_patch(rect)

    # 3. Clean up the spine
    ax.spines['bottom'].set_visible(False) # Hide the default black line
    
    # Force Equal Geometry for X-axis
    ax.set_xlim(0, num_chrs)
    tick_centers = np.arange(num_chrs) + 0.5
    ax.set_xticks(tick_centers)
    ax.set_xticklabels(chromosomes)
    
    # (Optional) Keep or remove vertical dividers based on preference
    for i in range(1, num_chrs):
        ax.axvline(x=i, color='grey', linestyle=':', linewidth=0.5, alpha=0.4, zorder=0)

    ax.set_xlabel('Chromosome', labelpad=15) # Add padding to avoid overlap with the bar
    ax.set_ylabel('Number of Unique Traits')
    ax.set_ylim(0, trait_y_level * 1.05)

    # 5. Draw Top Axis (Traits)
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(phen_info["trait_x"])
    ax2.set_xticklabels(phen_info["phenotype_code"], rotation=45, ha="left", fontsize=6)

    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.tick_params(axis="x", length=0)

    # --- NEW: Adapted rectangle logic for the top X-axis ---
    # Blends X-axis data coordinates (trait_x) with Y-axis relative coordinates (1.0 = top of plot)
    trans = mtransforms.blended_transform_factory(ax2.transData, ax2.transAxes)

    # Dynamic width based on the X-axis span so blocks touch or space correctly
    block_width = num_chrs * 0.008

    for _, row in phen_info.iterrows():
        rect = patches.Rectangle(
            (row["trait_x"] - (block_width / 2), 1.0),
            width=block_width,
            height=0.03,
            transform=trans,
            facecolor=row["color"],
            edgecolor="none",
            clip_on=False,
            alpha=0.8,
        )
        ax2.add_patch(rect)
    # -------------------------------------------------------

    # 6. Add Category Legend
    # --- CHANGED: Legend now uses your specific category_color_dict ---
    legend_patches = [
        patches.Patch(
            color=category_color_dict[cat],
            label=textwrap.fill(
                str(cat), width=30, break_long_words=False, break_on_hyphens=False
            ),
            linewidth=0.5,
        )
        for cat in unique_categories
    ]

    ax.legend(
        handles=legend_patches,
        title="Phenotype Categories",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )
    # ------------------------------------------------------------------

    plt.tight_layout()
    plt.savefig(f"{out_dir}/pleiotropy_network.pdf", bbox_inches="tight")
    plt.close()


if __name__ == "__main__":

    data_path = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/coloc_results/merged_coloc.tsv"
    phen_info_path = "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/phen_dict_renamed.xlsx"
    out_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/coloc_results"

    snp_info, phen_info = get_data(data_path, phen_info_path)
    set_style()
    plot_coloc_results(snp_info, phen_info, out_dir)
