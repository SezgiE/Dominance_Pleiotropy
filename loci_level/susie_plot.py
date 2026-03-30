import os
import sys
import colorsys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy.stats import norm

import glob

def get_susie_data(phen_code, susie_dir):
    """Finds all individual locus SuSiE files for a phenotype and merges them in memory."""
    
    # Search for all files that start with the phen_code and end with _susie_res.tsv
    # Adjust this search pattern if your file naming convention is slightly different!
    search_pattern = os.path.join(susie_dir, f"{phen_code}_susie_res.tsv")
    locus_files = glob.glob(search_pattern)
    
    if not locus_files:
        return None
        
    print(f"  -> Found {len(locus_files)} locus files for {phen_code}. Merging in memory...")
    
    # Read and concatenate all found files
    df_list = []
    for file in locus_files:
        try:
            temp_df = pd.read_csv(file, sep='\t')
            df_list.append(temp_df)
        except Exception as e:
            print(f"  -> Error reading {file}: {e}")
            
    if not df_list:
        return None
        
    df = pd.concat(df_list, ignore_index=True)
    
    # Calculate P-values directly from the Z-score for plotting
    if 'dom_z_score' in df.columns:
        df['dom_log10_pval'] = -np.log10(np.clip(2 * norm.sf(np.abs(df['dom_z_score'])), a_min=1e-300, a_max=1.0))
    
    # Apply your MHC mask
    mhc_mask = (
        (df['chr'] == 6) & 
        (df['pos'] >= 25000000) & 
        (df['pos'] <= 34000000)
    )
    df = df[~mhc_mask].reset_index(drop=True).copy()
    
    return df

def get_gene_annotations(gene_bed_path):
    """Loads a standard BED file containing: chr, start, end, gene_name"""
    genes_df = pd.read_csv(gene_bed_path, sep='\t', header=None, 
                           names=['chr', 'start', 'end', 'gene_name'])
    genes_df['chr'] = genes_df['chr'].astype(str).str.replace('chr', '')
    genes_df['chr'] = pd.to_numeric(genes_df['chr'], errors='coerce')
    return genes_df

def set_style():
    """Hardcodes matplotlib parameters to format the plot."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8,
        'axes.labelsize': 9,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'axes.linewidth': 1.0,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'pdf.fonttype': 42,
        'ps.fonttype': 42
    })

def plot_regional_association(plot_df, genes_df, phen_code, phen_name, output_dir, buffer_kb=50):
    """Generates 3-panel SuSiE LocusZoom plots for each locus."""
    
    set_style()
    os.makedirs(output_dir, exist_ok=True)

    # SuSiE outputs use locus_id instead of ld_id
    unique_blocks = plot_df['locus_id'].dropna().unique()
    n_blocks = len(unique_blocks)
    print(f"Found {n_blocks} unique loci to plot.")

    if n_blocks == 0:
        return None

    plot_df = plot_df.sort_values("dom_log10_pval", ascending=False).reset_index(drop=True)
    plot_df["-log10_plotting"] = np.where(
        plot_df["dom_log10_pval"] > 300, 
        300 + (plot_df["dom_log10_pval"] / 100),
        plot_df["dom_log10_pval"]
    )

    ncols = 2
    nrows = int(np.ceil(n_blocks / ncols))

    # Slightly increased height per row to accommodate the new 3rd PIP panel
    fig = plt.figure(figsize=(6.5 * ncols, 5.5 * nrows))
    master_gs = GridSpec(nrows, ncols, figure=fig, wspace=0.6, hspace=0.4)

    for idx, block in enumerate(unique_blocks):
        print(f"Plotting locus: {block}")
        
        locus_data = plot_df[plot_df['locus_id'] == block].copy()
        
        chrom = locus_data['chr'].iloc[0]
        # Use min/max positions of the variants to set the window
        ld_start = locus_data['pos'].min()
        ld_end = locus_data['pos'].max()
        
        plot_start = ld_start - (buffer_kb * 1000)
        plot_end = ld_end + (buffer_kb * 1000)
        
        window_df = plot_df[
            (plot_df['chr'] == chrom) & 
            (plot_df['pos'] >= plot_start) & 
            (plot_df['pos'] <= plot_end)
        ].copy()

        window_genes = genes_df[
            (genes_df['chr'] == chrom) & 
            (genes_df['end'] >= plot_start) & 
            (genes_df['start'] <= plot_end)
        ].copy()

        row_idx = idx // ncols
        col_idx = idx % ncols

        # --- THE 3-PANEL GRIDSPEC ---
        gs_inner = GridSpecFromSubplotSpec(3, 1, subplot_spec=master_gs[row_idx, col_idx], 
                                           height_ratios=[3, 2, 1], hspace=0.1)
        
        ax_gwas = fig.add_subplot(gs_inner[0])
        ax_pip  = fig.add_subplot(gs_inner[1], sharex=ax_gwas)
        ax_gene = fig.add_subplot(gs_inner[2], sharex=ax_gwas)

        # 1. Define LocusZoom Color Mapping for GWAS panel based on SuSiE lead_r2
        colors_gwas = []
        for r2 in window_df['lead_r2']:
            if pd.isna(r2): colors_gwas.append('#D3D3D3')      # Grey (Background)
            elif r2 >= 0.8: colors_gwas.append('#D43F3A')      # Red
            elif r2 >= 0.6: colors_gwas.append('#EEA236')      # Orange
            elif r2 >= 0.4: colors_gwas.append('#5CB85C')      # Green
            elif r2 >= 0.2: colors_gwas.append('#5BC0DE')      # Light Blue
            else: colors_gwas.append('#357EBD')                # Dark Blue

        window_df['color_gwas'] = colors_gwas
        
        # Sort so high LD SNPs are drawn on top
        window_df = window_df.sort_values('lead_r2', ascending=True, na_position='first')

        # --- PANEL 1: GWAS ---
        # Plot Background
        bg_mask = window_df['lead_r2'].isna() | (window_df['lead_r2'] < 0.2)
        ax_gwas.scatter(
            window_df.loc[bg_mask, 'pos'] / 1e6, 
            window_df.loc[bg_mask, '-log10_plotting'],
            c=window_df.loc[bg_mask, 'color_gwas'],
            edgecolors='none', s=20, alpha=0.5, zorder=2
        )

        # Plot Foreground
        fg_mask = window_df['lead_r2'] >= 0.2
        ax_gwas.scatter(
            window_df.loc[fg_mask, 'pos'] / 1e6, 
            window_df.loc[fg_mask, '-log10_plotting'],
            c=window_df.loc[fg_mask, 'color_gwas'],
            edgecolors='white', linewidths=0.2, s=30, alpha=0.9, zorder=3
        )

        ax_gwas.set_ylabel(r'$-\log_{10}(P)$', fontweight='bold')
        ax_gwas.axhline(-np.log10((5e-8)/1060), color='black', linestyle='--', linewidth=0.8, zorder=1)
        plt.setp(ax_gwas.get_xticklabels(), visible=False)


        # --- PANEL 2: PIPs and Credible Sets ---
        # Plot Background SNPs (Not in a CS)
        bg_cs_mask = (window_df['CS'] == 0) | (window_df['CS'].isna())
        ax_pip.scatter(
            window_df.loc[bg_cs_mask, 'pos'] / 1e6, 
            window_df.loc[bg_cs_mask, 'PIP'],
            c='#D3D3D3', edgecolors='none', s=20, alpha=0.5, zorder=2
        )

        # Assign high-contrast colors to each Credible Set
        cs_colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00']
        unique_cs = sorted([cs for cs in window_df['CS'].unique() if pd.notna(cs) and cs > 0])
        
        lead_snps = [] # Keep track of the highest PIP in each CS
        
        for i, cs in enumerate(unique_cs):
            cs_mask = window_df['CS'] == cs
            color = cs_colors[i % len(cs_colors)]
            
            # Plot the CS SNPs
            ax_pip.scatter(
                window_df.loc[cs_mask, 'pos'] / 1e6, 
                window_df.loc[cs_mask, 'PIP'],
                c=color, edgecolors='black', linewidths=0.3, s=40, alpha=0.9, zorder=3,
                label=f'CS {int(cs)}'
            )
            
            # Identify Lead SNP for this CS
            lead_idx = window_df.loc[cs_mask, 'PIP'].idxmax()
            lead_snps.append(lead_idx)

        # Highlight Lead SNPs in BOTH top and middle panels
        for lead_idx in lead_snps:
            # Top Panel Diamond
            ax_gwas.scatter(
                window_df.loc[lead_idx, 'pos'] / 1e6, 
                window_df.loc[lead_idx, '-log10_plotting'],
                c='#9400D3', marker='D', edgecolors='black', linewidths=0.8, s=40, zorder=5
            )
            ax_gwas.annotate(
                window_df.loc[lead_idx, 'rsid'],
                xy=(window_df.loc[lead_idx, 'pos'] / 1e6, window_df.loc[lead_idx, '-log10_plotting']),
                xytext=(4, 4), textcoords="offset points", fontsize=6, fontstyle='italic', zorder=6
            )
            
            # Middle Panel Diamond
            ax_pip.scatter(
                window_df.loc[lead_idx, 'pos'] / 1e6, 
                window_df.loc[lead_idx, 'PIP'],
                c='#9400D3', marker='D', edgecolors='black', linewidths=0.8, s=40, zorder=5
            )

        ax_pip.set_ylabel('PIP', fontweight='bold')
        ax_pip.set_ylim(-0.05, 1.05)
        if len(unique_cs) > 0:
            ax_pip.legend(loc='upper right', frameon=False, fontsize=6)
        plt.setp(ax_pip.get_xticklabels(), visible=False)


        # --- PANEL 3: GENES ---
        y_levels = [0, 1, 2, 3]
        current_level = 0
        min_gene_size_bp = 10000 
        filtered_genes = window_genes[(window_genes['end'] - window_genes['start']) > min_gene_size_bp]

        for _, gene in filtered_genes.iterrows():
            g_start = max(gene['start'], plot_start) / 1e6
            g_end = min(gene['end'], plot_end) / 1e6
            g_center = (g_start + g_end) / 2
            
            ax_gene.plot([g_start, g_end], [current_level, current_level], 
                         color='#2C308B', linewidth=2, solid_capstyle='butt')
            ax_gene.text(g_center, current_level + 0.2, gene['gene_name'], 
                         fontsize=5, fontstyle='italic', ha='center', va='bottom', color='black')
            current_level = y_levels[(y_levels.index(current_level) + 1) % len(y_levels)]

        ax_gene.set_xlabel(f'Chromosome {chrom} Position (Mb)', fontweight='bold')
        ax_gene.set_xlim(plot_start / 1e6, plot_end / 1e6)
        ax_gene.set_ylim(-1, 3.5)
        ax_gene.set_yticks([])
        ax_gene.spines['left'].set_visible(False)
        ax_gene.spines['right'].set_visible(False)
        ax_gene.spines['top'].set_visible(False)

        ax_gwas.ticklabel_format(useOffset=False, style='plain', axis='x')
        ax_pip.ticklabel_format(useOffset=False, style='plain', axis='x')
        ax_gene.ticklabel_format(useOffset=False, style='plain', axis='x')

        # Custom Legend for GWAS Panel
        ld_legend = [
            mpatches.Patch(color='#D43F3A', label=r'$r^2 \geq 0.8$'),
            mpatches.Patch(color='#EEA236', label=r'$0.6 \leq r^2 < 0.8$'),
            mpatches.Patch(color='#5CB85C', label=r'$0.4 \leq r^2 < 0.6$'),
            mpatches.Patch(color='#5BC0DE', label=r'$0.2 \leq r^2 < 0.4$'),
            mpatches.Patch(color='#357EBD', label=r'$0.0 \leq r^2 < 0.2$'),
            plt.Line2D([0], [0], color='black', linestyle='--', linewidth=0.8, 
                       label=r'$P \approx 4.7 \times 10^{-11}$')
        ]
        ax_gwas.legend(handles=ld_legend, loc='upper left', bbox_to_anchor=(1.02, 1),
                       frameon=False, title="LD to Lead PIP SNP", fontsize=6, title_fontsize=7)

    fig.suptitle(f"Phenotype: {phen_name}", x=0.05, y=0.98, ha='left', fontsize=14, fontweight='bold')
    plt.subplots_adjust(top=0.92)

    output_file = os.path.join(output_dir, f'susie_plot_{phen_code}.pdf')
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.close()
    
    print(f"-> Saved grid plot: {output_file}")


if __name__ == "__main__":
    
    # Update this to where you save your master TSVs
    susie_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/susie_results" 
    phen_dict_path = "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/phen_dict_renamed.xlsx"
    gene_bed_path = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/ld_info/human_genes_hg19.bed"
    out_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/susie_plots"

    traits = pd.read_excel(phen_dict_path, usecols=["phenotype_code", "description"])

    for index, row in traits.iterrows():
        phen_code = row['phenotype_code']
        phen_name = row['description']

        print(f"Processing: {phen_code}")
        
        plot_df = get_susie_data(phen_code, susie_dir)

        if plot_df is None or plot_df.empty:
            print(f"No SuSiE data for {phen_code}. Skipping...")
            continue

        genes_df = get_gene_annotations(gene_bed_path)
        
        plot_regional_association(plot_df, genes_df, phen_code, phen_name, out_dir)