import os
import sys
import colorsys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec


def get_data(phen_code, sumstat_dir, loci_dir):
    """Gets and merges locus data."""

    raw_df = pd.read_csv(f"{sumstat_dir}/{phen_code}_sig_SNPs.tsv.bgz", sep='\t', compression='gzip')


    loci_df = pd.read_csv(f"{loci_dir}/{phen_code}_sig_loci.tsv", sep='\t', 
                          usecols=["variant", "neg_log10_pval","indep_status", "indep_id", "r4", 
                                   "lead_status", "lead_id","r4_lead",  "indep_start", 
                                   "indep_end", "ld_start", "ld_end", "ld_id"])
    
    if loci_df.empty:
        print("No loci to plot.")
        sys.exit(0)
    
    plot_df = raw_df.merge(loci_df, on="variant", how="left")

    return plot_df


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


def lighten_color(hex_color, amount=0.5):
    try:
        c = mcolors.to_rgb(hex_color)
        # Convert RGB to HLS (Hue, Lightness, Saturation)
        h, l, s = colorsys.rgb_to_hls(*c)
        # Increase lightness: l + (1.0 - l) * amount
        return colorsys.hls_to_rgb(h, l + (1.0 - l) * amount, s)
    except:
        return hex_color


def plot_regional_association(plot_df, genes_df, phen_code, phen_name, output_dir, buffer_kb=50):
    """Generates LocusZoom-style plots for each unique independent locus."""
    
    set_style()
    os.makedirs(output_dir, exist_ok=True)

    # Identify the unique ld blocks
    unique_blocks = plot_df['ld_id'].dropna().unique()
    n_blocks = len(unique_blocks)
    print(f"Found {len(unique_blocks)} unique loci to plot.")

    plot_df = plot_df.sort_values("neg_log10_pval", ascending=False).reset_index(drop=True)
    
    lowest_valid_p = plot_df.loc[plot_df["dominance_pval"] != 0, "neg_log10_pval"].max()
    plot_df["-log10_plotting"] = np.where(
        plot_df["dominance_pval"] == 0, 
        lowest_valid_p + (plot_df["neg_log10_pval"] / 100),
        -np.log10(plot_df["dominance_pval"])
    )

    ncols = 2  # Set to 1 for a single vertical column, 2 or 3 for a grid
    nrows = int(np.ceil(n_blocks / ncols))

    # Scale the figure size dynamically. 
    # Width = 3.5 per col (+ padding for legends). Height = 4.5 per row.
    fig = plt.figure(figsize=(6.5 * ncols, 4.5 * nrows))
    
    # Create the master grid. wspace=0.7 guarantees your custom right-side legends won't overlap the next column!
    master_gs = GridSpec(nrows, ncols, figure=fig, wspace=0.6, hspace=0.4)

    # 3. Loop through each master locus and generate a subplot
    for idx, block in enumerate(unique_blocks):
        print(f"Plotting locus: {block}")
        
        locus_data = plot_df[plot_df['ld_id'] == block].copy()
        
        chrom = locus_data['chr'].iloc[0]
        ld_start = locus_data['ld_start'].min()
        ld_end = locus_data['ld_end'].max()
        
        # Define the plot window (add the buffer)
        plot_start = ld_start - (buffer_kb * 1000)
        plot_end = ld_end + (buffer_kb * 1000)
        
        # Slice the specific window for plotting
        window_df = plot_df[
            (plot_df['chr'] == chrom) & 
            (plot_df['pos'] >= plot_start) & 
            (plot_df['pos'] <= plot_end)
        ].copy()
        

        # SNPs not in locus_data will get NaN, which we fill with -1 (for background grey)
        window_df['r4'] = window_df['r4'].fillna(-1)

        window_genes = genes_df[
            (genes_df['chr'] == chrom) & 
            (genes_df['end'] >= plot_start) & 
            (genes_df['start'] <= plot_end)
        ].copy()

        row_idx = idx // ncols
        col_idx = idx % ncols

        # Create your 3:1 nested GridSpec exactly inside this specific grid cell
        gs_inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=master_gs[row_idx, col_idx], 
                                           height_ratios=[3, 1], hspace=0.05)
        
        ax_main = fig.add_subplot(gs_inner[0])
        ax_gene = fig.add_subplot(gs_inner[1], sharex=ax_main)


        # Define standard LocusZoom Color Mapping
        colors = []
        for r4 in window_df['r4']:
            if r4 >= 0.8: colors.append('#D43F3A')      # Red
            elif r4 >= 0.6: colors.append('#EEA236')    # Orange
            elif r4 >= 0.4: colors.append('#5CB85C')    # Green
            elif r4 >= 0.2: colors.append('#5BC0DE')    # Light Blue
            elif r4 >= 0.0: colors.append("#357EBD")    # Dark Blue
            else: colors.append('#357EBD')              # Grey (Background)

        window_df['color'] = colors
        
        # Sort so the highest correlation SNPs are drawn ON TOP of the grey background
        window_df = window_df.sort_values('r4')
        window_df = window_df.drop_duplicates(subset=['variant'], keep='first')

        # Plot Background and Proxies
        scatter_df = window_df[(window_df['indep_status'] != True) & (window_df['lead_status'] != True)]
        ax_main.scatter(
            scatter_df['pos'] / 1e6, # Convert to Megabases for clean X-axis labels
            scatter_df['-log10_plotting'],
            c=scatter_df['color'],
            edgecolors='white',
            linewidths=0.1,
            s=20,
            alpha=0.9,
            zorder=2
        )

         # Highlight the Lead SNPs
        lead_row = window_df[window_df['lead_status'] == True]
        if not lead_row.empty:
            ax_main.scatter(
                lead_row['pos'] / 1e6,
                lead_row['-log10_plotting'],
                c='#9400D3',
                marker='D',          
                edgecolors='black',
                linewidths=0.8,
                s=25,
                zorder=3
            )
        
        # Annotate the Lead SNPs with their variant ID
        for _, row in lead_row.iterrows():
            ax_main.annotate(
                row['rsid'],
                xy=(row['pos'] / 1e6, row['-log10_plotting']),
                xytext=(4, 4),   # Offset the text 4 points right and 4 points up
                textcoords="offset points",
                fontsize=6,
                fontstyle='italic',
                color='black',
                zorder=4         # Keep text on the very top layer
            )

        
        # Highlight the independent SNPs
        indep_row = window_df[(window_df['indep_status'] == True) & (window_df['lead_status'] != True)]

        colors_indep = []
        lead_color = '#9400D3'
        for r4_l in indep_row['r4_lead']:
            if r4_l >= 0.243: colors_indep.append(lighten_color(lead_color, amount=0.2))    
            elif 0.127 <= r4_l < 0.243: colors_indep.append(lighten_color(lead_color, amount=0.5))
            elif 0.01 <= r4_l < 0.127: colors_indep.append(lighten_color(lead_color, amount=0.8))

        indep_row['color_indep'] = colors_indep

        if not indep_row.empty:
            ax_main.scatter(
                indep_row['pos'] / 1e6,
                indep_row['-log10_plotting'],
                c=indep_row['color_indep'],      
                marker='^',          
                edgecolors='black',
                linewidths=0.8,
                s=25,
                zorder=3
            )

        # 8. Formatting the Axes
        ax_main.set_ylabel(r'$-\log_{10}(P)$', fontweight='bold')
        ax_main.set_xlim(plot_start / 1e6, plot_end / 1e6)
        
        # Add the genome-wide significance line
        ax_main.axhline(-np.log10((5e-8)/1060), color='black', linestyle='--', linewidth=0.8, zorder=1)

        # 9. Custom LD Legend
        ld_legend = [
            mpatches.Patch(color='#D43F3A', label=r'$r^4 \geq 0.8$'),
            mpatches.Patch(color='#EEA236', label=r'$0.6 \leq r^4 < 0.8$'),
            mpatches.Patch(color='#5CB85C', label=r'$0.4 \leq r^4 < 0.6$'),
            mpatches.Patch(color='#5BC0DE', label=r'$0.2 \leq r^4 < 0.4$'),
            mpatches.Patch(color='#357EBD', label=r'$0.0 \leq r^4 < 0.2$'),
            mpatches.Patch(color='none', label=' '),mpatches.Patch(color='none', label=' '),
            mpatches.Patch(color='none', label=' '),mpatches.Patch(color='none', label=' '),
            mpatches.Patch(color='none', label=' '),mpatches.Patch(color='none', label=' '),
            mpatches.Patch(color='none', label=' '),
            plt.Line2D([0], [0], color='black', linestyle='--', linewidth=0.8, 
                       label=r'$P = 4.7 \times 10^{-11}$')
        ]

        if indep_row.empty:
            ld_title = "LD to Lead SNP"

        elif not indep_row.empty:
            ld_title = 'LD to Independent SNP'

            # Lead SNP color legend
            cmap_colors = [
                lighten_color(lead_color, amount=0.8), 
                lighten_color(lead_color, amount=0.5), 
                lighten_color(lead_color, amount=0.2)
                ]
            custom_cmap = mcolors.LinearSegmentedColormap.from_list("indep_purple", cmap_colors)

            cax = ax_main.inset_axes([1.089, 0.39, 0.0475, 0.029]) 
            
            norm = mcolors.Normalize(vmin=0.01, vmax=0.30) 
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=custom_cmap),
                            cax=cax, orientation='horizontal')
            cb.set_ticks([])
            cb.outline.set_linewidth(0)
            cax.text(1.3, 0.5, r'$0.01 \leq r^4 < 0.36$', size=6, 
                    transform=cax.transAxes, ha='left', va='center')
            cax.text(-0.42, 1.8, 'LD to Lead SNP', fontsize=7, 
                    transform=cax.transAxes, ha='left', va='bottom')


        ax_main.legend( handles=ld_legend,loc='upper left',
                bbox_to_anchor=(1.05, 1),frameon=False,
                title=ld_title,
                fontsize=6,
                title_fontsize=7
                )

        plt.setp(ax_main.get_xticklabels(), visible=False)
        y_levels = [0, 1, 2, 3]
        current_level = 0


        # Define minimum size threshold in base pairs
        min_gene_size_bp = 10000 
        
        # Filter the dataframe to ONLY include genes larger than the threshold
        filtered_genes = window_genes[(window_genes['end'] - window_genes['start']) > min_gene_size_bp]

        for _, gene in filtered_genes.iterrows():
            g_start = max(gene['start'], plot_start) / 1e6
            g_end = min(gene['end'], plot_end) / 1e6
            g_center = (g_start + g_end) / 2
            
            # Draw the line and the text for the filtered genes
            ax_gene.plot([g_start, g_end], [current_level, current_level], 
                         color='#2C308B', linewidth=2, solid_capstyle='butt')
            ax_gene.text(g_center, current_level + 0.2, gene['gene_name'], 
                         fontsize=4, fontstyle='italic', ha='center', va='bottom', color='black')
            
            # Increment the level
            current_level = y_levels[(y_levels.index(current_level) + 1) % len(y_levels)]


        ax_gene.set_xlabel(f'Chromosome {chrom} Position (Mb)', fontweight='bold')
        ax_gene.set_xlim(plot_start / 1e6, plot_end / 1e6)
        ax_gene.set_ylim(-1, 3.5)
        ax_gene.set_yticks([])
        ax_gene.spines['left'].set_visible(False)
        ax_gene.spines['right'].set_visible(False)
        ax_gene.spines['top'].set_visible(False)


    # 10. Save as vector PDF
    fig.suptitle(f"Phenotype: {phen_name}", x=0.05, y=0.98, ha='left', fontsize=14, fontweight='bold')
    
    # Push the entire grid down slightly to make room for the header
    plt.subplots_adjust(top=0.92)

    output_file = os.path.join(output_dir, f'l_plot_{phen_code}.pdf')

    # Notice tight_layout() is removed! bbox_inches='tight' handles it perfectly.
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.close()
    
    print(f"-> Saved  grid plot: {output_file}")


if __name__ == "__main__":
    
    sumstat_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sumstats_QCed"
    phen_dict_path = "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/phen_dict_renamed.xlsx"
    gene_bed_path = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/ld_info/human_genes_hg19.bed"
    loci_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sig_loci"
    out_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/loci_plots"

    traits = pd.read_excel(phen_dict_path, usecols=["phenotype_code", "description"])
    #traits = traits[traits["phenotype_code"] == "1747_2"]



    for index, row in traits.iterrows():
        
        phen_code = row['phenotype_code']
        phen_name = row['description']

        try:
            plot_df = get_data(phen_code, sumstat_dir, loci_dir)
        except FileNotFoundError:
            # Skip to the next phenotype if the file/dir doesn't exist
            print(f"No significant loci data for{phen_code} ")
            continue
        
        print(f"Processing: {phen_code}")

        genes_df = get_gene_annotations(gene_bed_path)
        
        plot_regional_association(plot_df, genes_df, phen_code, phen_name, out_dir)