import os
import colorsys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerTuple
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec



def get_data(phen_code, sumstat_dir, loci_dir, susie_dir):
    """Gets and merges locus data."""

    raw_df = pd.read_csv(f"{sumstat_dir}/{phen_code}_sig_SNPs.tsv.bgz", sep='\t', compression='gzip')

    # Create a mask for  MHC region (Chr 6, 25 Mb - 34 Mb)
    mhc_mask = (
        (raw_df['chr'] == 6) & 
        (raw_df['pos'] >= 25000000) & 
        (raw_df['pos'] <= 34000000)
    )

    # Keep everything that is NOT in the mask
    masked_df = raw_df[~mhc_mask].reset_index(drop=True).copy()

    if masked_df.empty:
        return None

    loci_df = pd.read_csv(f"{loci_dir}/{phen_code}_sig_loci.tsv", sep='\t', 
                          usecols=["variant","indep_status", "indep_id", "r4", 
                                   "lead_status", "lead_id","r4_lead",  "indep_start", 
                                   "indep_end", "ld_start", "ld_end", "ld_id"])
    
    if loci_df.empty:
        print("No loci to plot.")
        return None
    
    susie_df = pd.read_csv(f"{susie_dir}/{phen_code}_susie_res.tsv", sep='\t', 
                        usecols=["variant", "locus_id", "pos", "rsid", "dom_z_score", "PIP", "CS", "CS_prob", "low_purity", 
                                 "lead_r2", "post_mean", "post_sd", "lambda"])
    
    if susie_df.empty:
        print("No SuSiE results to plot.")
        return None
    

    plot_df = masked_df.merge(loci_df, on="variant", how="left")
    
    # Remove duplicates
    indep_pvals = masked_df[['variant', 'dom_log10_pval']].rename(
        columns={'variant': 'indep_id', 'dom_log10_pval': 'indep_pval'}
    )
    
    plot_df = plot_df.merge(indep_pvals, on='indep_id', how='left')
    plot_df = plot_df[plot_df['indep_pval'].isna() | (plot_df['indep_pval'] >= plot_df['dom_log10_pval'])]

    group_sizes = plot_df.groupby(['variant', 'lead_id'])['variant'].transform('size')
    rows_to_drop = (group_sizes > 1) & (plot_df['lead_id'] == plot_df['indep_id'])

    plot_df = plot_df[~rows_to_drop]
    plot_df = plot_df.drop(columns=['indep_pval'])

    return plot_df, susie_df


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


def plot_regional_association(plot_number, plot_df, genes_df, susie_df, phen_code, phen_name, output_dir, buffer_kb=50):
    """Generates LocusZoom-style plots for each unique independent locus."""
    
    set_style()
    os.makedirs(output_dir, exist_ok=True)

    # Identify the unique ld blocks
    unique_blocks = plot_df['ld_id'].dropna().unique()
    n_blocks = len(unique_blocks)
    print(f"Found {len(unique_blocks)} unique loci to plot.")

    if n_blocks == 0:
        print(f"No locus to plot for ' {phen_name}'. Skipping...")
        return None
    else:
        plot_number += 1

    plot_df = plot_df.sort_values("dom_log10_pval", ascending=False).reset_index(drop=True)
    
    plot_df["-log10_plotting"] = np.where(
        plot_df["dom_log10_pval"] > 300, 
        300 + (plot_df["dom_log10_pval"] / 100),
        plot_df["dom_log10_pval"]
        )

    ncols = 2  # Set to 1 for a single vertical column, 2 or 3 for a grid
    nrows = int(np.ceil(n_blocks / ncols))

    # Scale the figure size
    fig = plt.figure(figsize=(6.5 * ncols, 6.5 * nrows))
    
    # Create the master grid.
    master_gs = GridSpec(nrows, ncols, figure=fig, wspace=0.6, hspace=0.4)

    # Loop through locus and generate a subplot
    for idx, block in enumerate(unique_blocks):
        print(f"Plotting locus: {block}")
        
        locus_data = plot_df[plot_df['ld_id'] == block].copy()
        pip_df = susie_df[susie_df["locus_id"] == block].copy()
        
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
        

        # SNPs not in locus_data will get NaN
        window_df['r4'] = window_df['r4'].fillna(-1)

        window_genes = genes_df[
            (genes_df['chr'] == chrom) & 
            (genes_df['end'] >= plot_start) & 
            (genes_df['start'] <= plot_end)
        ].copy()

        row_idx = idx // ncols
        col_idx = idx % ncols

        # 3:1 nested GridSpec
        gs_inner = GridSpecFromSubplotSpec(3, 1, subplot_spec=master_gs[row_idx, col_idx], 
                                    height_ratios=[3, 1.5, 0.75], hspace=0.2)
        
        ax_main = fig.add_subplot(gs_inner[0])
        ax_pip  = fig.add_subplot(gs_inner[1], sharex=ax_main)
        ax_gene = fig.add_subplot(gs_inner[2], sharex=ax_main)


        # Define Color Mapping
        colors = []
        for r4 in window_df['r4']:
            if r4 >= 0.9: colors.append('#D43F3A')      # Red
            elif r4 >= 0.8: colors.append('#EEA236')    # Orange
            elif r4 >= 0.7: colors.append('#5CB85C')    # Green
            elif r4 >= 0.6: colors.append('#5BC0DE')    # Light Blue
            elif r4 >= 0.0: colors.append("#357EBD")    # Dark Blue
            else: colors.append('#357EBD')              # Grey (Background)

        window_df['color'] = colors
        
        # 
        window_df = window_df.sort_values('r4', ascending=False)
        window_df = window_df.drop_duplicates(subset=['variant'], keep='first')

        # Plot Background and Proxies
        scatter_df = window_df[(window_df['indep_status'] != True) & (window_df['lead_status'] != True)]
        ax_main.scatter(
            scatter_df['pos'] / 1e6, # Convert to Megabases
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
                zorder=4         
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
        
        if window_df['-log10_plotting'].max() > 300:
            
            # Find the exact row with the highest value
            highest_snp = window_df.loc[window_df['dom_log10_pval'].idxmax()]
            
            ax_main.scatter(
                            highest_snp['pos'] / 1e6,
                            highest_snp['-log10_plotting'] + 20,
                            marker='*',
                            color='#780505',
                            s=15,           
                            zorder=5
                        )

        # Formatting the Axes
        ax_main.set_ylabel(r'$-\log_{10}(P)$', fontweight='bold')
        ax_main.set_xlim(plot_start / 1e6, plot_end / 1e6)
        
        # Add the genome-wide significance line
        ax_main.axhline(-np.log10((5e-8)/1060), color='black', linestyle='--', linewidth=0.8, zorder=1)

        # LD Legend
        lead_color = '#9400D3'
        
        # Base Legend (Colors + Lead SNP)
        ld_legend = [
            mpatches.Patch(color='#D43F3A', label=r'$r_D^2 \geq 0.9$'),
            mpatches.Patch(color='#EEA236', label=r'$0.8 \leq r_D^2 < 0.9$'),
            mpatches.Patch(color='#5CB85C', label=r'$0.7 \leq r_D^2 < 0.8$'),
            mpatches.Patch(color='#5BC0DE', label=r'$0.6 \leq r_D^2 < 0.7$'),
            mpatches.Patch(color='#357EBD', label=r'$0.0 \leq r_D^2 < 0.6$'),
            
            mpatches.Patch(color='none', label=' '),
            plt.Line2D([0], [0], marker='D', color='w', markerfacecolor=lead_color, 
                       markeredgecolor='black', markersize=5, label='Lead SNP')

        ]

        # Independent SNPs (Triangles) only if they exist in this plot
        if not indep_row.empty:
            ld_legend.extend([
                plt.Line2D([0], [0], marker='^', color='w', markerfacecolor=lead_color, 
                       markeredgecolor='black', markersize=5, label='Independent SNP'),

                mpatches.Patch(color='none', label=' '),

                mpatches.Patch(color='none', label="LD to Lead SNP")
            ])
            
            spacer = plt.Line2D([0], [0], marker='none', color='none')
            
            t1 = plt.Line2D([0], [0], marker='^', linestyle='none', 
                            markerfacecolor=lighten_color(lead_color, amount=0.8), markeredgecolor='black', markersize=6)
            t2 = plt.Line2D([0], [0], marker='^', linestyle='none', 
                            markerfacecolor=lighten_color(lead_color, amount=0.5), markeredgecolor='black', markersize=6)
            t3 = plt.Line2D([0], [0], marker='^', linestyle='none', 
                            markerfacecolor=lighten_color(lead_color, amount=0.2), markeredgecolor='black', markersize=6)
            
            ld_legend.extend([
                (spacer, spacer,spacer, spacer, t1, t2, t3) 
            ])

            ld_legend.extend([
            mpatches.Patch(color='none', label=r'$0.1 \leq r_D^2 < 0.6$')
        ])

        # Append the Threshold Line
        ld_legend.extend([
            mpatches.Patch(color='none', label=' '), 
            plt.Line2D([0], [0], color='black', linestyle='--', linewidth=0.8, 
                       label=r'$P \approx 4.7 \times 10^{-11}$'),
        ])
        
        # Append the Star (if compressed)
        max_pval = window_df['dom_log10_pval'].max()
        if max_pval > 300:
            ld_legend.append(
                plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='#780505', 
                           markeredgecolor='#780505', markersize=6, 
                           label=f'Max $-\\log_{{10}}(P) = {max_pval:.0f}$')
            )

        # Draw the Legend
        ld_labels = [h.get_label() if not isinstance(h, tuple) else h[0].get_label() for h in ld_legend]

        ax_main.legend(
            handles=ld_legend,
            labels=ld_labels,       
            loc='upper left',
            bbox_to_anchor=(1.05, 1), 
            frameon=False,
            title='LD to Independent SNP' if not indep_row.empty else 'LD to Lead SNP',
            fontsize=6,
            title_fontsize=7,
            handlelength=1.2,
            labelspacing=0.5,
            handler_map={tuple: HandlerTuple(ndivide=None, pad=16)} 
        )


        plt.setp(ax_main.get_xticklabels(), visible=False)
        y_levels = [0, 1, 2, 3]
        current_level = 0

        
        # ---------------- PANEL 2: PIPs and Credible Sets ------------------ #
        pip_df['abs_z'] = pip_df['dom_z_score'].abs()
        pip_df = pip_df.sort_values('abs_z', ascending=False)

        # Background SNPs (Not in a CS)
        bg_cs_mask = (pip_df['CS'] == 0) | (pip_df['CS'].isna())
        ax_pip.scatter(
            pip_df.loc[bg_cs_mask, 'pos'] / 1e6,  # Still using the jittered X
            pip_df.loc[bg_cs_mask, 'PIP'],
            c='#357EBD', edgecolors='white', linewidths=0.1, s=20, alpha=0.9, zorder=2
        )

        # colors to each Credible Set
        cs_colors = [
            '#E41A1C',
            '#377EB8',
            '#4DAF4A',
            '#984EA3',
            '#FF7F00',
            '#A65628',
            '#F781BF',
            '#17BECF', 
            '#800000',
            '#B8860B'
        ]
        pip_df =pip_df[pip_df["low_purity"] == False].copy()
        unique_cs = sorted([cs for cs in pip_df['CS'].unique() if pd.notna(cs) and cs > 0])
        
        pip_legend = []

        for i, cs in enumerate(unique_cs):
            cs_mask = pip_df['CS'] == cs
            color = cs_colors[i % len(cs_colors)]
            
            # Plot the CS SNPs
            ax_pip.scatter(
                pip_df.loc[cs_mask, 'pos'] / 1e6, 
                pip_df.loc[cs_mask, 'PIP'],
                c=color, edgecolors='black', linewidths=1, s=25, alpha=0.9, zorder=3,
                label=f'CS {int(cs)}'
            )

            # --- legend ---
            max_pip = pip_df.loc[cs_mask, 'PIP'].max()
            top_snps = pip_df.loc[cs_mask & (pip_df['PIP'] == max_pip), 'rsid'].tolist()
            rsid_str = ", ".join(top_snps)

            if max_pip >= 0.01:
                pip_str = f"{max_pip:.2f}"    
            elif max_pip >= 0.001:
                pip_str = f"{max_pip:.3f}"  
            elif max_pip >= 0.0001:
                pip_str = f"{max_pip:.4f}"  
            else:
                pip_str = f"{max_pip:.2e}"
            
            pip_legend.append(
                plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, 
                           markeredgecolor='black', markersize=5, 
                           label=f'CS {int(cs)}: {rsid_str} (PIP = {pip_str})')
            )

        if pip_legend:
            ax_pip.legend(
                handles=pip_legend,
                loc='upper left',
                bbox_to_anchor=(1.05, 1), 
                frameon=False,
                title='95% Credible Sets\n(Highest PIP SNP)',
                fontsize=5,
                title_fontsize=7,
                labelspacing=0.6,
                handletextpad=0.2
            )
            
        
        ax_pip.set_ylabel('PIP')

         # ----------------------------------------------------------- #


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

        # Turn off the automatic scientific offset
        ax_main.ticklabel_format(useOffset=False, style='plain', axis='x')
        ax_gene.ticklabel_format(useOffset=False, style='plain', axis='x')


    # Save
    fig.text(0.05, 0.98, f"Figure {plot_number}", ha='left', fontsize=14, fontweight='bold')
    fig.text(0.05, 0.95, f"Phenotype: {phen_name}", ha='left', fontsize=11, fontweight='bold')
    
    plt.subplots_adjust(top=0.92)

    output_file = os.path.join(output_dir, f'{plot_number}_plot_{phen_code}.pdf')

    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.close()
    
    print(f"-> Saved  grid plot: {output_file}")
    
    return plot_number


if __name__ == "__main__":
    
    sumstat_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sumstats_QCed"
    phen_dict_path = "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/phen_dict_renamed.xlsx"
    gene_bed_path = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/ld_info/human_genes_hg19.bed"
    loci_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sig_loci"
    susie_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/susie_results"
    out_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/loci_results/loci_plots"

    traits = pd.read_excel(phen_dict_path, usecols=["phenotype_code", "description"])
    #traits = traits[traits["phenotype_code"] == "1747_2"]

    plot_number = 0
    for index, row in traits.iterrows():
        phen_code = row['phenotype_code']
        phen_name = row['description']

        try:
            plot_df, susie_df = get_data(phen_code, sumstat_dir, loci_dir, susie_dir)

            if plot_df is None or plot_df.empty:
                print(f"No data returned for {phen_code}")
                continue
                
        except FileNotFoundError:
            print(f"No significant loci data for{phen_code} ")
            continue
        
        print(f"Processing: {phen_code}")

        genes_df = get_gene_annotations(gene_bed_path)
        
        plot_number = plot_regional_association(plot_number, plot_df, genes_df, susie_df, phen_code, phen_name, out_dir)