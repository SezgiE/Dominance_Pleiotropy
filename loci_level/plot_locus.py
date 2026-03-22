import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def get_data(phen_code, sumstat_dir, loci_dir):
    """Gets and merges locus data."""

    raw_df = pd.read_csv(f"{sumstat_dir}/{phen_code}_sig_SNPs.tsv.bgz", sep='\t', compression='gzip')
    loci_df = pd.read_csv(f"{loci_dir}/{phen_code}_sig_loci.tsv", sep='\t', 
                          usecols=["variant", "indep_status", "indep_id", "r4", 
                                   "lead_status", "lead_id","r4_lead",  "indep_start", 
                                   "indep_end", "ld_start", "ld_end", "ld_id"])
    
    if loci_df.empty:
        print("No loci to plot.")
        sys.exit(0)
    
    plot_df = raw_df.merge(loci_df, on="variant", how="left")

    return plot_df


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


def plot_regional_association(plot_df, output_dir, buffer_kb=250):
    """Generates LocusZoom-style plots for each unique independent locus."""
    
    set_style()
    os.makedirs(output_dir, exist_ok=True)

    # Identify the unique ld blocks
    unique_blocks = plot_df['ld_id'].dropna().unique()
    print(f"Found {len(unique_blocks)} unique loci to plot.")

    plot_df["-log10P"] = -np.log10(plot_df['dominance_pval'])

    # 3. Loop through each master locus and generate a plot
    for block in unique_blocks:
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

        # 4. Initialize the Figure
        fig, ax = plt.subplots(figsize=(3.5, 3.5)) # 89mm width is standard 1-column
        
        # Define standard LocusZoom Color Mapping
        colors = []
        for r4 in window_df['r4']:
            if r4 >= 0.8: colors.append('#D43F3A')      # Red
            elif r4 >= 0.6: colors.append('#EEA236')    # Orange
            elif r4 >= 0.4: colors.append('#5CB85C')    # Green
            elif r4 >= 0.2: colors.append('#5BC0DE')    # Light Blue
            elif r4 >= 0.0: colors.append('#357EBD')    # Dark Blue
            else: colors.append('#B0B0B0')              # Grey (Background)

        window_df['color'] = colors
        
        # Sort so the highest correlation SNPs are drawn ON TOP of the grey background
        window_df = window_df.sort_values('r4')
        window_df = window_df.drop_duplicates(subset=['variant'], keep='first')

        # Plot Background and Proxies
        ax.scatter(
            window_df['pos'] / 1e6, # Convert to Megabases for clean X-axis labels
            window_df['-log10P'],
            c=window_df['color'],
            edgecolors='white',
            linewidths=0.2,
            s=25,
            alpha=0.9,
            zorder=2
        )

        # Highlight the Lead SNP
        lead_row = window_df[window_df['indep_status'] == True]
        if not lead_row.empty:
            ax.scatter(
                lead_row['pos'] / 1e6,
                lead_row['-log10P'],
                c='#9400D3',          # Purple
                marker='D',           # Diamond
                edgecolors='black',
                linewidths=0.8,
                s=60,
                zorder=3
            )

        # 8. Formatting the Axes
        ax.set_xlabel(f'Chromosome {chrom} Position (Mb)', fontweight='bold')
        ax.set_ylabel(r'$-\log_{10}(P)$', fontweight='bold')
        ax.set_xlim(plot_start / 1e6, plot_end / 1e6)
        
        # Add the genome-wide significance line
        ax.axhline(-np.log10((5e-8)/1060), color='black', linestyle='--', linewidth=0.8, zorder=1)

        # 9. Custom LD Legend
        legend_elements = [
            mpatches.Patch(color='#D43F3A', label=r'$r^4 \geq 0.8$'),
            mpatches.Patch(color='#EEA236', label=r'$0.6 \leq r^4 < 0.8$'),
            mpatches.Patch(color='#5CB85C', label=r'$0.4 \leq r^4 < 0.6$'),
            mpatches.Patch(color='#5BC0DE', label=r'$0.2 \leq r^4 < 0.4$'),
            mpatches.Patch(color='#357EBD', label=r'$0.0 \leq r^4 < 0.2$'),
        ]
        ax.legend(handles=legend_elements, loc='upper right', frameon=False, title='LD to Lead SNP', fontsize=6, title_fontsize=7)

        # 10. Save as vector PDF
        safe_lead_name = block.replace(':', '_')
        output_file = os.path.join(output_dir, f'LocusZoom_{safe_lead_name}.pdf')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"-> Saved: {output_file}")


if __name__ == "__main__":
    
    
    phen_code = "1747_1"
    sumstat_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sumstats_QCed"
    loci_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/sig_loci"
    out_dir = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/loci_plots"

    plot_df = get_data(phen_code, sumstat_dir, loci_dir)
    
    plot_regional_association(plot_df, out_dir)