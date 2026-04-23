import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def run_permutation_test(pip_vals_path, n_iterations=1000, metric=lambda x: np.percentile(x, 95), seed=42):

    pip_df = pd.read_parquet(pip_vals_path)
    
    # Extract the arrays from the DataFrame
    parent_pop = pip_df[pip_df["source"] == "all"]["PIP"].dropna().values
    observed_subset = pip_df[pip_df["source"] == "pleio"]["PIP"].dropna().values
    
    n_subset = len(observed_subset)
    
    # Calculate the observed metrics
    observed_stat = metric(observed_subset)
    parent_stat = metric(parent_pop)
    
    # Generate the empirical null distribution
    null_distribution = np.zeros(n_iterations)
    
    for i in range(n_iterations):
        print(f"Iteration: {i}")
        # Draw a random sample
        random_sample = np.random.choice(parent_pop, size=n_subset, replace=False)
        null_distribution[i] = metric(random_sample)
        
    # Calculate the empirical p-value
    if observed_stat < parent_stat:
        p_value = np.sum(null_distribution <= observed_stat) / n_iterations
    else:
        p_value = np.sum(null_distribution >= observed_stat) / n_iterations
        
    # Multiply by 2 for a two-tailed test
    p_value_2tailed = min(1.0, p_value * 2)
    
    # Print summary
    print("--- Resampling Test Results ---")
    print(f"Parent Population Size: {len(parent_pop)}")
    print(f"Subset Size: {n_subset}")
    print(f"Parent Metric (Median): {parent_stat:.4f}")
    print(f"Observed Subset Metric (Median): {observed_stat:.4f}")
    print(f"P-value (2-tailed): {p_value_2tailed:.4f}")
    
    return p_value_2tailed


def set_style():
    """Hardcodes matplotlib parameters to format the plot."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 12,
        'axes.labelsize': 12,
        'xtick.labelsize': 9,
        'ytick.labelsize': 8,
        'axes.linewidth': 1.0,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'pdf.fonttype': 42,
        'ps.fonttype': 42
    })


def gtex_plot(fisher_results_path, pip_vals_path, perm_pip_pval, output_dir):

    fisher_all_df = pd.read_csv(fisher_results_path, sep="\t")
    fisher_all_df['name'] = fisher_all_df['name'].str.replace('_', ' ')
    
    set_style()
    colors = ['#3C5488', '#E64B35', '#4DBBD5', '#00A087', '#7E818D']
    

    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.2, 1], height_ratios=[1, 1], wspace=0.2, hspace=0.2)
    
    ax_a = fig.add_subplot(gs[:, 0])    # Panel A: Spans both rows in the left column
    ax_b = fig.add_subplot(gs[0, 1])    # Panel B: Top right row
    ax_c = fig.add_subplot(gs[1, 1])    # Panel C: Bottom right row
    

    # ---------------------------------------------------------
    # Panel A: Tissue Enrichment (Horizontal)
    # ---------------------------------------------------------
    tissue_df = fisher_all_df[fisher_all_df['type'] == 'tissue'].copy()
    tissue_df = tissue_df.sort_values('pleio_variants', ascending=False).reset_index(drop=True)
    
    y_pos_a = np.arange(len(tissue_df))
    # Swapped to barh for horizontal bars
    ax_a.barh(y_pos_a, tissue_df['pleio_variants'], color=colors[0], height=0.7, edgecolor='black', linewidth=1.0)
    
    # Significance markers (adjusted X/Y coordinates for horizontal alignment)
    for i, row in tissue_df.iterrows():
        if row['is_significant']:
            ax_a.text(row['pleio_variants'] + (tissue_df['pleio_variants'].max() * 0.01), i+0.35, 
                      '*', ha='left', va='center', fontsize=12, color='black')
            
    ax_a.plot([], [], marker='', color='black', linestyle='None', markersize=10, label='* FDR-corrected p < 0.05')
    ax_a.legend(frameon=False, loc='lower right')
            
    # Swapped ticks and labels to the Y-axis
    ax_a.set_xlim(left=0)
    ax_a.set_ylim(-0.5, len(tissue_df) - 0.5) 
    ax_a.set_yticks(y_pos_a)
    ax_a.set_yticklabels(tissue_df['name'])
    ax_a.invert_yaxis() # Flips the axis so the highest value is at the top
    ax_a.set_xlabel('Pleiotropic Variants (N)')
    
    # Shifted the 'A' label further left to account for the long tissue names on the Y-axis
    ax_a.text(-0.45, 1.045, 'A.', transform=ax_a.transAxes, fontsize=16, fontweight='bold', va='top')


    # ---------------------------------------------------------
    # Panel B: Biotype Enrichment
    # ---------------------------------------------------------
    bio_df = fisher_all_df[fisher_all_df['type'] == 'biotype'].copy()
    bio_df['name'] = bio_df['name'].str.title()
    bio_df = bio_df.sort_values('pleio_variants', ascending=False).reset_index(drop=True)
    
    x_pos_b = np.arange(len(bio_df))
    
    biotype_colors = {
        'Protein Coding':'#00A087',         # Dark Blue
        'Lncrna': '#E64B35',                 # Red
        'Processed Pseudogene': '#4DBBD5',   # Light Blue
        'Unprocessed Pseudogene': '#3C5488'  # Teal
    }

    for i in range(len(bio_df)):
        biotype_name = bio_df['name'].iloc[i]
        bar_color = biotype_colors.get(biotype_name, '#7E818D')
        
        ax_b.bar(x_pos_b[i], bio_df['pleio_variants'].iloc[i], color=bar_color, 
                 width=0.7, edgecolor='black', linewidth=1.0, label=bio_df['name'].iloc[i])
    
    main_legend = ax_b.legend(frameon=True, edgecolor='black', fancybox=False, framealpha=1.0, loc='upper right', ncol=1, fontsize=10)
    ax_b.add_artist(main_legend)
    
    star_line, = ax_b.plot([], [], marker='', color='black', linestyle='None', markersize=10)
    ax_b.legend([star_line], ['* FDR-corrected p < 0.05'], frameon=False, 
                loc='upper right', bbox_to_anchor=(1.0, 0.65), fontsize=9)

    ax_b.set_xticks(x_pos_b)
    ax_b.set_xticklabels("", rotation=0, ha="center")
    ax_b.set_ylabel('Pleiotropic Variants (N)')
    ax_b.set_xlabel('Gene Biotypes')
    ax_b.text(-0.15, 1.1, 'B.', transform=ax_b.transAxes, fontsize=16, fontweight='bold', va='top')


    # ---------------------------------------------------------
    # Panel C: PIP Distribution (Miami Plot with SymLog Scale)
    # ---------------------------------------------------------
    pip_df = pd.read_parquet(pip_vals_path)
    
    # Extract the PIP values, dropping any NaNs just to be safe
    pip_all = pip_df[pip_df['source'] == 'all']['PIP'].dropna().values
    pip_pleio = pip_df[pip_df['source'] == 'pleio']['PIP'].dropna().values

    median_all = np.median(pip_all)
    median_pleio = np.median(pip_pleio)

    perc_95_all = np.percentile(pip_all, 95)
    perc_95_pleio = np.percentile(pip_pleio, 95)
    
    # Calculate percentages. 
    # Multiply the pleio weights by -1 to force the histogram to draw downward.
    weight_all = np.ones_like(pip_all) / len(pip_all) * 100
    weight_pleio = (np.ones_like(pip_pleio) / len(pip_pleio) * 100) * -1
    
    # Plot 'Up' (All Causal Variants) - Filled then Outlined
    ax_c.hist(pip_all, bins=50, weights=weight_all, histtype='bar', 
              alpha=0.6, color=colors[2], label='All Causal Variants', edgecolor='black', linewidth=1.0)

    # Plot 'Down' (Pleiotropic Variants) - Filled then Outlined
    ax_c.hist(pip_pleio, bins=50, weights=weight_pleio, histtype='bar', 
              alpha=0.6, color=colors[3], label='Pleiotropic Variants', edgecolor='black', linewidth=1.0)
    
    # Draw a solid black line at y=0 to act as the mirror axis
    ax_c.axhline(0, color='black', linewidth=1)

    # Median lines
    ax_c.axvline(median_all, ymin=0.5, ymax=1.0, color='black', linestyle='--', linewidth=1.5)
    ax_c.axvline(median_pleio, ymin=0.0, ymax=0.5, color='black', linestyle='--', linewidth=1.5)
    ax_c.text(median_all + 0.03, 20, f'Median: {median_all:.3g}', color='black', va='center', fontsize = 8)
    ax_c.text(median_pleio + 0.03, -20, f'Median: {median_pleio:.3g}', color='black', va='center', fontsize = 8)

    # 95 Percentile lines
    ax_c.axvline(perc_95_all, ymin=0.5, ymax=0.76, color='black', linestyle='--', linewidth=1.5)
    ax_c.axvline(perc_95_pleio, ymin=0.5, ymax=0.2, color='black', linestyle='--', linewidth=1.5)
    ax_c.text(perc_95_all + 0.03, 0.76, f'95th percentile: {perc_95_all:.3g}', color='black', va='center', fontsize = 8)
    ax_c.text(perc_95_pleio + 0.03, -2.5, f'95th percentile: {perc_95_pleio:.3g}', color='black', va='center', fontsize = 8)
    
    # Apply Symmetric Log Scale to Y-axis to visualize the skewed tail
    ax_c.set_yscale('symlog', linthresh=0.1) 
    
    # Force specific symmetrical log ticks so the axis is readable
    ax_c.set_yticks([-100, -10, -1, -0.1, 0, 0.1, 1, 10, 100])
    
    # Custom formatter to handle absolute values AND clean up log formatting
    def symlog_formatter(x, pos):
        if x == 0: 
            return "0"
        return f"{abs(x):g}" # :g strips trailing zeros naturally
        
    ax_c.yaxis.set_major_formatter(ticker.FuncFormatter(symlog_formatter))
    
    # Set symmetrical Y-limits slightly past 100% to ensure top/bottom balance
    ax_c.set_ylim(-150, 150)

    # Lock X-axis to exact probability bounds (with a small buffer so bars aren't cut off)
    ax_c.set_xlim(-0.02, 1.02)
    ax_c.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    # Labels and formatted Panel Letter
    ax_c.set_xlabel('Posterior Inclusion Probability (PIP)')
    ax_c.set_ylabel('Percentage of Variants (%)')
    ax_c.text(-0.15, 1.02, 'C.', transform=ax_c.transAxes, fontsize=16, fontweight='bold', va='top')
    
    # Publication-ready legend with solid white background and sharp corners
    main_legend = ax_c.legend(frameon=True, edgecolor='black', fancybox=False, framealpha=1.0, loc='upper right', fontsize=10)
    ax_c.add_artist(main_legend)
    
    star_line, = ax_c.plot([], [], marker='', color='black', linestyle='None', markersize=10)
    ax_c.legend([star_line], [fr'$p_{{95th\ \%ile\ diff}} = {perm_pip_pval:.3g}$'], frameon=False, 
                loc='lower right', bbox_to_anchor=(1.0, 0.72), fontsize=9)

    # Save as highly scalable vector graphic
    plt.savefig(f"{output_dir}/gtex_enrich_all.pdf", bbox_inches='tight', transparent=True)
    plt.close()


if __name__ == "__main__":

    output_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/plots"
    fisher_results_path = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/res_all_pleio/fisher_results.tsv"
    pip_vals_path = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/res_all_pleio/merged_pip_values.parquet"

    # perm_pip_pval = run_permutation_test(pip_vals_path)
    # print(perm_pip_pval)
    
    perm_pip_pval = 0.2720
    gtex_plot(fisher_results_path, pip_vals_path, perm_pip_pval, output_dir)