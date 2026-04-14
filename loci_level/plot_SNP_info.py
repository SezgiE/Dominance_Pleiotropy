import os
import textwrap
import pandas as pd
import numpy as np
import upsetplot as us
import matplotlib.pyplot as plt


def set_style():
    """Hardcodes matplotlib parameters to format the plot."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 12,
        'axes.labelsize': 10,
        'xtick.labelsize': 7,
        'ytick.labelsize': 8,
        'axes.linewidth': 1.0,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'pdf.fonttype': 42,
        'ps.fonttype': 42
    })


def set_style2():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['font.size'] = 12
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams['xtick.major.width'] = 1.0
    plt.rcParams['ytick.major.width'] = 1.0
    plt.rcParams['pdf.fonttype'] = 42


def set_style3():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['font.size'] = 10
    plt.rcParams['xtick.labelsize'] = 8
    plt.rcParams['ytick.labelsize'] = 8
    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams['xtick.major.width'] = 1.0
    plt.rcParams['ytick.major.width'] = 1.0
    plt.rcParams['pdf.fonttype'] = 42


def upset_plot(coloc_filepath, vep_filepath, out_dir, filename="pleiotropy_upset.png"):
    """
    Loads raw Coloc and VEP files, merges them, and generates an UpSet plot 
    visualizing pleiotropic variant sharing among traits.

    Parameters:
    - coloc_filepath (str): Path to the coloc summary TSV file.
    - vep_filepath (str): Path to the cleaned VEP results TXT/TSV file.
    - out_dir (str): Directory where the plot will be saved.
    - filename (str): Name of the output image file.
    """
    print("[*] Loading and merging datasets...")
    
    # --- Load and Clean Data ---
    coloc_df = pd.read_csv(coloc_filepath, sep='\t') 
    
    # Check how many '##' metadata lines exist so Pandas can skip them
    with open(vep_filepath, 'r') as f:
        skip_lines = sum(1 for line in f if line.startswith('##'))
        
    vep_df = pd.read_csv(vep_filepath, sep='\t', skiprows=skip_lines)
    vep_df.rename(columns={'#Uploaded_variation': 'Variant_ID'}, inplace=True)

    
    # Convert comma-separated phen_codes into actual Python lists for the UpSet logic
    coloc_df['categories'] = coloc_df['phen_categories'].apply(
            lambda x: [textwrap.fill(code.strip(), width=25) for code in str(x).split(',')] if pd.notnull(x) else []
        )
    
    # Merge on variant identifiers
    merged_df = pd.merge(
        coloc_df, 
        vep_df, 
        left_on='variant', 
        right_on='Variant_ID', 
        how='inner'
    )
    
    if merged_df.empty:
        raise ValueError("The merged dataset is empty. Check that the variant IDs match between your Coloc and VEP files.")

    # --- Prepare Data for UpSet Plot ---
    print("[*] Preparing data matrix for UpSet plot...")
    
    # Get unique variants, their associated traits and consequences
    unique_vars_df = merged_df[['Variant_ID']].drop_duplicates()
    var_traits_map = merged_df.groupby('Variant_ID')['categories'].first().to_dict()
    var_consq_map = merged_df.groupby('Variant_ID')['Consequence'].first().fillna('unknown').apply(lambda x: str(x).split(',')[0]).to_dict()
    
    # Extract all unique traits across the entire dataset to build the matrix columns
    all_traits = set()
    for traits_list in merged_df['categories']:
        all_traits.update(traits_list)
    unique_traits = sorted(list(all_traits))

    # Build the binary indicator matrix (Rows = Variants, Columns = Traits, Values = True/False)
    indicator_data = []
    for variant_id in unique_vars_df['Variant_ID']:
        traits_list = var_traits_map.get(variant_id, [])

        raw_consq = var_consq_map.get(variant_id, 'Unknown')
        clean_consq = str(raw_consq).replace('_', ' ').capitalize()
        clean_consq = clean_consq.replace('utr', 'UTR').replace('Utr', 'UTR')
        
        # Consequence into the row data
        indicator_row = {
            'Variant_ID': variant_id, 
            'Consequence': clean_consq
        }
        
        for trait in unique_traits:
            indicator_row[trait] = True if trait in traits_list else False
        indicator_data.append(indicator_row)

    indicator_df = pd.DataFrame(indicator_data)
    indicator_df.set_index(unique_traits, inplace=True)


    # --- 3. Generate and Save Plot ---
    print("[*] Generating UpSet plot...")
    os.makedirs(out_dir, exist_ok=True)

    # Define plot size
    set_style()
    fig_size = (15, 8) 
    fig = plt.figure(figsize=fig_size)
    

    # Initialize UpSet
    upset = us.UpSet(
        indicator_df, 
        subset_size='count', 
        intersection_plot_elements=0, 
        element_size=None,
        facecolor="black",  
        shading_color=0.1,      
        show_counts=False
    )

    upset.add_stacked_bars(by='Consequence', colors='tab20', title="Variant Count", elements=4)
    axes_dict = upset.plot(fig=fig)


    # Legend setup for the stacked bars
    for ax_name, ax in axes_dict.items():
        if ax_name not in ['matrix', 'shading', 'totals']:
            legend = ax.get_legend()
            if legend:
                handles, labels = ax.get_legend_handles_labels()

                sorted_pairs = sorted(zip(handles, labels), key=lambda x: x[1])
                handles, labels = zip(*sorted_pairs)

                ax.legend(
                    handles, 
                    labels, 
                    loc='upper left', 
                    bbox_to_anchor=(-0.4, 1.15), 
                    fontsize=8,
                    labelspacing=0.6,   
                    frameon=True,
                    title="Variant Consequences",
                )
                plt.setp(legend.get_texts(), fontsize='6')
            
            # Add counts on top of the stacked bars
            from collections import defaultdict
            totals = defaultdict(float)
            for p in ax.patches:
                totals[p.get_x() + p.get_width()/2] += p.get_height()
            
            for x, total in totals.items():
                if total > 0:
                    ax.text(x, total, f'{int(total)}', ha='center', va='bottom', 
                            fontsize=8, fontweight='bold')


    # Add titles and labels to the current active figure
    plt.suptitle('Pleiotropic Variants across Trait Categories', y=0.95, fontsize=14, fontweight='bold',
                 ha='center', family='sans-serif')
    
    # Save plot with tight bbox to prevent cutoff of trait names
    out_path = os.path.join(out_dir, filename)
    plt.savefig(out_path, bbox_inches='tight', dpi=600) # Added 300 dpi for publication quality
    plt.close(fig) 
    
    print(f"[*] Success! Pleiotropy UpSet plot saved to: {out_path}")


def snp_3D_plot(result_df_path, output_path):
    
    # Data Prep
    result_df = pd.read_csv(result_df_path, sep='\t')
    result_df["std_dom_Bsq"] = result_df["std_dom_b"] ** 2 
    result_df = result_df.sort_values(by=["category", "std_dom_Bsq"], ascending=[True, True])
    
    categories = result_df["category"].unique()
    cat_mapping = {cat: i for i, cat in enumerate(categories)}
    result_df["cat_num"] = result_df["category"].map(cat_mapping)
    
    set_style2()
    
    # Initialize Multi-panel Figure
    fig = plt.figure(figsize=(16, 8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.4], wspace=0.1, hspace=0.3)

    # --- PANEL A: MAF Distribution ---
    ax1 = fig.add_subplot(gs[0, 0])
    weights1 = np.ones_like(result_df["maf"]) / len(result_df["maf"])
    ax1.hist(result_df["maf"], bins=30, weights=weights1, color="#2c7fb8", edgecolor='black', alpha=0.7)
    
    maf_median = result_df["maf"].median()
    ax1.axvline(maf_median, color='red', linestyle='--', linewidth=1)
    ax1.text(maf_median + 0.02, ax1.get_ylim()[1]*0.8, f'Median: {maf_median:.3f}', color='red', fontsize=9)
    
    ax1.set_xlim(-0.02, 0.5)
    ax1.set_ylim(0, 0.3)
    ax1.set_yticks([0, 0.1, 0.2, 0.3])
    ax1.set_xlabel('MAF')
    ax1.set_ylabel('Proportion of Variants')
    ax1.text(-0.12, 1.1, 'A', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top')


    # --- PANEL B: Effect Size Distribution ---
    ax2 = fig.add_subplot(gs[1, 0])
    clipped_data = result_df["std_dom_Bsq"].clip(upper=0.01)
    weights2 = np.ones_like(clipped_data) / len(clipped_data)
    ax2.hist(clipped_data, bins=30, weights=weights2, color="#7fcdbb", edgecolor='black', alpha=0.7)
    
    beta_median = result_df["std_dom_Bsq"].median()
    ax2.axvline(beta_median, color='red', linestyle='--', linewidth=1)
    ax2.text(beta_median + 0.0002, ax2.get_ylim()[1]*0.8, f'Median: {beta_median:.4f}', color='red', fontsize=9)
    
    ax2.set_ylim(0, 0.4)
    ax2.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
    ax2.set_xlim(-0.00035, 0.01)
    ax2.set_xticks([0, 0.002, 0.004, 0.006, 0.008, 0.010001])
    ax2.set_xticklabels(['0', '0.002', '0.004', '0.006', '0.008', '>0.01'])
    ax2.set_xlabel('Squared Std. Effect Size (\u03B2\u00B2)')
    ax2.set_ylabel('Proportion of Variants')
    ax2.text(-0.12, 1.1, 'B', transform=ax2.transAxes, fontsize=16, fontweight='bold', va='top')


    # --- PANEL C: 3D Visualization ---
    ax = fig.add_subplot(gs[:, 1], projection='3d')
    ax.set_box_aspect([2.0, 2.0, 1.8])
    ax.dist = 7
    colors_cmap = plt.get_cmap('tab10')
    
    # Scatter plot
    scatter = ax.scatter(result_df["maf"], 
                         result_df["cat_num"], 
                         result_df["std_dom_Bsq"], 
                         c=result_df["cat_num"], 
                         cmap='tab10', 
                         s=55, 
                         alpha=0.7, 
                         edgecolors='w', 
                         linewidth=0.3)

    # Drop lines to floor
    norm = plt.Normalize(vmin=result_df["cat_num"].min(), vmax=result_df["cat_num"].max())
    cmap = plt.get_cmap('tab10')
    for x, y, z, c_idx in zip(result_df["maf"], result_df["cat_num"], result_df["std_dom_Bsq"], result_df["cat_num"]):
        # Pass the normalized index to the cmap to lock the colors together
        ax.plot([x, x], [y, y], [0, z], color=cmap(norm(c_idx)), alpha=0.6, linewidth=1.5)

    # 3D Labels & Formatting
    ax.set_xlabel('MAF', labelpad=10)
    ax.set_zlabel('', labelpad=0) 
    ax.text2D(-0.04, 0.5, 'Squared Std. Effect Size (\u03B2\u00B2)', 
              transform=ax.transAxes, rotation=90, verticalalignment='center', fontsize=12)
    
    ax.set_yticks(range(len(categories)))
    ax.set_yticklabels([])
    ax.xaxis.pane.fill = ax.yaxis.pane.fill = ax.zaxis.pane.fill = False
    ax.view_init(elev=20, azim=38)
    
    # Panel C Label
    ax.text2D(-0.05, 1.042, 'C', transform=ax.transAxes, fontsize=16, fontweight='bold', va='top')

    # Legend with wrapping
    legend = ax.legend(*scatter.legend_elements(), title="Trait Categories", 
                       loc="center left", bbox_to_anchor=(1, 0.5), frameon=False,
                       prop={'size': 10}, title_fontsize=12)
    
    for i, text in enumerate(legend.get_texts()):
        text.set_text("\n".join(textwrap.wrap(categories[i], width=25)))

    plt.savefig(output_path, dpi=600, bbox_inches='tight', format='pdf')


def plot_effect_direction(df_path, all_snps_df, output_filename):
    """
    Generates a 2-panel figure: 
    Panel A: Grouped bar chart showing count of significant Additive vs Dominant traits per SNP.
    Panel B: Ranked variance range plot showing the spread of dominance ratios.
    """

    pval_thresh = -np.log10((5e-8)/1060) 
    df = pd.read_csv(df_path, sep='\t')
    all_snps = pd.read_csv(all_snps_df, compression="gzip", sep='\t', 
                           usecols=['variant', "rsid", "add_sig_total"])

    # Dominance Ratio for valid additive traits
    df_valid = df[df['add_log10_pval'] > pval_thresh].copy()
    df_valid['dom_ratio'] = df_valid['std_dom_b'] / df_valid['std_add_b']
    df_valid['dom_ratio_clipped'] = df_valid['dom_ratio'].clip(lower=-2, upper=2)

    # Aggregate data per SNP for Panel B
    snp_stats = df_valid.groupby('variant')['dom_ratio_clipped'].agg(
        ['min', 'max', 'median']
    ).reset_index()

    # Calculate INDEPENDENT counts for Panel A
    sig_add_mask = df["add_log10_pval"] > pval_thresh

    add_counts = df[sig_add_mask].groupby('variant').size().reset_index(name='add_sig_count')
    dom_counts = df.groupby('variant').size().reset_index(name='dom_sig_count')

    # Merge counts into snp_stats
    snp_stats = pd.merge(snp_stats, add_counts, on='variant', how='left').fillna(0)
    snp_stats = pd.merge(snp_stats, dom_counts, on='variant', how='left')
    snp_stats = pd.merge(snp_stats, all_snps[['variant', 'add_sig_total']], on='variant', how='inner')
    snp_stats["sig_ratio"] = snp_stats["add_sig_count"] / (snp_stats["dom_sig_count"])

    # Sort SNPs strictly by their median dominance ratio (The S-curve)
    snp_stats = snp_stats.sort_values('median').reset_index(drop=True)

    # Initialize the multi-panel figure
    set_style()
    fig, (ax_A, ax_B, ax_C) = plt.subplots(
        nrows=3, ncols=1, figsize=(10, 7), sharex=True, 
        gridspec_kw={'height_ratios': [0.75, 2.5, 1.25]} # Middle panel stays the largest
    )
    colors_palette = ['#3C5488', '#E64B35', '#4DBBD5', '#00A087', '#7E818D']
    plt.subplots_adjust(hspace=0.3)
    x_pos = np.arange(len(snp_stats))


    # ==========================================
    # PANEL A: Grouped Bar Chart for Trait Counts
    # ==========================================
    ratio_colors = np.where(snp_stats['sig_ratio'] >= 1, '#00A087', '#FF9F1C')
    
    ax_A.bar(
        x_pos, 
        snp_stats['sig_ratio'], 
        width=1.0,
        color=ratio_colors,
        alpha=0.9,
        zorder=2,
        edgecolor='black', linewidth=0.1,
    )
    
    ax_A.set_ylabel(r'$N_{add}  /  N_{dom}$', labelpad=10)
    ax_A.spines['bottom'].set_visible(False)
    ax_A.tick_params(axis='x', length=0)
    
    # Add panel label
    ax_A.text(-0.1, 1.32, 'A', transform=ax_A.transAxes, fontsize=14, fontweight='bold', va='top')


    # ==========================================
    # PANEL B: Dominance Ratio Variance
    # ==========================================
    
    bar_colors = np.where((snp_stats['min'] < 0) & (snp_stats['max'] > 0), "#E64B35", '#3C5488')

    ax_B.vlines(
        x=x_pos, ymin=snp_stats['min'], ymax=snp_stats['max'], 
        color=bar_colors, linewidth=1, alpha=0.5, zorder=1
    )

    ax_B.scatter(
        x=x_pos, y=snp_stats['median'], 
        color='#3C5488', s=15, edgecolors='none', alpha=1.0, zorder=2
    )

    # Reference lines
    ax_B.axhline(0, color='black', linewidth=1.0, linestyle='-', zorder=0)
    ax_B.axhline(1, color='#A0A0A0', linewidth=0.8, linestyle='--', zorder=0)
    ax_B.axhline(-1, color='#A0A0A0', linewidth=0.8, linestyle='--', zorder=0)
    
    # Text annotations (dynamically placed at the end of the X-axis)
    label_x = len(snp_stats) - 5
    ax_B.text(label_x, 0.05, 'Additive', color='black', fontsize=9, va='bottom', ha='right')
    ax_B.text(label_x, 1.05, 'Dominant', color='#505050', fontsize=9, va='bottom', ha='right')
    ax_B.text(label_x, -0.95, 'Recessive', color='#505050', fontsize=9, va='bottom', ha='right')

    # Add label
    ax_B.text(-0.1, 1.08, 'B', transform=ax_B.transAxes, fontsize=14, fontweight='bold', va='top')

    # Format Axes
    ax_B.set_xticks([]) 
    ax_B.set_xlim(-2, len(snp_stats) + 2) 
    ax_B.set_ylabel(r'Dominance Ratio ($\beta_{dom} / \beta_{add}$)', labelpad=10)
    ax_B.set_ylim(-2.2, 2.2)


    # ==========================================
    # PANEL C: Total Additive Traits
    # ==========================================
    
    ax_C.bar(
        x_pos, snp_stats['add_sig_total'], width=0.7, color='#4DBBD5', 
        edgecolor='black', linewidth=0.2, alpha=1.0, zorder=2, label = "Additive"
    )

    ax_C.bar(
        x_pos+0.1, snp_stats['dom_sig_count'], 
        width=0.5, color='#E64B35', alpha=0.8, 
        edgecolor='black', linewidth=0.2, label='Dominant', zorder=3
    )
    
    ax_C.axhline(0, color='black', linewidth=1.0, zorder=1)
    ax_C.set_ylabel('Total Trait\nCount', labelpad=10)
    ax_C.legend(frameon=False, loc='upper left', fontsize=10)
    ax_C.set_ylim(0, 61)
    ax_C.set_yticks([0, 20, 40, 60])

    ax_C.grid(axis='y', color='#D0D0D0', linestyle='--', linewidth=0.8, zorder=0)
    ax_C.spines['bottom'].set_visible(False)
    ax_C.tick_params(axis='x', length=0)
    ax_C.text(-0.1, 1.2, 'C', transform=ax_C.transAxes, fontsize=14, fontweight='bold', va='top')

    # The master X-axis formatting is now anchored to Panel C
    ax_C.set_xlim(-2, len(snp_stats) + 2) 
    ax_C.set_xlabel(f'{snp_stats["variant"].nunique()} Pleiotropic SNPs (Ranked by Median Dominance Ratio)', labelpad=10)
    
    # Save
    plt.savefig(output_filename, format='pdf', bbox_inches='tight')


    # ==========================================
    # Descriptive Info
    # ==========================================
    is_not_in_stats = ~df["variant"].isin(snp_stats["variant"])
    unique_variants_df = df[is_not_in_stats][['variant', "phen_code"]]
    unique_variants_df = unique_variants_df.merge(all_snps[['variant', 'rsid']], on='variant', how='inner')
    print(unique_variants_df)

    # specific dominance effect patterns among the SNPs
    only_negative = len(snp_stats[snp_stats['max'] < 0])
    only_positive = len(snp_stats[snp_stats['min'] > 0])
    crosses_zero = len(snp_stats[(snp_stats['min'] < 0) & (snp_stats['max'] > 0)])
    print(f"1 - SNPs with only negative dominance ratio: {only_negative}")
    print(f"2 - SNPs with only positive dominance ratio: {only_positive}")
    print(f"3 - SNPs with both positive and negative ratio: {crosses_zero}")
    print("Crossed zero:")
    zero_snps = snp_stats[(snp_stats['min'] < 0) & (snp_stats['max'] > 0)]
    zero_snps = zero_snps.merge(all_snps[['variant', 'rsid']], on='variant', how='inner')
    for row in zero_snps.itertuples():
        print(f" - {row.variant} ({row.rsid})")

    # comparisons between additive and dominant trait counts
    equal_counts = len(snp_stats[snp_stats['add_sig_total'] == snp_stats['dom_sig_count']])
    add_higher = len(snp_stats[snp_stats['add_sig_total'] > snp_stats['dom_sig_count']])
    dom_higher = len(snp_stats[snp_stats['add_sig_total'] < snp_stats['dom_sig_count']])
    print(f"4 - SNPs where Additive count EQUALS Dominance count: {equal_counts}")
    print(f"5 - SNPs where Additive count is HIGHER than Dominance count: {add_higher}")
    print(f"6 - SNPs where Additive count is LOWER than Dominance count: {dom_higher}")


# ==========================================
if __name__ == "__main__":

    all_snps_df = "/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/significant_SNPs/all_sig_SNPs.tsv.gz"
    coloc_snps = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/coloc_results/coloc_snps.tsv"
    coloc_snps_info = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/coloc_results/coloc_snp_info.tsv"
    out_dir="/Users/sezgi/Documents/dominance_pleiotropy/loci_level/loci_results" 

    #upset_plot( coloc_snps, f"{out_dir}/vep_res.txt", out_dir)
    #snp_3D_plot(f"{out_dir}/snp_info.tsv", f"{out_dir}/snp_maf.pdf")
    plot_effect_direction(coloc_snps_info, all_snps_df, f"{out_dir}/snps_effect_direction.pdf")