import os
import textwrap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.transforms as mtransforms
from matplotlib.colors import LinearSegmentedColormap 


def plot_chromosome_density(merged_df, out_dir, bin_size=1_000_000):
    """
    Generates a high-resolution chromosome density map of similarity scores.
    Aggregates data into specified genomic bins (default 1 Mb).
    """
    print(f"Binning data into {bin_size / 1_000_000:.0f} Mb windows...")
    
    # 1. Bin the positions
    merged_df["dom_sig_total"] = np.where(merged_df["dom_sig_total"] == 1, 0, merged_df["dom_sig_total"])
    merged_df['bin_start'] = (merged_df['pos'] // bin_size) * bin_size
    
    # 2. Aggregate the total_similarity scores (using mean to find hotspots)
    density_df = (
    merged_df.groupby(['chr', 'bin_start'])['dom_sig_total']
    .agg(high_count = lambda x: (x > 1).sum(),
         total_count = 'count').reset_index())
    
    density_df['proportion'] = density_df['high_count'] / density_df['total_count']
    density_df['proportion'] = density_df['proportion'].fillna(0)

    # Calculate the maximum position per chromosome to draw the track backgrounds
    chr_lengths = merged_df.groupby('chr')['pos'].max().to_dict()
    
    # 3. Global Plot Styling (Publication Standards)
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica'],
        'font.size': 7,
        'axes.linewidth': 0.5,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'pdf.fonttype': 42,  # Ensures text is editable in Illustrator
        'ps.fonttype': 42
    })
    
    # Create figure: 183mm is the standard double-column width (~7.2 inches)
    fig, ax = plt.subplots(figsize=(7.2, 5), dpi=300)
    
    # Set up the colormap and normalization based on score percentiles to prevent outliers from washing out the colors
    cmap = LinearSegmentedColormap.from_list('custom_gradient', [(0.443,	0.776,	0.545),
                                                                 (0.988,	0.824,	0.024),
                                                                 (0.945,	0.298,	0.133)
                                                                 ])
    vmin = density_df['proportion'].quantile(0.01)
    vmax = density_df['proportion'].quantile(0.99)
    norm = Normalize(vmin=vmin, vmax=vmax)
    
    chromosomes = range(1, 23)
    track_height = 0.6
    
    print("Drawing chromosome tracks...")
    
    # 4. Draw each chromosome
    for chrm in chromosomes:
        chr_data = density_df[density_df['chr'] == chrm]
        max_len = chr_lengths.get(chrm, 0)
        
        # Draw a light grey background track for the entire chromosome length
        # Adjust the 'facecolor' hex code to change the color of the gap regions (e.g., '#FFFFFF' for white, '#D3D3D3' for grey)
        bg_rect = patches.Rectangle((0, chrm - track_height/2), max_len, track_height, 
                                    facecolor="#9D9D9D", edgecolor='none', zorder=1)
        ax.add_patch(bg_rect)
        
        # Draw the colored bins
        for _, row in chr_data.iterrows():
            color = cmap(norm(row['proportion']))
            rect = patches.Rectangle((row['bin_start'], chrm - track_height/2), 
                                     bin_size, track_height, 
                                     facecolor=color, edgecolor='none', zorder=2)
            ax.add_patch(rect)
            
    # 5. Axis Formatting
    ax.set_yticks(list(chromosomes))
    ax.set_yticklabels([f"Chr {c}" for c in chromosomes])
    ax.set_ylim(22.5, 0.5) # Invert y-axis so Chr 1 is at the top
    
    # Format x-axis to show Megabases (Mb)
    ax.set_xlim(0, 250_000_000)

    # Added a descriptive main title with appropriate padding so it doesn't crowd the top chromosome
    ax.set_title('Density of d-Pleiotropic Loci Across Chromosomes', fontsize=12, fontweight='bold', pad=10)
    
    # Create x-ticks every 50 Mb
    xticks = np.arange(0, 250_000_001, 50_000_000)
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{int(x / 1_000_000)}" for x in xticks])
    ax.set_xlabel('Genomic Position (Mb)', fontweight='bold')
    
    # Remove top and right spines for a clean look
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # 6. Add Colorbar
    cbar_ax = fig.add_axes([0.92, 0.3, 0.02, 0.4])
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_ticks([norm.vmin, norm.vmax])
    cbar.set_ticklabels(['Low', 'High'], fontsize=6)
    cbar.ax.tick_params(size=2.5)
    cbar.set_label('Variant Density within 1Mb Window', fontsize=7, fontweight='bold', labelpad=8, rotation=270)
    cbar.outline.set_linewidth(0.5)
    cbar.ax.tick_params(width=0.5, labelsize=6)
    
    gap_box = patches.Rectangle((0, -0.15), 1, 0.08, transform=cbar_ax.transAxes, 
                                facecolor='#808080', edgecolor='black', linewidth=0.5, clip_on=False)
    cbar_ax.add_patch(gap_box)
    
    # Tell matplotlib to move the main colorbar ticks and their values to the left side
    cbar.ax.yaxis.set_ticks_position('left')
    cbar.ax.yaxis.set_label_position('right')
    
    # Move the custom GAP tick to the left edge (0.0) and extend it slightly outward to the left (-0.2)
    cbar_ax.plot([0.0, -0.2], [-0.11, -0.11], color='black', linewidth=0.5, transform=cbar_ax.transAxes, clip_on=False)

    # Move the 'GAP' text to the left side, right-aligning it so it sits perfectly next to the new tick
    cbar_ax.text(-0.3, -0.11, 'Gap', transform=cbar_ax.transAxes, fontsize=6, va='center', ha='right')

    # 7. Save High-Resolution Plot
    output_file = os.path.join(out_dir, "density_plot.png")

    plt.subplots_adjust(left=0.1, right=0.88, top=0.95, bottom=0.1)
    plt.savefig(output_file, bbox_inches='tight', transparent=False)
    plt.close()
    
    print(f"Density map successfully saved to: {output_file}")


def plot_manhattan(merged_df, out_dir):
    # 1. Prepare Dataset (Cumulative Positions)
    merged_df = merged_df.sort_values(['chr', 'pos'])
    
    # Calculate chromosome lengths and cumulative offsets
    chr_lengths = merged_df.groupby('chr')['pos'].max()
    offsets = chr_lengths.cumsum().shift(1).fillna(0)
    merged_df['BPcum'] = merged_df['pos'] + merged_df['chr'].map(offsets)
    
    # Calculate centers for X-axis labels
    axis_df = merged_df.groupby('chr')['BPcum'].agg(['min', 'max']).mean(axis=1)

    # 2. Publication Style Settings
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica'],
        'font.size': 7,
        'axes.linewidth': 0.5,
        'pdf.fonttype': 42
    })

    fig, ax = plt.subplots(figsize=(8, 4), dpi=600)

    # 3. Plotting
    colors = ['#808080', '#87CEEB'] # Grey and Skyblue
    for i, (chrm, group) in enumerate(merged_df.groupby('chr')):
        ax.scatter(group['BPcum'], group['dom_sig_total'], 
                   c=colors[i % 2], s=12, alpha=0.8, edgecolors='none')
    
    ax.text(1, 0.92, r'$p < 4.72 \times 10^{-11}$', 
            transform=ax.transAxes, fontsize=7, ha='right')

    # Titles
    ax.set_title('d-Pleiotropic Variants', fontsize=10, fontweight='bold', pad=10)
    ax.set_xlabel('Chromosome', fontsize=7, fontweight='bold')
    ax.set_ylabel('Number of Significant Traits', fontsize=7, fontweight='bold')
    
    # Y-axis styling
    y_ticks = [0, 5, 10, 15, 20]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([str(y) for y in y_ticks])
    ax.set_ylim(0 - 0.5, merged_df['dom_sig_total'].max() + 5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', color='gray', linestyle='--', linewidth=0.3, alpha=0.5)
    ax.tick_params(axis='y', labelsize=6)

    # X-axis styling
    max_pos = merged_df['BPcum'].max()
    buffer = max_pos * 0.02 
    ax.set_xlim(-buffer, max_pos + buffer)   
    ax.set_xticks(axis_df.values)
    ax.set_xticklabels(axis_df.index)
    ax.tick_params(axis='x', labelsize=6)
    ax.set_xticks(axis_df.values)
    ax.set_xticklabels(axis_df.index)
    
    
    output_file = os.path.join(out_dir, "man_plot.png")

    plt.subplots_adjust(left=0.1, right=0.88, top=0.95, bottom=0.1)
    plt.savefig(output_file, dpi=600, bbox_inches='tight', transparent=False)
    plt.close()
    print(f"Manhattan plot saved to: {output_file}")


def plot_desc_percentages(desc_file_path, out_dir):
    df = pd.read_excel(desc_file_path)

    total_QCed = df.loc[df['chr'] == 'Total', 'QCed variants'].values[0]

    plot_df = df[(df['chr'] != 'Total') & (df['chr'] != 'X')].copy()
    
    # Calculate percentages
    all_cols = ['add_sig', 'add_pleiotropy', 'dom_sig', 'dom_pleiotropy', 'dom_pleiotropy_category']
    for col in all_cols:
        plot_df[f'{col}_pct'] = (plot_df[col] / plot_df['QCed variants']) * 100
        
    # Setup subplots: 2 rows, 1 column, shared X-axis
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica'],
        'font.size': 10, 'axes.linewidth': 0.5, 'pdf.fonttype': 42
    })
        
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(14, 8), sharex=True, gridspec_kw={'height_ratios': [1, 1]})
    fig.subplots_adjust(hspace=0.25) # Close the gap between plots
    
    x = np.arange(len(plot_df['chr']))
    width = 0.3 # Wider bars since we are splitting them up
    
    # --- PANEL A: Additive (Top) ---
    ax1.bar(x - width/2, plot_df['add_sig_pct'], width, label='Add. Sig.', color="#3C5488", edgecolor='black', linewidth=0.5)
    ax1.bar(x + width/2, plot_df['add_pleiotropy_pct'], width, label='Add. Pleiotropy', color="#DC0000", edgecolor='black', linewidth=0.5)
    
    ax1.set_ylabel('Variants (%)')
    ax1.spines[['top', 'right']].set_visible(False)
    ax1.legend(frameon=False, bbox_to_anchor=(1.01, 1), loc='upper left')
    
    # Optional: Annotate N= total on the top plot
    text_y = ax1.get_ylim()[1] * 1.05
    ax1.vlines(x, ymin=0, ymax=text_y, colors='gray', linestyles='dashed', linewidth=0.5, alpha=0.5, zorder=0)
    for pos, n in zip(x, plot_df['QCed variants']):
        ax1.text(pos, text_y, f"N={n:,}", ha='center', va='bottom', fontsize=8, rotation=45)

    # --- PANEL B: Dominant (Bottom) ---
    width_dom = 0.2
    dom_max = plot_df[['dom_sig_pct', 'dom_pleiotropy_pct', 'dom_pleiotropy_category_pct']].max().max() + 0.02
    ax2.bar(x - width_dom, plot_df['dom_sig_pct'], width_dom, label='Dom. Sig.', color="#00A087", edgecolor='black', linewidth=0.5)
    ax2.bar(x, plot_df['dom_pleiotropy_pct'], width_dom, label='Dom. Pleiotropy', color="#B09C85", edgecolor='black', linewidth=0.5)
    ax2.bar(x + width_dom, plot_df['dom_pleiotropy_category_pct'], width_dom, label='Dom. Pleio. Category', color="#8491B4", edgecolor='black', linewidth=0.5)
    ax2.vlines(x, ymin=0, ymax=dom_max, colors='gray', linestyles='dashed', linewidth=0.5, alpha=0.5, zorder=0)
    
    # Titles
    ax1.set_title("A. Additive variants", loc='left', fontweight='bold', fontsize=12, pad=40)
    ax2.set_title("B. Dominance variants", loc='left', fontweight='bold', fontsize=12, pad= 15)

    ax2.set_ylabel('Variants (%)')
    ax2.set_xlabel('Chromosome')
    ax2.set_xticks(x)
    ax2.set_xticklabels(plot_df['chr'])
    ax2.spines[['top', 'right']].set_visible(False)
    ax2.legend(frameon=False, bbox_to_anchor=(1.01, 1), loc='upper left')
    
    plot_path = f"{out_dir}/desc_plot.png"
    plt.savefig(plot_path, dpi=600, bbox_inches='tight')
    plt.close()


def plot_pleiotropy_matrix(merged_df, phen_names, out_dir, chromosomes=list(range(1, 23)), bin_size=1_000_000):
    print("1. Loading trait dictionary...")
    names_list = pd.read_excel(phen_names, usecols=["phenotype_code", "description", "category"])
    trait_dict = dict(zip(names_list['phenotype_code'].astype(str), names_list['description']))
    category_dict = dict(zip(names_list['description'], names_list['category']))
    
    print(f"2. Preparing data for {len(chromosomes)} chromosomes...")
    df = merged_df[merged_df['chr'].isin(chromosomes)].copy()
    
    # Standardize density calculation
    df["dom_sig_total"] = np.where(df["dom_sig_total"] == 1, 0, df["dom_sig_total"])
    df['bin_start'] = (df['pos'] // bin_size) * bin_size
    
    # Calculate continuous X-axis (Cumulative Position)
    chr_lengths = df.groupby('chr')['pos'].max()
    offsets = chr_lengths.cumsum().shift(1).fillna(0)
    df['BPcum'] = df['pos'] + df['chr'].map(offsets)
    
    # Calculate centers for X-axis labels
    axis_df = df.groupby('chr')['BPcum'].agg(['min', 'max']).mean(axis=1)
    
    # Calculate density per bin
    density_df = (df.groupby(['chr', 'bin_start'])
                  .agg(
                      high_count = ('dom_sig_total', lambda x: (x > 1).sum()), 
                      total_count = ('dom_sig_total', 'count'),
                      BPcum = ('BPcum', 'min')
                    ).reset_index())
    
    density_df['proportion'] = density_df['high_count'] / density_df['total_count']
    density_df['proportion'] = density_df['proportion'].fillna(0)
    
    # Isolate strictly the Hubs (proportion > 0)
    hub_bins = density_df[density_df['proportion'] > 0].copy()
    
    print("3. Extracting and mapping traits for all hubs...")
    records = []
    for _, row in hub_bins.iterrows():
        chrm = row['chr']
        b_start = row['bin_start']
        b_cum = row['BPcum']
        prop = row['proportion']
        
        # Get all variants in this specific 1Mb hub
        vars_in_bin = df[(df['chr'] == chrm) & (df['bin_start'] == b_start)]
        
        # Count how many variants map to each trait in this hub
        trait_counts = {}
        for t_str in vars_in_bin['sig_dom_traits'].dropna():
            codes = [c.strip() for c in str(t_str).split(',')]
            for c in codes:
                mapped_trait = trait_dict.get(c, c) 
                trait_counts[mapped_trait] = trait_counts.get(mapped_trait, 0) + 1
                
        # Store a record for every trait connected to this hub, including the SNP count
        for t, count in trait_counts.items():
            cat = category_dict.get(t, "Other")
            records.append({'chr': chrm, 'BPcum': b_cum, 'proportion': prop, 'trait': t, 'n_snps': count, 'category': cat})

    matrix_df = pd.DataFrame(records)

    # Calculate bin-specific normalization for n_snps (0 to 1 scale within each hub)
    # If a hub only has traits with the exact same SNP count, it defaults to 1.0 (red)
    matrix_df['n_snps_norm'] = matrix_df.groupby(['chr', 'BPcum'])['n_snps'].transform(
        lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() > x.min() else 1.0
    )
    
    # Sort traits alphabetically for the Y-axis (or by frequency)
    unique_traits_df = matrix_df[['category', 'trait']].drop_duplicates().sort_values(by=['category', 'trait'], ascending=[False, False])
    unique_traits = unique_traits_df['trait'].tolist()
    
    trait_to_y = {t: i for i, t in enumerate(unique_traits)}
    matrix_df['y_pos'] = matrix_df['trait'].map(trait_to_y)


    print("4. Rendering  matrix...")
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica'],
        'font.size': 7, 'axes.linewidth': 0.5, 'pdf.fonttype': 42
    })
    
    # Create dual-panel figure (Top = Density, Bottom = Matrix)
    # Height ratios: Matrix panel is 5 times taller to fit all traits
    fig, (ax_top, ax_bottom) = plt.subplots(nrows=2, ncols=1, figsize=(7.2, 9), dpi=600, 
                                            sharex=True, gridspec_kw={'height_ratios': [1, 20]})

    fig.subplots_adjust(hspace=0.05) # Bring panels very close together
    
    # Global Density Track
    cmap = LinearSegmentedColormap.from_list('custom_gradient', [(0.443,	0.776,	0.545),
                                                                 (0.988,	0.824,	0.024),
                                                                 (0.945,	0.298,	0.133)
                                                                 ])
    
    vmin = density_df['proportion'].quantile(0.01)
    vmax = density_df['proportion'].quantile(0.99)
    norm = Normalize(vmin=vmin, vmax=vmax)
    
    track_height = 1

    # Draw the continuous base track in the lighter color for the entire genome
    ax_top.add_patch(patches.Rectangle((0, 0), df['BPcum'].max(), track_height, 
                                       facecolor=(0.722, 0.888, 0.773), edgecolor='none'))
    

    # Overlay the darker intensity color specifically for odd chromosomes (1, 3, 5...)
    base_color = (0.443, 0.776, 0.545)
    for chrm in chromosomes:
        if chrm % 2 != 0:
            c_min = df[df['chr'] == chrm]['BPcum'].min()
            c_max = df[df['chr'] == chrm]['BPcum'].max()
            ax_top.add_patch(patches.Rectangle((c_min, 0), c_max - c_min, track_height, 
                                               facecolor=base_color, edgecolor='none'))
        
    # Draw colored density bins
    for _, row in hub_bins.iterrows():
        color = cmap(norm(row['proportion']))
        visual_width = bin_size * 1.5 
        
        rect = patches.Rectangle((row['BPcum'], 0), 
                                 visual_width, 
                                 track_height, 
                                 facecolor=color, 
                                 edgecolor=color,
                                 linewidth=0.5) 
        ax_top.add_patch(rect)
        
    ax_top.set_ylim(0, track_height)
    ax_top.set_yticks([])
    ax_bottom.tick_params(axis='y', length=5, width=0.8, pad=4)
    ax_top.set_title('Variant Density in 1Mb Window', fontweight='bold', fontsize=7, loc='center', pad=6)
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    ax_top.spines['left'].set_visible(False)
    
    # --- BOTTOM PANEL: Trait Matrix ---
    # Add subtle alternating background shading to group chromosomes visually
    for i, chrm in enumerate(chromosomes):
        if i % 2 == 0:
            c_min = df[df['chr'] == chrm]['BPcum'].min()
            c_max = df[df['chr'] == chrm]['BPcum'].max()
            ax_bottom.axvspan(c_min, c_max, facecolor="#F4F4F4", edgecolor='none', zorder=0)

    # Plot the traits as a dense scatter matrix
    # Plot using the bin-normalized values (0 to 1) with 'coolwarm' (blue to red)
    ax_bottom.scatter(matrix_df['BPcum'], matrix_df['y_pos'], 
                      marker='|', c=matrix_df['n_snps_norm'], cmap='coolwarm', 
                      vmin=0, vmax=1, s=20, lw=1, zorder=2)

    # Format Y-axis with Trait Names
    ax_bottom.set_yticks(range(len(unique_traits)))
    ax_bottom.set_yticklabels(unique_traits, fontsize=5) # 5pt font to fit them all cleanly
    ax_bottom.set_ylim(-1, len(unique_traits))

    #create a color palette for categories
    unique_categories = unique_traits_df['category'].unique()
    custom_palette = ["#3C5488", '#DC0000','#00A087',"#89603D",'#8491B4', '#91D1C2',
                       '#631879',"#B09C85",'#00A05B',"#E64B35","#C59316", '#4DBBD5']
    
    # 3. Apply the 70% opacity to all colors at once
    new_palette = [mcolors.to_rgba(c, alpha=1) for c in custom_palette]
    
    # 4. Reverse the categories so the first color applies to the visual top of the plot
    top_to_bottom_cats = list(reversed(unique_categories))
    
    # 5. Map them together safely (the % len(new_palette) ensures it loops if you have >11 categories)
    category_color_dict = {cat: new_palette[i % len(new_palette)] for i, cat in enumerate(top_to_bottom_cats)}

    trans = mtransforms.blended_transform_factory(ax_bottom.transAxes, ax_bottom.transData)

    # Loop through the exact Y-positions and traits
    for y_pos, trait in enumerate(unique_traits):
        
        # Look up the color for this trait's category
        cat = unique_traits_df[unique_traits_df['trait'] == trait]['category'].values[0]
        color = category_color_dict[cat]
        
        # Draw a colored block exactly where the tick mark sits
        # X starts at -0.03 (just outside the plot) and width=0.03 brings it flush to the axis line
        # Y starts at y_pos - 0.5 with height=1.0 so the blocks touch each other perfectly without gaps
        rect = patches.Rectangle((-0.008, y_pos - 0.5), width=0.008, height=1.0, 
                                 transform=trans, facecolor=color, edgecolor='none', 
                                 clip_on=False, alpha=0.8)
        ax_bottom.add_patch(rect)

    # Format X-axis with Chromosome Centers
    ax_bottom.set_xlim(0, df['BPcum'].max())
    ax_bottom.set_xticks(axis_df.values)
    ax_bottom.set_xticklabels(axis_df.index)
    ax_bottom.set_xlabel('Chromosome', fontweight='bold', labelpad=8)

    # Clean up matrix borders
    ax_bottom.spines['top'].set_visible(False)
    ax_bottom.spines['right'].set_visible(False)
    ax_bottom.grid(axis='y', color='gray', linestyle=':', linewidth=0.3, alpha=0.5)

    # y-axis title
    ax_bottom.text(-0.02, 1.02, 'Phenotypes', 
                   transform=ax_bottom.transAxes, 
                   fontsize=7, fontweight='bold', 
                   ha='right', va='bottom')
    
    # number of hits legend
    cbar_ax_bottom = fig.add_axes([1.03, 0.2, 0.02, 0.2]) 
    
    # Map the coolwarm colormap (from 0 to 1) to match your bin-specific normalization
    sm_bottom = plt.cm.ScalarMappable(cmap='coolwarm', norm=Normalize(vmin=0, vmax=1))
    sm_bottom.set_array([])
    
    # Hit Colorbar Legend
    cbar_bottom = fig.colorbar(sm_bottom, cax=cbar_ax_bottom)
    cbar_bottom.set_ticks([0, 1])
    cbar_bottom.set_ticklabels(['Low','High'], fontsize=6)
    cbar_bottom.ax.set_title('Hit Density\nin 1Mb Window', fontsize=7, fontweight='bold', pad=10, loc='center')
    cbar_bottom.outline.set_linewidth(0.5)
    cbar_bottom.ax.tick_params(width=0.5, size=2.5)

    # Variant Density legend
    cbar_ax_top = fig.add_axes([1.0, 0.8442, 0.08, 0.01]) 
    sm_top = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm_top.set_array([])
    cbar_top = fig.colorbar(sm_top, cax=cbar_ax_top, orientation='horizontal')
    cbar_top.ax.xaxis.set_ticks_position('top')
    cbar_top.set_ticks([norm.vmin, norm.vmax])
    cbar_top.set_ticklabels(['Low', 'High'], fontsize=6)
    cbar_top.outline.set_linewidth(0.5)
    cbar_top.ax.tick_params(width=0.5, size=2.5)


    # Categories Legend
    visual_order_categories = reversed(unique_traits_df['category'].unique())
    legend_elements = [patches.Patch(facecolor=category_color_dict[cat], edgecolor='black', 
                                     label=textwrap.fill(str(cat), width=30, break_long_words=False,
                                                         break_on_hyphens=False), 
                                     linewidth=0.5) 
                       for cat in visual_order_categories]
    
    
    ax_bottom.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.03, 0.93), 
                     ncol=1, fontsize=7, frameon=False, labelspacing=0.8,
                     title="Phenotype Categories", title_fontproperties={'weight': 'bold', 'size': 7})

    # Output
    output_file = os.path.join(out_dir, "pleiotropy_matrix.png")
    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    plt.close()
    print(f"Matrix successfully saved to: {output_file}")


if __name__ == "__main__":
   
    sig_SNPs_path = "/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/significant_SNPs/all_sig_SNPs.tsv.gz"
    output_dir = "/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/results"
    phen_names = "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/phen_dict_renamed.xlsx"
    desc_file_path = "/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/results/SNP_descriptives_plotting.xlsx"


    sig_SNPs_df = pd.read_csv(sig_SNPs_path, sep="\t", compression="gzip", 
                            usecols=["variant", "chr", "pos","rsid", 
                                     "dom_sig_total", "sig_dom_traits"],
                                     dtype={"sig_dom_traits": str})

    plot_manhattan(sig_SNPs_df, output_dir)
    plot_chromosome_density(sig_SNPs_df, output_dir)
    plot_pleiotropy_matrix(sig_SNPs_df, phen_names, output_dir)
    plot_desc_percentages(desc_file_path, output_dir)