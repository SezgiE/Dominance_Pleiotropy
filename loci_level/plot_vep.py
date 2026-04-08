import os
import textwrap
import pandas as pd
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


# ==========================================
if __name__ == "__main__":

    upset_plot(
        coloc_filepath="/Users/sezgi/Documents/dominance_pleiotropy/loci_level/coloc_results/snp_info.tsv",
        vep_filepath="/Users/sezgi/Documents/dominance_pleiotropy/loci_level/coloc_results/vep_res.txt",
        out_dir="/Users/sezgi/Documents/dominance_pleiotropy/loci_level/coloc_results/vep_plots"
    )