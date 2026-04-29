import os
import glob
import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D



def set_style():
    """Hardcodes matplotlib parameters to format the plot."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 12,
        'axes.labelsize': 12,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'axes.linewidth': 1.0,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'pdf.fonttype': 42,
        'ps.fonttype': 42
    })


def std_expression(gtex_med_TPM_path, gtex_pleio_res_path, tissues_dir):
    
    
    gtex_med_tpm = pd.read_csv(gtex_med_TPM_path, sep='\t', skiprows=2)
    gtex_pleio = pd.read_csv(gtex_pleio_res_path, sep='\t')

    gene_list = gtex_pleio["gene_id"].unique()

    tissues = [
        os.path.basename(f).split(".v11")[0] 
        for f in glob.glob(os.path.join(tissues_dir, "*.parquet"))
        ]

    cols_to_select = ["Name", "Description"] + tissues
    
    gtex_tpm_pleio = gtex_med_tpm[gtex_med_tpm['Name'].isin(gene_list)].copy()


    # Standardize row-wise for numeric columns
    numeric_cols = gtex_tpm_pleio.select_dtypes(include=[np.number]).columns

    gtex_tpm_pleio_std = gtex_tpm_pleio.copy()
    gtex_tpm_pleio_std[numeric_cols] = gtex_tpm_pleio[numeric_cols].apply(
        lambda row: (row - row.mean()) / row.std(), 
        axis=1
        )
    
    return gtex_tpm_pleio_std


def hallmark_enrichment_analysis(df, output_dir):
    gene_symbols = df['Description'].dropna().unique().tolist()
    print(len(gene_symbols))

    enr = gp.enrichr(
        gene_list=gene_symbols,
        gene_sets=['MSigDB_Hallmark_2020','KEGG_2026', 
                   "Reactome_Pathways_2024", "GO_Biological_Process_2025",
                   "GO_Molecular_Function_2025", "GTEx_Tissues_V8_2023"],
        organism='human',
        outdir=None
    )
    
    # Extract results and filter for significance (FDR < 0.05)
    results = enr.results
    sig_results = results[results['Adjusted P-value'] < 0.05].copy()
    
    # Calculate -log10(FDR) for plotting
    sig_results['log10(FDR)'] = -np.log10(sig_results['Adjusted P-value'])
    
    # Sort by significance
    sig_results = sig_results.sort_values('log10(FDR)', ascending=True)
    sig_results[["GoIs", "Total"]] = sig_results["Overlap"].str.split("/", n=2, expand=True).iloc[:, :2].astype(int)
    sig_results["Overlap"] = sig_results["GoIs"] / sig_results["Total"]

    sig_results.to_csv(f'{output_dir}/gene_enrich_res.tsv', index=False, header = True, sep="\t")
    
    return sig_results


def plot_heatmap(df, output_dir):

    heatmap_data = df.set_index('Description').drop(columns=['Name']).T
    heatmap_data.index = heatmap_data.index.str.replace('_', ' ')
    heatmap_data = heatmap_data.sort_index(axis=0).sort_index(axis=1)
    data_matrix = heatmap_data.values
    
    set_style()
    fig, ax = plt.subplots(figsize=(len(heatmap_data.columns) * 0.25, len(heatmap_data) * 0.15 + 4))

    # Calculate symmetric bounds for the diverging colormap
    v_max = np.nanmax(np.abs(data_matrix))
    
    # Use pcolormesh instead of seaborn
    # Use pcolormesh with the Red-Yellow-Blue (reversed) colormap
    c = ax.pcolormesh(
        data_matrix, 
        cmap='RdYlBu_r', 
        vmin=-v_max, 
        vmax=v_max, 
        edgecolors='white', 
        linewidth=0.5
    )

    for i in range(5, data_matrix.shape[0], 5):
        ax.axhline(i, color='black', linewidth=0.5, linestyle='--', alpha=0.2, zorder=2)
        
    for j in range(5, data_matrix.shape[1], 5):
        ax.axvline(j, color='black', linewidth=0.5, linestyle='--', alpha=0.2, zorder=2)

    # Invert y-axis to match dataframe order (top-to-bottom)
    ax.invert_yaxis()

    # Center the ticks
    ax.set_xticks(np.arange(data_matrix.shape[1]) + 0.5)
    ax.set_yticks(np.arange(data_matrix.shape[0]) + 0.5)

    # Move x-axis to the top
    ax.xaxis.tick_top()

    # Apply labels (changed ha='right' to ha='left' so text doesn't overlap the plot when on top)
    ax.set_xticklabels(heatmap_data.columns, rotation=45, ha='left')
    ax.set_yticklabels(heatmap_data.index)

    ax.tick_params(axis='x', top=True, bottom=False, length=5, width=0.6)
    ax.tick_params(axis='y', left=True, right=False, length=3, width=0.5)

    for spine in ax.spines.values():
        spine.set_visible(False)

    # Add colorbar at the bottom
    cbar = fig.colorbar(c, ax=ax, orientation='horizontal', shrink=0.4, aspect=30, pad=0.02)
    cbar.outline.set_visible(False)
    cbar.set_label('Median Expression (Z-score)')
    cbar.outline.set_linewidth(0.2)

    # Clean axes
    ax.set_ylabel('')
    ax.set_xlabel('')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/expression_heatmap.pdf', dpi=600, bbox_inches='tight')


def plot_publication_mirror(df, output_dir):
    
    # Data Prep
    df['Term'] = df['Term'].str.replace(r'\s*\([^)]*\)', '', regex=True)
    df['Term'] = df['Term'].str.title()
    df_sorted = df.sort_values(['Gene_set', 'Term'], ascending=[False, False]).reset_index(drop=True)
    
    colors = ['#3C5488', '#E64B35', '#4DBBD5', '#00A087', '#7E818D']
    unique_dbs = df_sorted['Gene_set'].unique()
    db_color_map = {db: colors[i % len(colors)] for i, db in enumerate(unique_dbs)}
    row_colors = df_sorted['Gene_set'].map(db_color_map)

    all_genes = set()

    for genes_str in df_sorted['Genes']:
        all_genes.update(str(genes_str).split(';'))
    unique_genes = sorted(list(all_genes))

    # Initialize the plot
    set_style()
    fig, (ax_l, ax_r, ax_m) = plt.subplots(
        1, 3, 
        figsize=(11, len(df_sorted) * 0.3 + 1.5), 
        sharey=True,
        gridspec_kw={'width_ratios': [2, 2, 1.2], 'wspace': 0} 
    )


    # Overlap histogram
    ax_l.barh(range(len(df_sorted)), -df_sorted['Overlap'], 
              color="#E64B35", edgecolor='black', linewidth=0.5, height=0.7, zorder=3)
    ax_l.set_xlabel('Overlap Proportion', fontsize=12)
    ax_l.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "" if x == 0 else f"{abs(x):.2f}"))


    # -log10(FDR) plot
    ax_r.barh(range(len(df_sorted)), df_sorted['log10(FDR)'], 
              color="#3C5488", edgecolor='black', linewidth=0.5, height=0.7, zorder=3)
    ax_r.set_xlabel(r'$-\log_{10}(\text{FDR})$', fontsize=12)


    # Y-labels
    ax_l.set_yticks(range(len(df_sorted)))
    ax_l.set_yticklabels(df_sorted['Term'], ha='right', fontsize=9)
    
    # Hide the shared spine in the middle
    ax_l.spines['right'].set_visible(True)
    ax_r.spines['left'].set_visible(True)
    
    # Ensure a vertical line exists at x=0 for both
    ax_l.axvline(0, color='black', linewidth=0.5, zorder=4)
    ax_r.axvline(0, color='black', linewidth=0.5, zorder=4)


    for ax in [ax_l, ax_r, ax_m]:
        for spine in ['top', 'right', 'left', 'bottom']:
            if ax == ax_m and spine == 'top': continue
            if spine != 'bottom': ax.spines[spine].set_visible(False)
        
        if ax != ax_m:
            ax.grid(axis='x', linestyle='--', alpha=0.3, zorder=0)
        
        ax.tick_params(axis='y', left=False)


    plt.tight_layout()
    plt.savefig(f'{output_dir}/mirror_enrichment_plot.pdf', dpi=300, bbox_inches='tight')



if __name__ == "__main__":

    gtex_med_TPM_path = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_expression/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct"
    gtex_pleio_res_path = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/gtex_res/gtex_susie_pleio_snps.tsv"
    tissues_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_v11_Susie"

    output_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/plots"

    std_exp_data = std_expression(gtex_med_TPM_path, gtex_pleio_res_path, tissues_dir)
    plot_heatmap(std_exp_data, output_dir)

    enrichment_df = hallmark_enrichment_analysis(std_exp_data, output_dir)
    plot_publication_mirror(enrichment_df, output_dir)