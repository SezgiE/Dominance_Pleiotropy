import os
import textwrap
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles

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


def plot_intersection(all_genes_dir, output_dir):

    all_genes_df = pd.read_csv(all_genes_dir, sep='\t')

    eqtl_genes = set(all_genes_df["gene_id_eqtl"].dropna())
    pos_genes = set(all_genes_df["gene_id_pos"].dropna())

    intersected_genes = eqtl_genes.intersection(pos_genes)
    eqtl_only_genes = eqtl_genes - pos_genes
    pos_only_genes = pos_genes - eqtl_genes

    gene_to_cat = {}
    for col in ["gene_id_eqtl", "gene_id_pos"]:
        subset = all_genes_df[[col, "category"]].dropna()
        gene_to_cat.update(dict(zip(subset[col], subset["category"])))


    int_cats = pd.Series([gene_to_cat[g] for g in intersected_genes if g in gene_to_cat]).value_counts()
    eqtl_only_cats = pd.Series([gene_to_cat[g] for g in eqtl_only_genes if g in gene_to_cat]).value_counts()
    pos_only_cats = pd.Series([gene_to_cat[g] for g in pos_only_genes if g in gene_to_cat]).value_counts()

   
    # Combine into a dataframe and convert to percentages for a clean comparison
    cat_df = pd.DataFrame({'Intersected': int_cats, 'eQTL-only': eqtl_only_cats, 'Positional-only': pos_only_cats}).fillna(0)
    cat_df.index = [textwrap.fill(str(cat), width=25) for cat in cat_df.index]
    cat_df_pct = cat_df.div(cat_df.sum(axis=0), axis=1) * 100

    set_style()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), gridspec_kw={'width_ratios': [2, 2.5], 'wspace': 0.3})
    
    # Panel A: Venn Diagram
    v = venn2([eqtl_genes, pos_genes], 
              set_labels=('eQTL Mapped Genes', 'Positionally Mapped Genes'),
              set_colors=('#E64B35', '#4DBBD5'), 
              alpha=0.8, 
              ax=axes[0])
    current_x, current_y = v.set_labels[1].get_position()
    v.set_labels[1].set_position((current_x - 0.3, current_y - 0.05))
    
    c = venn2_circles([eqtl_genes, pos_genes], 
                      linestyle='solid', linewidth=1.0, color='black', ax=axes[0])
    axes[0].set_title('A.', loc='left', fontweight='bold', fontsize=14, x=-0.1)

    # Panel B: Category Distribution (Stacked Bar)
    npg_palette = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', 
                   '#8491B4', '#91D1C2', '#DC0000', '#7E6148', '#B09C85']
    
    # Transpose so x-axis represents the Groups, and stacks represent Categories
    cat_df_pct.T.plot(kind='bar', stacked=True, ax=axes[1], color=npg_palette[:len(cat_df_pct)], edgecolor='black', linewidth=0.5)
    
    axes[1].set_title('B.', loc='left', fontweight='bold', fontsize=14, x=-0.1)
    axes[1].set_ylabel('Percentage of genes (%)')
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=0)
    
    # Clean up axes and move legend outside the plot area
    axes[1].legend(title='Gene Category', bbox_to_anchor=(0.98, 1), loc='upper left', frameon=False)

    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, 'gene_intersection_panels.pdf')
    plt.savefig(out_path, format='pdf', bbox_inches='tight', transparent=True)
    plt.close()


if __name__ == "__main__":

    all_genes_path = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/genes_all/genes_all.tsv"
    output_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/plots"

    plot_intersection(all_genes_path, output_dir)