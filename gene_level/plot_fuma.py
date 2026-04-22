import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import plotly.graph_objects as go


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


def network_plot(fuma_dir, output_dir):
    # Load positional mapping from genes.txt
    genes_df = pd.read_csv(f"{fuma_dir}/genes.txt", sep='\t')

    # Load eQTL mapping from eqtl.txt
    eqtl_df = pd.read_csv(f"{fuma_dir}/eqtl.txt", sep='\t')

    # ==========================================
    # 2. Filter for the "Hub" Genes (Avoid the Hairball)
    # ==========================================
    # Count how many total SNPs (Positional + eQTL) map to each gene
    gene_hits = {}

    # Count Positional
    for idx, row in genes_df.dropna(subset=['posMapSNPs']).iterrows():
        snps = str(row['posMapSNPs']).split(';')
        gene_hits[row['symbol']] = gene_hits.get(row['symbol'], 0) + len(snps)

    # Count eQTL
    for idx, row in eqtl_df.iterrows():
        gene_hits[row['symbol']] = gene_hits.get(row['symbol'], 0) + 1

    # Get the top 25 most highly connected genes
    top_genes = sorted(gene_hits, key=gene_hits.get, reverse=True)[:25]

    # ==========================================
    # 3. Build the Network Graph
    # ==========================================
    G = nx.Graph()

    # Add Edges for Positional Mapping (Solid Grey Lines)
    for idx, row in genes_df.dropna(subset=['posMapSNPs']).iterrows():
        gene = row['symbol']
        if gene in top_genes:
            snps = str(row['posMapSNPs']).split(';')
            for snp in snps:
                G.add_node(snp, bipartite=0, type='SNP')
                G.add_node(gene, bipartite=1, type='Gene')
                G.add_edge(snp, gene, weight=1.5, type='Positional', color='#B0B0B0', style='solid')

    # Add Edges for eQTL Mapping (Dashed Colored Lines)
    for idx, row in eqtl_df.iterrows():
        gene = row['symbol']
        snp = row['uniqID']
        if gene in top_genes:
            G.add_node(snp, bipartite=0, type='SNP')
            G.add_node(gene, bipartite=1, type='Gene')
            # Using a bright Tangerine/Red for eQTLs to highlight functional evidence
            G.add_edge(snp, gene, weight=2.0, type='eQTL', color='#FF9F1C', style='dashed')

    # ==========================================
    # 4. Nature Genetics Style Visualization
    # ==========================================
    set_style()  # Apply the custom style
    fig, ax = plt.subplots(figsize=(12, 14))

    # Separate nodes for the bipartite layout
    snp_nodes = {n for n, d in G.nodes(data=True) if d['type'] == 'SNP'}
    gene_nodes = {n for n, d in G.nodes(data=True) if d['type'] == 'Gene'}

    # Create the Bipartite Layout (SNPs on left, Genes on right)
    pos = nx.bipartite_layout(G, snp_nodes, align='vertical')

    # Extract edge attributes for plotting
    edges = G.edges(data=True)
    edge_colors = [d['color'] for u, v, d in edges]
    edge_styles = [d['style'] for u, v, d in edges]
    edge_weights = [d['weight'] for u, v, d in edges]

    # Draw Positional Edges first (background)
    pos_edges = [(u, v) for u, v, d in edges if d['type'] == 'Positional']
    nx.draw_networkx_edges(G, pos, edgelist=pos_edges, width=1.5, alpha=0.6, edge_color='#B0B0B0', style='solid', ax=ax)

    # Draw eQTL Edges on top (foreground)
    eqtl_edges = [(u, v) for u, v, d in edges if d['type'] == 'eQTL']
    nx.draw_networkx_edges(G, pos, edgelist=eqtl_edges, width=2.0, alpha=0.8, edge_color='#FF9F1C', style='dashed', ax=ax)

    # Draw Nodes
    # SNPs: Small Royal Blue dots
    nx.draw_networkx_nodes(G, pos, nodelist=snp_nodes, node_size=60, node_color='#4361EE', edgecolors='black', linewidths=0.5, ax=ax)
    # Genes: Larger dark dots
    nx.draw_networkx_nodes(G, pos, nodelist=gene_nodes, node_size=120, node_color='#2B2D42', edgecolors='black', linewidths=1.0, ax=ax)

    # Add Labels (Only labeling the Genes to keep it clean)
    gene_labels = {n: n for n in gene_nodes}
    # Shift the gene labels slightly to the right so they don't overlap the node
    label_pos = {k: [v[0] + 0.05, v[1]] for k, v in pos.items() if k in gene_nodes}
    nx.draw_networkx_labels(G, label_pos, labels=gene_labels, font_size=11, font_weight='bold', font_family='sans-serif', horizontalalignment='left', ax=ax)

    # Add a custom legend
    import matplotlib.lines as mlines
    legend_elements = [
        mlines.Line2D([0], [0], marker='o', color='w', label='Pleiotropic SNP', markerfacecolor='#4361EE', markeredgecolor='black', markersize=8),
        mlines.Line2D([0], [0], marker='o', color='w', label='Mapped Gene', markerfacecolor='#2B2D42', markeredgecolor='black', markersize=10),
        mlines.Line2D([0], [0], color='#B0B0B0', lw=2, linestyle='solid', label='Positional Mapping (< 10kb)'),
        mlines.Line2D([0], [0], color='#FF9F1C', lw=2, linestyle='dashed', label='eQTL Mapping (GTEx v8)')
    ]
    ax.legend(handles=legend_elements, loc='upper right', frameon=False, fontsize=11, bbox_to_anchor=(1.15, 1))

    # Clean up axes
    ax.axis('off')
    plt.title('Pleiotropic Hubs: SNP-to-Gene Mapping Network', fontsize=16, fontweight='bold', pad=20)
    plt.tight_layout()

    # Save as high-res PDF for publication
    saving_path = f"{output_dir}/pleiotropy_network.pdf"
    plt.savefig(saving_path, dpi=600, bbox_inches='tight')





if __name__ == "__main__":

    output_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/plots"
    fuma_result_dir = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/fuma_out"

    #network_plot(fuma_result_dir, output_dir)