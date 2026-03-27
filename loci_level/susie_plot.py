import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import norm

# ==========================================
# 1. NATURE GENETICS PLOT SETTINGS
# ==========================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'pdf.fonttype': 42,       
    'ps.fonttype': 42,
    'axes.linewidth': 1.0,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'axes.spines.top': False, 
    'axes.spines.right': False, 
    'xtick.direction': 'out',
    'ytick.direction': 'out'
})

def calculate_pval_from_z(z_scores):
    """Calculates strict -log10 P-values directly from Z-scores."""
    # norm.sf(abs(Z)) gives the one-tailed p-value, multiply by 2 for two-tailed
    p_values = 2 * norm.sf(np.abs(z_scores))
    p_values = np.clip(p_values, a_min=1e-400, a_max=1.0) # Prevent log10(0) crash
    return -np.log10(p_values)

def get_ld_color(r2):
    """Maps R^2 to the classic LocusZoom color palette."""
    if pd.isna(r2): return '#7F7F7F' # Grey
    if r2 >= 0.8: return '#FF0000'   # Red
    if r2 >= 0.6: return '#FFA500'   # Orange
    if r2 >= 0.4: return '#00FF00'   # Green
    if r2 >= 0.2: return '#87CEFA'   # Light Blue
    return '#000080'                 # Navy

def plot_stacked_locus(df, out_pdf):
    """
    Generates a publication-ready 2-panel locus plot based on your exact structure.
    """
    # 1. Prepare the math
    # Calculate -log10(P) using your dom_z_score column
    if 'log10_pval' not in df.columns:
        df['log10_pval'] = calculate_pval_from_z(df['dom_z_score'])
        
    df['pos_Mb'] = df['pos'] / 1_000_000
    
    # 2. Setup the Canvas
    fig = plt.figure(figsize=(8, 7))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1.5, 1], hspace=0.15)
    ax_gwas = fig.add_subplot(gs[0])
    ax_pip = fig.add_subplot(gs[1], sharex=ax_gwas)

    # ==========================================
    # TOP PANEL: GWAS (-log10 P-value)
    # ==========================================
    colors_gwas = df['lead_r2'].apply(get_ld_color)
    
    # Background noise SNPs
    bg_mask = df['lead_r2'] < 0.2
    ax_gwas.scatter(df.loc[bg_mask, 'pos_Mb'], df.loc[bg_mask, 'log10_pval'], 
                    c=colors_gwas[bg_mask], s=30, alpha=0.5, edgecolors='none')
    
    # Foreground high LD SNPs
    fg_mask = df['lead_r2'] >= 0.2
    ax_gwas.scatter(df.loc[fg_mask, 'pos_Mb'], df.loc[fg_mask, 'log10_pval'], 
                    c=colors_gwas[fg_mask], s=40, alpha=0.9, edgecolors='black', linewidths=0.5)

    # Lead SNP Diamond
    lead_idx = df['PIP'].idxmax()
    ax_gwas.scatter(df.loc[lead_idx, 'pos_Mb'], df.loc[lead_idx, 'log10_pval'], 
                    c='#9400D3', marker='D', s=80, edgecolors='black', linewidths=1, zorder=5)

    ax_gwas.set_ylabel(r'$-\log_{10}(P)$')
    ax_gwas.tick_params(labelbottom=False)

    # ==========================================
    # BOTTOM PANEL: SuSiE PIPs & Credible Sets
    # ==========================================
    # Background SNPs (CS == 0 or NA)
    bg_cs = (df['CS'] == 0) | (df['CS'].isna())
    ax_pip.scatter(df.loc[bg_cs, 'pos_Mb'], df.loc[bg_cs, 'PIP'], 
                   c='#D3D3D3', s=30, alpha=0.5, edgecolors='none', label='Background')

    # Credible Sets
    cs_colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00']
    unique_cs = sorted([cs for cs in df['CS'].unique() if pd.notna(cs) and cs > 0])
    
    for i, cs in enumerate(unique_cs):
        cs_mask = df['CS'] == cs
        color = cs_colors[i % len(cs_colors)]
        ax_pip.scatter(df.loc[cs_mask, 'pos_Mb'], df.loc[cs_mask, 'PIP'], 
                       c=color, s=50, alpha=0.9, edgecolors='black', linewidths=0.5, 
                       label=f'Credible Set {int(cs)}')

    # Lead PIP Diamond in bottom plot
    ax_pip.scatter(df.loc[lead_idx, 'pos_Mb'], df.loc[lead_idx, 'PIP'], 
                   c='#9400D3', marker='D', s=80, edgecolors='black', linewidths=1, zorder=5)

    ax_pip.set_ylabel('Posterior Inclusion\nProbability (PIP)')
    ax_pip.set_xlabel(f'Chromosome {df["chr"].iloc[0]} Position (Mb)')
    ax_pip.set_ylim(-0.05, 1.05)
    ax_pip.legend(loc='upper right', frameon=False, prop={'size': 9})

    # Save output
    plt.savefig(out_pdf, format='pdf', dpi=600, bbox_inches='tight')
    plt.close()
    print(f"Saved publication plot: {out_pdf}")

# Example usage:
df = pd.read_csv("/Users/sezgi/Documents/dominance_pleiotropy/loci_level/susie_results/1747_2_susie_res.tsv", sep='\t')
plot_stacked_locus(df, "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/susie_results/Nature_Genetics_Locus1.pdf")