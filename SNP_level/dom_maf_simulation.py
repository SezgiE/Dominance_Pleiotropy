import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt

"""
    Max absolute dominance beta across all traits: ('30790_irnt', '6:160985526:G:A', 0.39655)
    Min absolute dominance beta across all traits: ('M13_FIBROBLASTIC', '7:38007174:G:A', 0.00101925)
    Median SNP heritability across all traits: 0.022
    Median dominance heritability across all traits: 0.0002367807
"""

# Parameters
a = 0
d_values = [0.001, 0.4] 
h2_d = 0.0002              # Median dominance heritability

mafs = np.linspace(0.1, 0.5, 21)
Ns = [10000, 50000, 100000, 200000, 350000]
iterations = 1000


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


def dom_maf_simulation():
    # Parameters

    results = []

    for d in d_values:
        for N in Ns:
            for maf in mafs:
                p = maf
                q = 1 - p
                prob = [q**2, 2*p*q, p**2]
                
                # Calculate Theoretical Broad-Sense Heritability (H2) based on MAF and median h2_d
                x = ((q - p)**2) / (2 * p * q)
                H2 = h2_d * (1 + x)
                
                pvals = []
                for _ in range(iterations):
                    print(f"Running simulation for d={d}, N={N}, MAF={maf:.2f} (Iteration {_+1}/{iterations})")
                    # Simulate Genotypes
                    g = np.random.choice([0, 1, 2], size=N, p=prob, replace=True)
                    
                    # Simulate Genetic Value
                    y_gen = np.zeros(N)
                    y_gen[g == 1] = a + d
                    y_gen[g == 2] = 2 * a

                    # Add Environmental Noise 
                    var_g = np.var(y_gen)
                    if var_g == 0: var_g = 1e-10
                    
                    var_e = var_g * (1 - H2) / H2

                    # (Safeguard in case var_g is somehow > 1)
                    if var_e <= 0: var_e = 1e-10 

                    y = y_gen + np.random.normal(0, np.sqrt(var_e), N)
                    
                    # Fit Linear Regression
                    slope, intercept, r_value, p_value, std_err = stats.linregress(g, y)
                    
                    #  -log10(p-value) using the t-statistic
                    t_stat = slope / std_err if std_err > 0 else 0
                    if t_stat == 0:
                        log10_p = 0
                    else:
                        log_p_base_e = stats.t.logsf(np.abs(t_stat), df=N-2)
                        log10_p = - (np.log10(2) + log_p_base_e / np.log(10))
                    
                    pvals.append(log10_p)
                    
                results.append({
                    "d": d,
                    'N': N,
                    'MAF': maf,
                    'logP': np.mean(pvals),
                    'Theoretical_H2': H2,     # Your formula's expectation
                    'Empirical_Var_G': var_g  # What the simulation actually generated
                })

    df = pd.DataFrame(results)
    df.to_csv('/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/results/dom_simulation_results.tsv', index=False, sep="\t")


def plot_dom_maf_simulation(df):
    
    set_style()
    fig, axes = plt.subplots(1, 2, figsize=(10, 5)) # Adjusted for dual panels

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(Ns)))
    panels = ['A', 'B']

    for idx, (ax, current_d) in enumerate(zip(axes, d_values)):
        subset_d = df[df['d'] == current_d]
        
        for i, N in enumerate(Ns):
            subset = subset_d[subset_d['N'] == N]
            ax.plot(subset['MAF'], subset['logP'], marker='o', markersize=4, 
                    linewidth=1.5, label=f'N = {N:,}', color=colors[i])
        
        # Genome-wide significance line
        gwas_sig = -np.log10((5e-8)/1060)
        ax.axhline(gwas_sig, color='#7E818D', linestyle='--', linewidth=1, zorder=0)
        
        # Adjusted text position to prevent overlap with data
        ax.text(0.8, 0.24, r'$p \lesssim 4.72 \times 10^{-11}$', transform=ax.transAxes, 
            fontsize=8, ha='left', color='black')

        
        ax.set_xlabel('Minor Allele Frequency (MAF)')
        if idx == 0:
            ax.set_ylabel(r'$-\log_{10}(P)$ of Additive Model')
        
        ax.set_title('')
    
        ax.text(-0.1, 1.05, panels[idx], transform=ax.transAxes, 
                fontsize=12, fontweight='bold', va='bottom')

    # Only place legend on the right panel to save space
    axes[1].legend(title='Sample Size (N)', loc='upper left', bbox_to_anchor=(0.8, 1), frameon=False)

    axes[1].text(1.00, 0.55, f'$\\mathbf{{Panel\\ A}}$:\n  $d={d_values[0]}$\n\n$\\mathbf{{Panel\\ B}}$:\n  $d={d_values[1]}$', 
                transform=axes[1].transAxes, fontsize=9, va='top', ha='left', color='#333333',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#333333', alpha=0.8))

    plt.tight_layout()
    plt.savefig('/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/results/dom_simulation.pdf', 
                dpi=600, bbox_inches='tight')


if __name__ == "__main__":
    # Run the simulation and save results
    dom_maf_simulation()

    # Load results and plot
    df_results = pd.read_csv('/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/results/dom_simulation_results.tsv', sep="\t")
    plot_dom_maf_simulation(df_results)