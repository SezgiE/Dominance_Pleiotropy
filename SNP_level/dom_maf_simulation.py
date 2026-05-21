import json
import numpy as np
import seaborn as sns
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Patch


def get_emp_dvals(all_snps_path):

    df = pd.read_csv(
        all_snps_path, sep="\t", compression="gzip", usecols=["variant", "dom_betas"]
    )

    parsed_rows = []
    for row in df.itertuples():
        betas_dict = json.loads(row.dom_betas)
        for pheno_code, beta_val in betas_dict.items():
            parsed_rows.append(
                {"variant": row.variant, "phenotype_code": pheno_code, "dom_betas": abs(beta_val)}
            )

    df_long = pd.DataFrame(parsed_rows)

    min_val = df_long["dom_betas"].min()
    q1_val = df_long["dom_betas"].quantile(0.25)
    q3_val = df_long["dom_betas"].quantile(0.75)
    max_val = df_long["dom_betas"].max()

    beta_summary_list = [min_val, q1_val, q3_val, max_val]

    return beta_summary_list


all_snps_path = "/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/significant_SNPs/all_sig_SNPs.tsv.gz"

# Parameters
a = 0
#d_values = get_emp_dvals(all_snps_path)
d_values = np.linspace(0, 0.5, 11)

mafs = np.linspace(0.1, 0.5, 21)
Ns = [10000, 50000, 100000, 500000]
iterations = 100


def dom_maf_simulation():

    results = []

    for d in d_values:
        for N in Ns:
            for maf in mafs:
                p = maf
                q = 1 - p
                prob = [q**2, 2*p*q, p**2]
                
                v_g_theoretical = (2*p*q*(a+d*(q-p))**2) + ((2*p*q*d)**2)
                var_e = 1 - v_g_theoretical
                
                pvals = []
                h2_vals = []
                for _ in range(iterations):
                    print(f"Running simulation for d={d}, N={N}, MAF={maf:.2f} (Iteration {_+1}/{iterations})")
                    # Simulate Genotypes
                    g = np.random.choice([0, 1, 2], size=N, p=prob, replace=True)
                    
                    # Simulate Genetic Value
                    y_gen = np.zeros(N)
                    y_gen[g == 1] = a + d
                    y_gen[g == 2] = 2 * a

                    # Add Environmental Noise 
                    y_raw = y_gen + np.random.normal(0, np.sqrt(var_e), N)
                    y = (y_raw - np.mean(y_raw)) / np.std(y_raw)
                    
                    var_g = np.var(y_gen)
                    if var_g == 0: var_g = 1e-10
                    empirical_H2 = var_g / (var_g + var_e)  
                    h2_vals.append(empirical_H2)
                    
                    # Fit Linear Regression
                    slope, intercept, r_value, p_value, std_err = stats.linregress(g, y)

                    if q == p: 
                        std_d = 0
                    else:
                        std_d = slope / (q - p)
                    
                    #  -log10(p-value) using the t-statistic
                    t_stat = slope / std_err if std_err > 0 else 0
                    if t_stat == 0:
                        log10_p = 0
                    else:
                        log_p_base_e = stats.norm.logsf(np.abs(t_stat))
                        log10_p = - (np.log10(2) + log_p_base_e / np.log(10))
                        
                        # if it hits inf, cap it at 300
                        if np.isinf(log10_p):
                            log10_p = 300
                    
                    pvals.append(log10_p)
                    
                results.append({
                    "d": d,
                    "std_d": std_d,
                    'N': N,
                    'MAF': maf,
                    'logP': np.mean(pvals),
                    'Emprical_H2': np.mean(h2_vals)    
                })

    df = pd.DataFrame(results)
    df.to_csv('/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/results/dom_simulation_results.tsv', index=False, sep="\t")


# plotting the results
def set_style():
    """Hardcodes matplotlib parameters to format the plot."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 12,
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'axes.linewidth': 1.0,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'pdf.fonttype': 42,
        'ps.fonttype': 42
    })


def plot_3d(df):
    set_style()
    df = df[df['logP'] < 250]
    
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection='3d')
    
    # Select a subset of Ns to prevent complete occlusion
    Ns_to_plot = [10000, 50000, 100000, 500000]
    colors = ['#440154', '#31688e', '#35b779', '#fde725']
    
    for N, color in zip(Ns_to_plot, colors):
        sub_df = df[df['N'] == N]
        
        # Pivot data to create a 2D matrix for the Z-axis
        Z_matrix = sub_df.pivot(index='d', columns='MAF', values='logP').values
        
        # Create X and Y meshgrids
        X, Y = np.meshgrid(sub_df['MAF'].unique(), sub_df['d'].unique())
        
        # Plot the surface
        ax.plot_surface(X, Y, Z_matrix, color=color, alpha=0.6, 
                        edgecolor='k', linewidth=0.3, antialiased=True)
        
    ax.tick_params(axis='x', pad=2)  # Distance for MAF numbers
    ax.tick_params(axis='y', pad=2)  # Distance for d numbers
    ax.tick_params(axis='z', pad=2)  # Distance for logP numbers

    ax.set_xlabel('Minor Allele Frequency (MAF)', labelpad=7)
    ax.set_ylabel('Dominance Effect Size (d)', labelpad=7)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel('Mean $-log_{10}(P)$', labelpad=6, rotation=90)
    
    # Add a custom legend
    legend_elements = [Patch(facecolor=c, edgecolor='k', alpha=0.6, label=f'N = {n:,}') 
                       for c, n in zip(colors, Ns_to_plot)]
    
    ax.legend(handles=legend_elements, loc='upper center', 
              bbox_to_anchor=(0.5, -0.18), frameon=False, 
              ncol=len(Ns_to_plot),
              fontsize=10)
    
    #(Elevation, Azimuth)
    ax.view_init(elev=25, azim=-135)
    ax.set_box_aspect(None, zoom=1.3)
    
    plt.tight_layout()
    plt.savefig('/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/results/dom_sim_plot.pdf')


def plot_dom_maf_simulation(df):
    
    set_style()
    fig, axes = plt.subplots(2, 2, figsize=(10, 7)) 
    
    # Flatten the 2x2 array into a 1D list of 4 axes so it matches the 4 d_values
    axes_flat = axes.flatten()

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(Ns)))
    panels = ['A', 'B', "C", "D"]

    for idx, (ax, current_d) in enumerate(zip(axes_flat, d_values)):
        subset_d = df[df['d'] == current_d]
        
        for i, N in enumerate(Ns):
            subset = subset_d[subset_d['N'] == N]
            ax.plot(subset['MAF'], subset['logP'], marker='o', markersize=4, 
                    linewidth=1.5, label=f'N = {N:,}', color=colors[i])
        
        # Genome-wide significance line
        gwas_sig = -np.log10((5e-8)/1060)
        ax.axhline(gwas_sig, color='#7E818D', linestyle='--', linewidth=1, zorder=0)
        
        # Adjusted text position to prevent overlap with data
        ax.text(0.8, 0.24, r'$p < 4.72 \times 10^{-11}$', transform=ax.transAxes, 
            fontsize=8, ha='left', color='black')

        # Nature Genetics style: outer edges only for axis labels
        if idx in [2, 3]:
            ax.set_xlabel('Minor Allele Frequency (MAF)')
        if idx in [0, 2]:
            ax.set_ylabel(r'Mean $-\log_{10}(P)$ based on Additive Model')
        
        ax.set_title('')
    
        ax.text(-0.1, 1.05, panels[idx], transform=ax.transAxes, 
                fontsize=12, fontweight='bold', va='bottom')

    # Place legend and text box outside the grid on the top right panel (Panel B / index 1)
    axes_flat[1].legend(title='Sample Size (N)', loc='upper left', bbox_to_anchor=(1.05, 1), frameon=False)

    # Expanded to include all 4 summary statistics (rounded to 4 decimal places)
    axes_flat[3].text(1.1, -0.15, 
                f'$\\mathbf{{Panel\\ A}}$ (Min):\n  $d={d_values[0]:.4f}$\n\n'
                f'$\\mathbf{{Panel\\ B}}$ (Q1):\n  $d={d_values[1]:.4f}$\n\n'
                f'$\\mathbf{{Panel\\ C}}$ (Q3):\n  $d={d_values[2]:.4f}$\n\n'
                f'$\\mathbf{{Panel\\ D}}$ (Max):\n  $d={d_values[3]:.4f}$', 
                transform=axes_flat[1].transAxes, fontsize=9, va='top', ha='left', color='#333333',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#333333', alpha=0.8))

    plt.tight_layout()
    plt.savefig('/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/results/dom_simulation.pdf', 
                dpi=600, bbox_inches='tight')


if __name__ == "__main__":
    # Run the simulation and save results
    dom_maf_simulation()

    # Load results and plot
    df_results = pd.read_csv('/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/results/dom_simulation_results.tsv', sep="\t")
    #plot_dom_maf_simulation(df_results)
    plot_3d(df_results)