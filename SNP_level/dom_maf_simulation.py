import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt

# Nature Genetics style settings
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['legend.frameon'] = False

# Parameters
a = 0
d = 0.0012
h2 = 0.05

mafs = np.linspace(0.1, 0.5, 21)
Ns = [10000, 50000, 100000, 200000, 350000]
iterations = 20

results = []

for N in Ns:
    for maf in mafs:
        p = maf
        q = 1 - p
        prob = [q**2, 2*p*q, p**2]
        
        pvals = []
        for _ in range(iterations):
            # 1. Simulate Genotypes
            g = np.random.choice([0, 1, 2], size=N, p=prob, replace=True)
            
            # 2. Simulate Genetic Value
            y_gen = np.zeros(N)
            y_gen[g == 1] = a + d
            y_gen[g == 2] = 2 * a
            
            # 3. Add Environmental Noise to maintain h2 = 0.05
            var_g = np.var(y_gen)
            if var_g == 0: var_g = 1e-10
            var_e = var_g * (1 - h2) / h2
            y = y_gen + np.random.normal(0, np.sqrt(var_e), N)
            
            # 4. Fit Linear Regression
            slope, intercept, r_value, p_value, std_err = stats.linregress(g, y)
            
            # 5. Robust -log10(p-value) using the t-statistic (prevents underflow to 0)
            t_stat = slope / std_err if std_err > 0 else 0
            if t_stat == 0:
                log10_p = 0
            else:
                log_p_base_e = stats.t.logsf(np.abs(t_stat), df=N-2)
                log10_p = - (np.log10(2) + log_p_base_e / np.log(10))
            
            pvals.append(log10_p)
            
        results.append({
            'N': N,
            'MAF': maf,
            'logP': np.mean(pvals)
        })

df = pd.DataFrame(results)

# 6. Plot the results
plt.figure(figsize=(9, 6))

colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(Ns)))

for i, N in enumerate(Ns):
    subset = df[df['N'] == N]
    plt.plot(subset['MAF'], subset['logP'], marker='o', markersize=6, 
             linewidth=2, label=f'N = {N:,}', color=colors[i])

# Genome-wide significance line
gwas_sig = -np.log10((5e-8)/1060)
plt.axhline(gwas_sig, color='red', linestyle='--', linewidth=1.5, zorder=0)
plt.gca().text(0.15, 0.1, r'$p < 4.72 \times 10^{-11}$', transform=plt.gca().transAxes, fontsize=7, ha='right', color='red')

y_max = df['logP'].max()

plt.xlabel('Minor Allele Frequency (MAF)')
plt.ylabel(r'$-\log_{10}(P)$ of Additive Model')
plt.suptitle('Power of Additive Model for Detecting Dominance Effect\nAcross MAF and Sample Size', fontsize=16)
plt.title(f'Underlying Model: a=0, {d}, h2=0.05', fontsize=12, loc="left")
plt.legend(title='Sample Size (N)', loc='upper right')

plt.tight_layout()
plt.savefig('/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/results/dom_simulation2.png', 
            dpi=600, bbox_inches='tight')