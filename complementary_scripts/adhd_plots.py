from ast import For
import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pyreadr

def set_style():
    """Hardcodes matplotlib parameters for Nature Genetics standards."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8,
        'axes.labelsize': 8,
        'axes.titlesize': 9,
        'xtick.labelsize': 5,
        'ytick.labelsize': 7,
        'legend.fontsize': 8,
        'legend.title_fontsize': 9,
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'pdf.fonttype': 42,
        'ps.fonttype': 42
    })


def set_style2():
    """Hardcodes matplotlib parameters for Nature Genetics standards."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8,
        'axes.labelsize': 10,
        'axes.titlesize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'legend.title_fontsize': 12,
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'pdf.fonttype': 42,
        'ps.fonttype': 42
    })


def plot_item_level_trends(trend_res_path, output_dir):
    file_list = glob.glob(os.path.join(trend_res_path, "*.rds"))
    
    plot_data_list = []
    
    # ------------------ DATA PROCESSING ------------------
    for file in file_list:
        filename = os.path.basename(file)
        q_name = re.sub(r'^(q\d+).*', r'\1', filename)
        
        # Read .rds file
        result = pyreadr.read_r(file)
        df = result[None] # pyreadr returns a dictionary of objects
        
        # Drop NAs for complete cases
        cols_to_check = ["FID", "PID", "mult_id_fam", "sex", "yob", "zyg", "assigned_age"]
        df = df.dropna(subset=[c for c in cols_to_check if c in df.columns]).copy()
        
        # Create 6-year bins for birth year
        yob_min = int(df['yob'].min())
        yob_max = int(df['yob'].max())
        
        # Define breaks and labels
        breaks = list(range(yob_min, yob_max + 7, 6))
        labels = [f"{breaks[i]}-{breaks[i+1]-1}" for i in range(len(breaks)-1)]
        
        df['yob_bin'] = pd.cut(df['yob'], bins=breaks, labels=labels, right=False, include_lowest=True)
        
        # For plotting: if group n < 200, bin it as "2004-2009"
        df['group_count'] = df.groupby(['yob_bin', 'respondent', 'assigned_age'], observed=True)['yob'].transform('count')
        df['yob_bin_plot'] = np.where(df['group_count'] < 200, "2004-2009", df['yob_bin'].astype(str))
        
        # Descriptive stats
        desc_stats = df.groupby(['yob_bin_plot', 'respondent', 'assigned_age'], observed=True).agg(
            mean_response=('response', 'mean'),
            sd_response=('response', 'std'),
            n=('response', 'count')
        ).reset_index()
        
        desc_stats['sem_response'] = desc_stats['sd_response'] / np.sqrt(desc_stats['n'])
        desc_stats['q_name'] = q_name
        
        plot_data_list.append(desc_stats)
        print(f"Processed: {q_name}")
        
    # Combine all items into one dataframe
    all_plot_data = pd.concat(plot_data_list, ignore_index=True)
    
    # Sort items numerically for the grid
    items = all_plot_data['q_name'].unique()
    items = sorted(items, key=lambda x: int(re.search(r'\d+', x).group()))
    
    # ------------------ PLOTTING ------------------
    set_style()
    
    # Grid math: 4 columns, dynamic rows (usually 3 for 12 items)
    n_cols = 4
    n_rows = int(np.ceil(len(items) / n_cols))
    
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(9, 9*4.5/8))
    axes = axes.flatten()
    fig.subplots_adjust(left=0.08, right=0.98, top=0.92, bottom=0.20, hspace=0.7, wspace=0.4)
    
    # NPG Color Palette & Styles
    rater_colors = {"m": "#3C5488", "f": "#F39B7F", "t": "#00A087"}
    age_linestyles = {"7": "-", "10": "--", "12": ":"}
    
    for i, item in enumerate(items):
        ax = axes[i]
        subset = all_plot_data[all_plot_data['q_name'] == item].copy()
        
        # Ensure x-axis categorical sorting
        unique_bins = sorted(subset['yob_bin_plot'].unique())
        
        y_max = round(subset['mean_response'].max(), 2) + 0.05
        
        for (rater, age), grp in subset.groupby(['respondent', 'assigned_age']):
            # Sort by the x-axis categorical order
            grp['yob_bin_plot'] = pd.Categorical(grp['yob_bin_plot'], categories=unique_bins, ordered=True)
            grp = grp.sort_values('yob_bin_plot')
            
            # Line + Points
            ax.plot(grp['yob_bin_plot'].astype(str), grp['mean_response'], 
                    color=rater_colors.get(rater, 'black'), 
                    linestyle=age_linestyles.get(str(int(age)), '-'), 
                    linewidth=1.2, marker='o', markersize=3)
            
        # Axes formatting
        ax.set_title(f"Item {item.replace('q', '')}", fontweight='bold', fontsize=8)
        ax.grid(True, axis='y', color='gray', linestyle='-', linewidth=0.2, alpha=0.5)
        ax.set_ylim(0, y_max)
        ax.set_yticks(np.linspace(0, y_max, 4).round(2))
        ax.tick_params(axis='x', rotation=20)
        
        # Minimalist labels
        if i % n_cols == 0:
            ax.set_ylabel("Mean Response")
        if i >= len(items) - n_cols:
            ax.set_xlabel("Birth Cohort")
            
    # Hide empty subplots
    for j in range(len(items), len(axes)):
        axes[j].set_visible(False)
        
    # Unified custom legends
    custom_lines_rater = [Line2D([0], [0], color=c, lw=1.5) for c in rater_colors.values()]
    custom_lines_age = [Line2D([0], [0], color='black', lw=1.5, ls=ls) for ls in age_linestyles.values()]
    
    l1 = fig.legend(custom_lines_rater, ['Mother', 'Father', 'Teacher'], 
                     loc='upper center', bbox_to_anchor=(0.35, 0.08), ncol=3, frameon=False)
    l2 = fig.legend(custom_lines_age, ['7', '10', '12'], 
                     loc='upper center', bbox_to_anchor=(0.7, 0.08), ncol=3, frameon=False)

    # Save out
    out_path = os.path.join(f"{output_dir}/descriptives_IL", "sex_agg_trend.pdf")
    plt.savefig(out_path, format='pdf', transparent=True)
    plt.close(fig)
    print(f"Plot saved successfully to {out_path}")


def plot_beta_distributions2(trend_res_path, trend_res_path_corr, output_dir):

    df = pd.read_excel(trend_res_path, usecols=["q_name","rater","rated_age","significance", 
                                                "beta_bc", "low_conf", "up_conf", "lower_thresh", "upper_thresh"])
    
    df_corr = pd.read_excel(trend_res_path_corr, usecols=["q_name","rater","rated_age","significance", 
                                                "beta_bc", "low_conf", "up_conf", "lower_thresh", "upper_thresh"])
    
    df = df[df['significance'] == True].copy()
    df = df.drop_duplicates(subset=['q_name', 'rater', 'rated_age']).reset_index(drop=True)
    

    df['err_low'] = df['beta_bc'] - df['low_conf']
    df['err_up'] = df['up_conf'] - df['beta_bc']
    

    df['q_num'] = df['q_name'].str.extract(r'(\d+)').astype(float)
    df = df.sort_values(by=['q_num', 'rated_age', 'rater'])
    items = df['q_name'].unique()
    
    
    set_style2()
    fig, ax = plt.subplots(figsize=(8, 6))

    for i in range(len(items)):
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color='gray', alpha=0.05, zorder=0)
    

    
    rater_colors = {"m": "#3C5488", "f": "#F39B7F", "t": "#00A087"}

    age_linestyles = {7: '-', 10: '--', 12: ':'}  # Solid, Dashed, Dotted
    
    for i, item in enumerate(items):
        subset = df[df['q_name'] == item]
        n_points = len(subset)

        offsets = np.linspace(-0.3, 0.3, n_points) if n_points > 1 else [0]
        
        for j, (_, row) in enumerate(subset.iterrows()):
            color = rater_colors.get(str(row['rater']).lower(), '#333333')
            age_val = int(row['rated_age'])
            
            ls = age_linestyles.get(age_val, '-')
            
            eb = ax.errorbar(row['beta_bc'], i + offsets[j], 
                        xerr=[[row['err_low']], [row['err_up']]],
                        marker='o', linestyle='none', color=color, markersize=3.5, 
                        elinewidth=1.0, capsize=0)
            
            eb[2][0].set_linestyle(ls)
            

    ax.axvline(0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    

    ax.set_yticks(range(len(items)))
    ax.set_yticklabels([f"Item {item.replace('q', '')}" for item in items])
    ax.set_xlabel("Beta Coefficient (95% CI)", labelpad=20)
    

    custom_lines_rater = [Line2D([0], [0], color=c, marker='o', lw=0, markersize=4) for c in rater_colors.values()]
    l1 = ax.legend(custom_lines_rater, ['Mother', 'Father', 'Teacher'], title="Respondent", 
              frameon=False, bbox_to_anchor=(1.02, 1), loc='upper left')
              
    
    custom_lines_age = [Line2D([0], [0], color='gray', linestyle=ls, lw=1.5) for ls in age_linestyles.values()] 
    l2 = ax.legend(custom_lines_age, ['7', '10', '12'], title="Assessment Age", 
              frameon=False, bbox_to_anchor=(1.02, 0.65), loc='upper left')
              
    ax.add_artist(l1) 


    fig.subplots_adjust(left=0.12, right=0.80, top=0.92, bottom=0.15)

    
    out_path = os.path.join(f"{output_dir}/model_results_IL", "beta_distributions_corr.pdf")
    plt.savefig(out_path, format='pdf', transparent=True)
    plt.close(fig)


def plot_beta_distributions(trend_res_path, trend_res_path_corr, output_dir):

    df = pd.read_excel(trend_res_path, usecols=["q_name","rater","rated_age","significance",
                                                "beta_bc", "low_conf", "up_conf", "lower_thresh", "upper_thresh"])
    
    df_corr = pd.read_excel(trend_res_path_corr, usecols=["q_name","rater","rated_age","significance",
                                                          "beta_bc", "low_conf", "up_conf", "lower_thresh", "upper_thresh"])

    # Process Baseline (df)
    df = df[df['significance'] == True].copy()
    df = df.drop_duplicates(subset=['q_name', 'rater', 'rated_age']).reset_index(drop=True)
    df['err_low'] = df['beta_bc'] - df['low_conf']
    df['err_up'] = df['up_conf'] - df['beta_bc']
    df['q_num'] = df['q_name'].str.extract(r'(\d+)').astype(float)

    baseline_sig_combos = df.set_index(['q_name', 'rater', 'rated_age']).index
    
    # Process Corrected (df_corr)
    corr_combos = df_corr.set_index(['q_name', 'rater', 'rated_age']).index
    df_corr = df_corr[(df_corr['significance'] == True) | (corr_combos.isin(baseline_sig_combos))].copy()
    df_corr = df_corr.drop_duplicates(subset=['q_name', 'rater', 'rated_age']).reset_index(drop=True)
    df_corr['err_low'] = df_corr['beta_bc'] - df_corr['low_conf']
    df_corr['err_up'] = df_corr['up_conf'] - df_corr['beta_bc']
    df_corr['q_num'] = df_corr['q_name'].str.extract(r'(\d+)').astype(float)

    # Get unique items from BOTH dataframes to ensure a perfectly shared y-axis
    all_items_df = pd.concat([df[['q_name', 'q_num']], df_corr[['q_name', 'q_num']]]).drop_duplicates()
    all_items_df = all_items_df.sort_values(by=['q_num'])
    items = all_items_df['q_name'].tolist()
    
    
    set_style2()
    # Create 1x2 grid sharing the y-axis
    fig, axes = plt.subplots(1, 2, figsize=(11, 7), sharey=True)
    
    rater_colors = {"m": "#3C5488", "f": "#F39B7F", "t": "#00A087"}
    age_linestyles = {7: '-', 10: '--', 12: ':'}
    
    panels = [
        (axes[0], df, 'Unadjusted', 'A.'),
        (axes[1], df_corr, 'Adjusted for Parental Age', 'B.')
    ]
    
    for ax, panel_df, title, panel_letter in panels:
        # Alternating background stripes
        for i in range(len(items)):
            if i % 2 == 0:
                ax.axhspan(i - 0.5, i + 0.5, color='gray', alpha=0.05, zorder=0)
                
        for i, item in enumerate(items):
            subset = panel_df[panel_df['q_name'] == item]
            n_points = len(subset)
            offsets = np.linspace(-0.35, 0.35, n_points) if n_points > 1 else [0]
            
            for j, (_, row) in enumerate(subset.iterrows()):
                color = rater_colors.get(str(row['rater']).lower(), '#333333')
                age_val = int(row['rated_age'])
                ls = age_linestyles.get(age_val, '-')
                
                is_sig = row.get('significance', True) 
                point_alpha = 1.0 if is_sig else 0.6

                eb = ax.errorbar(row['beta_bc'], i + offsets[j], 
                            xerr=[[row['err_low']], [row['err_up']]],
                            marker='o', linestyle='none', color=color, markersize=3.5, 
                            elinewidth=1.0, capsize=0, alpha=point_alpha)
                
                eb[2][0].set_linestyle(ls)
                
        ax.axvline(0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
        ax.set_title(title, fontweight='bold', fontsize=12)
        ax.set_xlabel("Beta Coefficient (99.86% CI)", labelpad=10)
        
        # Panel lettering
        ax.text(-0.05 if panel_letter == 'B.' else -0.15, 1.05, panel_letter, 
                transform=ax.transAxes, size=12, weight='bold', va='bottom', ha='right')

    # Set y-ticks ONLY on the left axis (axes[0])
    axes[0].set_yticks(range(len(items)))
    axes[0].set_yticklabels([f"Item {item.replace('q', '')}" for item in items])
    axes[1].tick_params(left=False) # Remove inner tick marks between the plots
    axes[1].spines['left'].set_visible(False)
    axes[0].invert_yaxis()
    axes[0].set_xlim(axes[1].get_xlim())

    # Unified legends attached to the figure at the bottom
    custom_lines_rater = [Line2D([0], [0], color=c, marker='o', lw=0, markersize=4) for c in rater_colors.values()]
    l1 = fig.legend(custom_lines_rater, ['Mother', 'Father', 'Teacher'], title="Respondent", 
                    frameon=False, bbox_to_anchor=(0.40, 0.08), loc='upper center', ncol=3,
                    fontsize=10, title_fontsize=12, handletextpad=0.4, columnspacing=1.0)
              
    custom_lines_age = [Line2D([0], [0], color='gray', linestyle=ls, lw=1.5) for ls in age_linestyles.values()] 
    l2 = fig.legend(custom_lines_age, ['7', '10', '12'], title="Assessment Age", 
                    frameon=False, bbox_to_anchor=(0.65, 0.08), loc='upper center', ncol=3,
                    fontsize=10, title_fontsize=12, handletextpad=0.4, columnspacing=1.0)

    # Adjust layout to fit both panels (pushing them together via wspace) and the bottom legend
    fig.subplots_adjust(left=0.1, right=0.98, top=0.9, bottom=0.18, wspace=0.1)

    
    out_path = os.path.join(f"{output_dir}/model_results_IL", "beta_distributions.pdf")
    plt.savefig(out_path, format='pdf', transparent=True)
    plt.close(fig)


def scale_combined_figure(trend_res_path, scale_trend_path, scale_trend_path_corr, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    # ==========================================
    # 1. DATA PROCESSING: TOP ROW (TRENDS)
    # ==========================================
    file_path = os.path.join(trend_res_path, "sumscore_Data.RDS")
    result = pyreadr.read_r(file_path)
    df_trends = result[None] 
    
    # Keep complete cases
    cols_to_check = ["FID", "PID", "mult_id_fam", "sex", "yob", "zyg", "assigned_age"]
    df_trends = df_trends.dropna(subset=[c for c in cols_to_check if c in df_trends.columns]).copy()
    
    # Create 6-year bins for birth year
    yob_min = int(df_trends['yob'].min())
    yob_max = int(df_trends['yob'].max())
    breaks = list(range(yob_min, yob_max + 7, 6))
    labels = [f"{breaks[i]}-{breaks[i+1]-1}" for i in range(len(breaks)-1)]
    
    df_trends['yob_bin'] = pd.cut(df_trends['yob'], bins=breaks, labels=labels, right=False, include_lowest=True)
    
    # Bin groups < 200
    df_trends['group_count'] = df_trends.groupby(['scale', 'yob_bin', 'respondent', 'assigned_age'], observed=True)['yob'].transform('count')
    df_trends['yob_bin_plot'] = np.where(df_trends['group_count'] < 200, "2004-2009", df_trends['yob_bin'].astype(str))
    
    # Descriptive stats for plotting
    desc_stats = df_trends.groupby(['scale', 'yob_bin_plot', 'respondent', 'assigned_age'], observed=True).agg(
        mean_response=('score', 'mean'),
        sd_response=('score', 'std'),
        n=('score', 'count')
    ).reset_index()
    desc_stats['sem_response'] = desc_stats['sd_response'] / np.sqrt(desc_stats['n'])


    # ==========================================
    # 2. DATA PROCESSING: BOTTOM ROW (BETAS)
    # ==========================================
    df = pd.read_excel(scale_trend_path, usecols=["scale","rater","rated_age","significance",
                                                "beta_bc", "low_conf", "up_conf", "lower_thresh", "upper_thresh"])
    
    df_corr = pd.read_excel(scale_trend_path_corr, usecols=["scale","rater", "rated_age", "significance",
                                                          "beta_bc", "low_conf", "up_conf", "lower_thresh", "upper_thresh"])

    # Process Baseline (df)
    df = df[df['significance'] == True].copy()
    df = df.drop_duplicates(subset=['scale', 'rater', 'rated_age']).reset_index(drop=True)
    df['err_low'] = df['beta_bc'] - df['low_conf']
    df['err_up'] = df['up_conf'] - df['beta_bc']

    baseline_sig_combos = df.set_index(['scale', 'rater', 'rated_age']).index
    
    # Process Corrected (df_corr)
    corr_combos = df_corr.set_index(['scale', 'rater', 'rated_age']).index
    df_corr = df_corr[(df_corr['significance'] == True) | (corr_combos.isin(baseline_sig_combos))].copy()
    df_corr = df_corr.drop_duplicates(subset=['scale', 'rater', 'rated_age']).reset_index(drop=True)
    df_corr['err_low'] = df_corr['beta_bc'] - df_corr['low_conf']
    df_corr['err_up'] = df_corr['up_conf'] - df_corr['beta_bc']

    # Define scales and labels manually to control order
    scales = ['emp', 'dsm']
    scale_labels = {'emp': 'AP', 'dsm': 'ADHP'}


    # ==========================================
    # 3. PLOTTING: COMBINED 2x2 GRID
    # ==========================================
    set_style2()
    fig, axes = plt.subplots(2, 2, figsize=(9, 9))
    fig.subplots_adjust(left=0.1, right=0.98, top=0.93, bottom=0.15, wspace=0.2, hspace=0.4)
    
    # NPG Color Palette & Styles
    rater_colors = {"m": "#3C5488", "f": "#F39B7F", "t": "#00A087"}
    age_linestyles = {"7": "-", "10": "--", "12": ":"}
    
    # --- Plot Top Row (Trends) ---
    panels_trend = [
        (axes[0, 0], 'emp', 'AP Scale', 'A.'),
        (axes[0, 1], 'dsm', 'ADHP Scale', 'B.')
    ]
    
    for ax, scale_type, title, panel_letter in panels_trend:
        subset = desc_stats[desc_stats['scale'] == scale_type].copy()
        unique_bins = sorted(subset['yob_bin_plot'].unique())
        
        for (rater, age), grp in subset.groupby(['respondent', 'assigned_age']):
            grp['yob_bin_plot'] = pd.Categorical(grp['yob_bin_plot'], categories=unique_bins, ordered=True)
            grp = grp.sort_values('yob_bin_plot')
            
            ax.plot(grp['yob_bin_plot'].astype(str), grp['mean_response'], 
                    color=rater_colors.get(rater, 'black'), 
                    linestyle=age_linestyles.get(str(int(age)), '-'), 
                    linewidth=1.2, marker='o', markersize=3.5)
            
        ax.set_title(title, fontweight='bold', fontsize=10)
        ax.grid(True, axis='y', color='gray', linestyle='-', linewidth=0.2, alpha=0.5)
        ax.tick_params(axis='x', labelsize = 9)
        ax.set_xlabel("Birth Cohort", size=10, labelpad=10)
        
        ax.set_ylim(0, 4)
        ax.set_yticks([0, 1, 2, 3, 4])
        ax.set_yticklabels([0, 1, 2, 3, 4])
        if panel_letter == 'A.':
            ax.set_ylabel("Mean Sum Score", size=10, labelpad=10)
            
        ax.text(-0.05 if panel_letter == 'B.' else -0.1, 1.05, panel_letter, 
                transform=ax.transAxes, size=12, weight='bold', va='bottom', ha='right')

    # --- Plot Bottom Row (Betas) --- 
    rater_colors = {"m": "#3C5488", "f": "#F39B7F", "t": "#00A087"}
    age_linestyles = {7: '-', 10: '--', 12: ':'}
    
    # Define the panels by scale rather than by adjustment status
    panels_beta = [
        (axes[1, 0], 'emp', 'AP Scale', 'C.'),
        (axes[1, 1], 'dsm', 'ADHP Scale', 'D.')
    ]
    
    # Set the y-axis categories and map them to their respective dataframes
    model_labels = ['Unadjusted', 'Adjusted for\nParental Age']
    model_dfs = [df, df_corr]
    
    for ax, scale_name, title, panel_letter in panels_beta:
        # Alternating background stripes for the two y-axis positions
        for i in range(len(model_labels)):
            if i % 2 == 0:
                ax.axhspan(i - 0.5, i + 0.5, color='gray', alpha=0.05, zorder=0)
                
        for i, model_df in enumerate(model_dfs):
            # Filter the current dataframe for the current panel's scale
            subset = model_df[model_df['scale'] == scale_name]
            n_points = len(subset)
            offsets = np.linspace(-0.3, 0.3, n_points) if n_points > 1 else [0]
            
            for j, (_, row) in enumerate(subset.iterrows()):
                color = rater_colors.get(str(row['rater']).lower(), '#333333')
                age_val = int(row['rated_age'])
                ls = age_linestyles.get(age_val, '-')
                
                is_sig = row.get('significance', True) 
                point_alpha = 1.0 if is_sig else 0.6

                eb = ax.errorbar(row['beta_bc'], i + offsets[j], 
                            xerr=[[row['err_low']], [row['err_up']]],
                            marker='o', linestyle='none', color=color, markersize=3.5, 
                            elinewidth=1.0, capsize=0, alpha=point_alpha)
                
                eb[2][0].set_linestyle(ls)
                
        ax.axvline(0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
        ax.set_title(title, fontweight='bold', fontsize=10)
        ax.set_xlabel("Beta Coefficient (99.17% CI)", labelpad=10)
        ax.tick_params(axis='x', labelsize=9)

        ax.set_xlim(-0.022, 0.002) 
        ax.set_xticks([-0.020, -0.015, -0.010, -0.005, 0.000])
        ax.grid(True, axis='x', color='gray', linestyle='--', linewidth=0.5, alpha=0.3)

        # Invert so 'Unadjusted' (index 0) stays at the top
        ax.invert_yaxis()
        
        # Panel lettering
        ax.text(-0.05 if panel_letter == 'D.' else -0.1, 1.05, panel_letter, 
                transform=ax.transAxes, size=12, weight='bold', va='bottom', ha='right')

    # Apply the unadjusted/adjusted text labels to the left panel
    axes[1, 0].set_yticks(range(len(model_labels)))
    axes[1, 0].set_yticklabels(model_labels, rotation=90, va='center')
    
    # Sync limits to ensure stripes perfectly align across the two bottom plots
    axes[1, 1].set_ylim(axes[1, 0].get_ylim()) 
    axes[1, 0].set_xlim(axes[1, 1].get_xlim())
    
    # Remove borders and numbers between C and D
    axes[1, 1].tick_params(left=False) 
    axes[1, 1].spines['left'].set_visible(False)
    axes[1, 1].set_yticks(axes[1, 0].get_yticks())
    axes[1, 1].set_yticklabels([])
    

    # Unified legends attached to the figure at the bottom
    custom_lines_rater = [Line2D([0], [0], color=c, marker='o', lw=0, markersize=4) for c in rater_colors.values()]
    l1 = fig.legend(custom_lines_rater, ['Mother', 'Father', 'Teacher'], title="Respondent", 
                    frameon=False, bbox_to_anchor=(0.38, 0.08), loc='upper center', ncol=3,
                    fontsize=9, title_fontsize=10, handletextpad=0.4, columnspacing=1.0)
              
    custom_lines_age = [Line2D([0], [0], color='gray', linestyle=ls, lw=1.5) for ls in age_linestyles.values()] 
    l2 = fig.legend(custom_lines_age, ['7', '10', '12'], title="Assessment Age", 
                    frameon=False, bbox_to_anchor=(0.67, 0.08), loc='upper center', ncol=3,
                    fontsize=9, title_fontsize=10, handletextpad=0.4, columnspacing=1.0)

    
    out_path = os.path.join(output_dir, "beta_distributions.pdf")
    plt.savefig(out_path, format='pdf', transparent=True)
    plt.close(fig)


def plot_parental_trends(scale_data_dir, parental_output):
    """
    Reads sumscore_Data.RDS, filters distinct twins prioritizing non-missing 
    parental age and EA, and plots the trends across Year of Birth.
    """
    # ------------------ DATA PROCESSING ------------------
    file_path = os.path.join(scale_data_dir, "sumscore_Data.RDS")
    result = pyreadr.read_r(file_path)
    df = result[None]
    
    # Keep complete cases for base columns
    cols_to_check = ["FID", "PID", "mult_id_fam", "sex", "yob", "zyg", "assigned_age"]
    df = df.dropna(subset=[c for c in cols_to_check if c in df.columns]).copy()
    
    # Filter out teachers
    df = df[df['respondent'] != 't'].copy()


    def process_parent(data, resp_type):
        d = data[data['respondent'] == resp_type].copy()
        
        d['na_age'] = d['parent_Age'].isna()
        d['na_ea'] = d['parent_EA'].isna()
        
        d = d.sort_values(by=['twin_id', 'na_age', 'na_ea'])
        
        d = d.drop_duplicates(subset=['twin_id'], keep='first').copy()
        
        d['parental_age'] = d['parent_Age'] - d['measured_age']

        return d[['participant_id', 'twin_id', 'yob', 'assigned_age', 'respondent', 'parent_EA', 'parental_age']]


    mother_data = process_parent(df, "m")
    father_data = process_parent(df, "f")
    
    # Combine for plotting
    plot_df = pd.concat([mother_data, father_data], ignore_index=True)
    
    # Calculate yearly means, standard deviations, N, and SEM
    trend_stats = plot_df.groupby(['yob', 'respondent']).agg(
        mean_age=('parental_age', 'mean'),
        sd_age=('parental_age', 'std'),
        n_age=('parental_age', 'count'),
        mean_ea=('parent_EA', 'mean'),
        sd_ea=('parent_EA', 'std'),
        n_ea=('parent_EA', 'count')
    ).reset_index()

    trend_stats['sem_age'] = trend_stats['sd_age'] / np.sqrt(trend_stats['n_age'])
    trend_stats['sem_ea'] = trend_stats['sd_ea'] / np.sqrt(trend_stats['n_ea'])


    # ------------------ PLOTTING ------------------
    set_style()
    
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    rater_colors = {"m": "#3C5488", "f": "#F39B7F"}
    
    panels = [
        (axes[0], 'mean_age', 'sem_age', 'Parental Age (Years)', 'A.'),
        (axes[1], 'mean_ea', 'sem_ea', 'Educational Attainment', 'B.')
    ]
    
    for ax, mean_col, sem_col, ylabel, letter in panels:
        for rater in ['m', 'f']:
            subset = trend_stats[trend_stats['respondent'] == rater].sort_values('yob')
            # Drop years where the mean is NaN to prevent broken lines
            subset = subset.dropna(subset=[mean_col])
            
            color = rater_colors[rater]
            
            # Plot line and shaded standard error region
            ax.plot(subset['yob'], subset[mean_col], color=color, linewidth=1.5)
            ax.fill_between(subset['yob'], 
                            subset[mean_col] - subset[sem_col], 
                            subset[mean_col] + subset[sem_col], 
                            color=color, alpha=0.2, linewidth=0)
            
        # Axes formatting
        ax.set_xlabel("Birth Cohorts")
        ax.set_ylabel(ylabel)
        ax.grid(True, axis='y', color='gray', linestyle='-', linewidth=0.2, alpha=0.5)
        ax.tick_params(axis='x', labelsize=7, colors='black')
        ax.tick_params(axis='y', labelsize=7, colors='black')
        
        # Panel lettering
        ax.text(-0.15, 1.05, letter, transform=ax.transAxes, 
                size=12, weight='bold', va='bottom', ha='right')

    # Add single legend at the bottom
    custom_lines = [Line2D([0], [0], color=c, lw=1.5) for c in rater_colors.values()]
    fig.legend(custom_lines, ['Mother', 'Father'], title="Respondent",
               loc='upper center', bbox_to_anchor=(0.5, 0.1), ncol=2, frameon=False)

    # Adjust layout to fit everything and ensure exact size
    fig.subplots_adjust(left=0.1, right=0.95, top=0.90, bottom=0.20, wspace=0.3)
    
    # Save plot
    plt.savefig(parental_output, format='pdf', dpi=600, transparent=True)
    plt.close(fig)

    excel_path = os.path.join(os.path.dirname(parental_output), "parental_trends_descriptives.xlsx")
    trend_stats.to_excel(excel_path, index=False)


def plot_hertiability_trends(parameters_df_path, output_dir):

    df = pd.read_excel(f"{parameters_df_path}/ADE_parameters_corr_SA.xlsx", usecols=["Scale", "Rater", "Age", "Model", "Birth.cohort", "std_varA", "std_varD"])
    df['Age'] = df['Age'].astype(int) 

    # Pre-define models
    selection_rules = [
        # Mother
        {"Scale": "AP", "Rater": "Mother", "Age": 7, "Model": "OM"},
        {"Scale": "AP", "Rater": "Mother", "Age": 10, "Model": "OM"},
        {"Scale": "AP", "Rater": "Mother", "Age": 12, "Model": "ADE"},
        
        {"Scale": "ADHP", "Rater": "Mother", "Age": 7, "Model": "E"},
        {"Scale": "ADHP", "Rater": "Mother", "Age": 10, "Model": "E"},
        {"Scale": "ADHP", "Rater": "Mother", "Age": 12, "Model": "E"},
        
        # Father
        {"Scale": "AP", "Rater": "Father", "Age": 7, "Model": "E"},
        {"Scale": "AP", "Rater": "Father", "Age": 10, "Model": "E"},
        {"Scale": "AP", "Rater": "Father", "Age": 12, "Model": "OM"},
        
        {"Scale": "ADHP", "Rater": "Father", "Age": 7, "Model": "E"},
        {"Scale": "ADHP", "Rater": "Father", "Age": 10, "Model": "OM"},
        {"Scale": "ADHP", "Rater": "Father", "Age": 12, "Model": "OM"},
        
        # Teacher
        {"Scale": "AP", "Rater": "Teacher", "Age": 7, "Model": "OM"},
        {"Scale": "AP", "Rater": "Teacher", "Age": 10, "Model": "OM"},
        {"Scale": "AP", "Rater": "Teacher", "Age": 12, "Model": "OM"},
        
        {"Scale": "ADHP", "Rater": "Teacher", "Age": 7, "Model": "OM"},
        {"Scale": "ADHP", "Rater": "Teacher", "Age": 10, "Model": "OM"},
        {"Scale": "ADHP", "Rater": "Teacher", "Age": 12, "Model": "AD"}
    ]

    # Convert the rules into a DataFrame
    rules_df = pd.DataFrame(selection_rules)

    # Filter
    filtered_df = pd.merge(df, rules_df, on=["Scale", "Rater", "Age", "Model"], how="inner")

    # Calculate Broad Sense Heritability (H^2)
    filtered_df['H2'] = filtered_df['std_varA'] + filtered_df['std_varD']

    # Styling parameters
    rater_colors = {"mother": "#3C5488", "father": "#F39B7F", "teacher": "#00A087", 
                    "m": "#3C5488", "f": "#F39B7F", "t": "#00A087"}
    age_linestyles = {7: '-', 10: '--', 12: ':'}
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={'wspace': 0.15})
    
    # Define the panels by scale
    panels_h2 = [
        (axes[0], 'AP', 'AP Scale', 'A.'),
        (axes[1], 'ADHP', 'ADHP Scale', 'B.')
    ]

    for ax, scale_name, title, panel_letter in panels_h2:
        # Filter dataframe for the current scale
        subset = filtered_df[filtered_df['Scale'].astype(str).str.upper() == scale_name.upper()]
        
        # Group by Rater and Age to plot continuous lines across birth cohorts
        if not subset.empty:
            for (rater, age), group in subset.groupby(['Rater', 'Age']):
                # chronological order
                group = group.sort_values(by='Birth.cohort') 
                
                color = rater_colors.get(str(rater).lower(), '#333333')
                ls = age_linestyles.get(int(age), '-')
                
                # Determine prominence based on the selected Model
                current_model = group['Model'].iloc[0]
                
                if current_model == "OM":
                    line_width = 1.0
                    marker_size=2.5
                    plot_color = "lightgrey"  # Strip color from background models
                    line_zorder = 1           
                else:
                    line_width = 1.5
                    marker_size=3.5
                    plot_color = rater_colors.get(str(rater).lower(), '#333333') # Full color for target models
                    line_zorder = 3           
                
                # Plot the line and points using plot_color
                ax.plot(group['Birth.cohort'], group['H2'], 
                        color=plot_color, linestyle=ls, linewidth=line_width, 
                        marker='o', markersize=marker_size, alpha=0.9, zorder=line_zorder)

        # Formatting
        ax.set_title(title, fontweight='bold', fontsize=10)
        ax.set_xlabel("Birth Cohort", labelpad=10)
        ax.tick_params(axis='both', labelsize=9)
        
        # Clean up borders
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Subtle horizontal grid
        ax.grid(True, axis='y', color='gray', linestyle='--', linewidth=0.5, alpha=0.3)
        ax.grid(False, axis='x')

        # Y-axis labeling and syncing
        if panel_letter == 'A.':
            ax.set_ylabel("Broad Sense Heritability ($H^2$)", labelpad=10)
        else:
            ax.set_ylabel("")
            
        # Panel lettering
        ax.text(-0.15 if panel_letter == 'A.' else -0.05, 1.05, panel_letter, 
                transform=ax.transAxes, size=12, weight='bold', va='bottom', ha='right')


    for ax in axes:
        ax.set_ylim(0.45, 1)
        ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])


    # Unified legends attached to the figure at the bottom
    custom_lines_rater = [Line2D([0], [0], color=c, marker='o', lw=1.5, markersize=4) for c in ["#3C5488", "#F39B7F", 
                                                                                                "#00A087", "lightgrey"]]
    l1 = fig.legend(custom_lines_rater, ['Mother', 'Father', 'Teacher', "No Trend"], title="Respondent", 
                    frameon=False, bbox_to_anchor=(0.35, -0.05), loc='upper center', ncol=4,
                    fontsize=8, title_fontsize=9, handletextpad=0.4, columnspacing=1.0)
              
    custom_lines_age = [Line2D([0], [0], color='gray', linestyle=ls, lw=1.5) for ls in age_linestyles.values()] 
    l2 = fig.legend(custom_lines_age, ['7', '10', '12'], title="Assessment Age", 
                    frameon=False, bbox_to_anchor=(0.65, -0.05), loc='upper center', ncol=3,
                    fontsize=8, title_fontsize=9, handletextpad=0.4, columnspacing=1.0)

    # Adjust layout so legends are not clipped
    plt.tight_layout(rect=[0, 0.08, 1, 1]) 
    
    out_path = os.path.join(output_dir, "heritability_trend.pdf")
    plt.savefig(out_path, format='pdf', transparent=True, bbox_inches='tight')


#dorret
def scale_panel_B(scale_data_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    # ==========================================
    # 1. DATA PROCESSING: TRENDS
    # ==========================================
    file_path = os.path.join(scale_data_dir, "sumscore_Data.RDS")
    result = pyreadr.read_r(file_path)
    df_trends = result[None] 
    
    # Keep complete cases
    cols_to_check = ["FID", "PID", "mult_id_fam", "sex", "yob", "zyg", "assigned_age"]
    df_trends = df_trends.dropna(subset=[c for c in cols_to_check if c in df_trends.columns]).copy()
    
    # Create 6-year bins for birth year
    yob_min = int(df_trends['yob'].min())
    yob_max = int(df_trends['yob'].max())
    breaks = list(range(yob_min, yob_max + 7, 6))
    labels = [f"{breaks[i]}-{breaks[i+1]-1}" for i in range(len(breaks)-1)]
    
    df_trends['yob_bin'] = pd.cut(df_trends['yob'], bins=breaks, labels=labels, right=False, include_lowest=True)
    
    # Bin groups < 200
    df_trends['group_count'] = df_trends.groupby(['scale', 'yob_bin', 'respondent', 'assigned_age'], observed=True)['yob'].transform('count')
    df_trends['yob_bin_plot'] = np.where(df_trends['group_count'] < 200, "2004-2009", df_trends['yob_bin'].astype(str))
    
    # Descriptive stats for plotting
    desc_stats = df_trends.groupby(['scale', 'yob_bin_plot', 'respondent', 'assigned_age'], observed=True).agg(
        mean_response=('score', 'mean'),
        sd_response=('score', 'std'),
        n=('score', 'count')
    ).reset_index()
    desc_stats['sem_response'] = desc_stats['sd_response'] / np.sqrt(desc_stats['n'])


    # ==========================================
    # 2. PLOTTING: SINGLE PANEL (ADHP Scale)
    # ==========================================
    # Assumes set_style2() is defined elsewhere in your environment
    if 'set_style2' in globals():
        set_style2()
    
    # Sized for a single-column layout (approx 89mm/3.5 inches)
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # NPG Color Palette & Styles
    rater_colors = {"m": "#3C5488", "f": "#F39B7F", "t": "#00A087"}
    age_linestyles = {"7": "-", "10": "--", "12": ":"}
    
    # Filter only for the ADHP (dsm) scale
    subset = desc_stats[desc_stats['scale'] == 'dsm'].copy()
    unique_bins = sorted(subset['yob_bin_plot'].unique())
    
    for (rater, age), grp in subset.groupby(['respondent', 'assigned_age']):
        grp['yob_bin_plot'] = pd.Categorical(grp['yob_bin_plot'], categories=unique_bins, ordered=True)
        grp = grp.sort_values('yob_bin_plot')
        
        ax.plot(grp['yob_bin_plot'].astype(str), grp['mean_response'], 
                color=rater_colors.get(rater, 'black'), 
                linestyle=age_linestyles.get(str(int(age)), '-'), 
                linewidth=1.2, marker='o', markersize=3.5)
        
    #ax.set_title('ADHP Scale', fontweight='bold', fontsize=10)
    ax.grid(True, axis='y', color='gray', linestyle='-', linewidth=0.2, alpha=0.5)
    
    # Formatting axes
    ax.tick_params(axis='x', labelsize=8) 
    ax.tick_params(axis='y', labelsize=9)
    ax.set_xlabel("Geboortecohort", size=10, labelpad=8)
    ax.set_ylabel("Gemiddelde Somscore", size=10, labelpad=8)
    ax.set_ylim(0, 4)
    ax.set_yticks([0, 1, 2, 3, 4])
    # Add overarching figure title
    fig.suptitle("Trends in aandachts en hyperactiviteits problemen over 30 jaar", 
                 fontsize=10, fontweight='bold', y=1)
    
    # Despine top and right for a cleaner aesthetic
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Unified legends stacked at the bottom
    custom_lines_rater = [Line2D([0], [0], color=c, marker='o', lw=0, markersize=4) for c in rater_colors.values()]
    fig.legend(custom_lines_rater, ['Moeder', 'Vader', 'Leerkracht'], title="Respondent", 
               frameon=False, bbox_to_anchor=(0.42, 0.12), loc='upper center', ncol=3,
               fontsize=8, title_fontsize=9, handletextpad=0.2, columnspacing=0.8)
              
    custom_lines_age = [Line2D([0], [0], color='gray', linestyle=ls, lw=1.5) for ls in age_linestyles.values()] 
    fig.legend(custom_lines_age, ['7', '10', '12'], title="Leeftijdsclassificatie", 
               frameon=False, bbox_to_anchor=(0.70, 0.12), loc='upper center', ncol=3,
               fontsize=8, title_fontsize=9, handletextpad=0.2, columnspacing=0.8)
    
    # Adjust layout to fit the main panel
    fig.subplots_adjust(left=0.15, right=0.95, top=0.90, bottom=0.25)
    
    out_path = os.path.join(output_dir, "ADHP_trends.pdf")
    # bbox_inches='tight' prevents the exterior legends from being clipped off the page
    plt.savefig(out_path, format='pdf', transparent=True, bbox_inches='tight', dpi=600)


if __name__ == "__main__":
    
    # ------------------ SETUP ------------------
    item_data_dir = "/Users/sezgi/Library/Mobile Documents/com~apple~CloudDocs/ADHD_paper/analyses/item_level_analyses/data_IL"
    trend_res_path ="/Users/sezgi/Library/Mobile Documents/com~apple~CloudDocs/ADHD_paper/analyses/item_level_analyses/model_results_IL/trend_plot_df.xlsx"
    trend_res_path_corr ="/Users/sezgi/Library/Mobile Documents/com~apple~CloudDocs/ADHD_paper/analyses/item_level_analyses/model_results_IL/trend_plot_df_corr.xlsx"
    output_dir = "/Users/sezgi/Library/Mobile Documents/com~apple~CloudDocs/ADHD_paper/analyses/item_level_analyses"
    
    #plot_item_level_trends(item_data_dir, output_dir)
    #plot_beta_distributions(trend_res_path, trend_res_path_corr, output_dir)
    
    scale_trend_path = "/Users/sezgi/Library/Mobile Documents/com~apple~CloudDocs/ADHD_paper/analyses/sum_score_analyses/model_results_SA/openmx_parameters_SA.xlsx"
    scale_trend_path_corr = "/Users/sezgi/Library/Mobile Documents/com~apple~CloudDocs/ADHD_paper/analyses/sum_score_analyses/model_results_SA/corr/openmx_parameters_corr_SA.xlsx"
    scale_data_dir = "/Users/sezgi/Library/Mobile Documents/com~apple~CloudDocs/ADHD_paper/analyses/sum_score_analyses/data_SA"
    scale_out_dir = "/Users/sezgi/Library/Mobile Documents/com~apple~CloudDocs/ADHD_paper/analyses/sum_score_analyses/model_results_SA"
    #scale_combined_figure(scale_data_dir, scale_trend_path, scale_trend_path_corr, scale_out_dir)
    #scale_panel_B(scale_data_dir, scale_out_dir)

    # Parental EA and Age plotting
    parental_output = "/Users/sezgi/Library/Mobile Documents/com~apple~CloudDocs/ADHD_paper/EA_QualityControl/trend_age/parental_trend.pdf"
    #plot_parental_trends(scale_data_dir, parental_output)

    # Heritability plotting
    parameters_dir = "/Users/sezgi/Library/Mobile Documents/com~apple~CloudDocs/ADHD_paper/analyses/sum_score_analyses/ADE_model_Results_SA/corr"
     
    plot_hertiability_trends(parameters_dir, parameters_dir)

    