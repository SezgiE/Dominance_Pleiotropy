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


def plot_beta_distributions(trend_res_path, output_dir):

    df = pd.read_excel(trend_res_path, usecols=["q_name","rater","rated_age","significance", 
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

    age_markers = {7: 'o', 10: '^', 12: 's'}  # Circle, Triangle, Square
    

    for i, item in enumerate(items):
        subset = df[df['q_name'] == item]
        n_points = len(subset)
        

        offsets = np.linspace(-0.3, 0.3, n_points) if n_points > 1 else [0]
        
        for j, (_, row) in enumerate(subset.iterrows()):
            color = rater_colors.get(str(row['rater']).lower(), '#333333')
            

            age_val = int(row['rated_age'])
            marker_shape = age_markers.get(age_val, 'o')
            

            ax.errorbar(row['beta_bc'], i + offsets[j], 
                        xerr=[[row['err_low']], [row['err_up']]], 
                        marker=marker_shape, linestyle='none', color=color, markersize=3.5, 
                        elinewidth=1.0, capsize=0)
            

    ax.axvline(0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    

    ax.set_yticks(range(len(items)))
    ax.set_yticklabels([f"Item {item.replace('q', '')}" for item in items])
    ax.set_xlabel("Beta Coefficient (95% CI)", labelpad=20)
    

    custom_lines_rater = [Line2D([0], [0], color=c, marker='o', lw=0, markersize=4) for c in rater_colors.values()]
    l1 = ax.legend(custom_lines_rater, ['Mother', 'Father', 'Teacher'], title="Respondent", 
              frameon=False, bbox_to_anchor=(1.02, 1), loc='upper left')
              
    custom_lines_age = [Line2D([0], [0], color='gray', marker=m, lw=0, markersize=4) for m in age_markers.values()]
    l2 = ax.legend(custom_lines_age, ['7', '10', '12'], title="Assessment Age", 
              frameon=False, bbox_to_anchor=(1.02, 0.65), loc='upper left')
              
    ax.add_artist(l1) 


    fig.subplots_adjust(left=0.12, right=0.80, top=0.92, bottom=0.15)

    
    out_path = os.path.join(f"{output_dir}/model_results_IL", "beta_distributions.pdf")
    plt.savefig(out_path, format='pdf', transparent=True)
    plt.close(fig)
    print(f"Beta distribution plot saved successfully to {out_path}")


def scale_combined_figure(trend_res_path, scale_trend_path, output_dir):
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
    df_betas = pd.read_excel(scale_trend_path, usecols=["scale","rater","rated_age","significance", 
                                                "beta_bc", "low_conf", "up_conf", "lower_thresh", "upper_thresh"])
    df_betas = df_betas[df_betas['significance'] == True].copy()
    df_betas = df_betas.drop_duplicates(subset=['scale', 'rater', 'rated_age']).reset_index(drop=True)

    ap_scale_df = df_betas[df_betas['scale'] == 'emp'].copy()
    adhp_scale_df = df_betas[df_betas['scale'] == 'dsm'].copy()
    
    # Calculate error margins
    for tmp_df in [ap_scale_df, adhp_scale_df]:
        tmp_df['err_low'] = tmp_df['beta_bc'] - tmp_df['low_conf']
        tmp_df['err_up'] = tmp_df['up_conf'] - tmp_df['beta_bc']
        tmp_df.sort_values(by=['rated_age', 'rater'], inplace=True)


    # ==========================================
    # 3. PLOTTING: COMBINED 2x2 GRID
    # ==========================================
    set_style2()
    
    # Initialize 2x2 grid
    fig, axes = plt.subplots(2, 2, figsize=(8, 6))
    
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
            
        ax.set_title(title, fontweight='bold', fontsize=9)
        ax.grid(True, axis='y', color='gray', linestyle='-', linewidth=0.2, alpha=0.5)
        ax.tick_params(axis='x', labelsize=8)
        ax.set_xlabel("Birth Cohort", size=10)
        
        ax.set_ylim(0, 4)
        ax.set_yticks([0, 1, 2, 3, 4])
        ax.set_yticklabels([0, 1, 2, 3, 4], fontsize=8)
        if panel_letter == 'A.':
            ax.set_ylabel("Mean Sum Score")
            
        ax.text(-0.15, 1.05, panel_letter, transform=ax.transAxes, 
                size=10, weight='bold', va='bottom', ha='right')

    # --- Plot Bottom Row (Betas) ---
    panels_beta = [
        (axes[1, 0], ap_scale_df, 'AP Scale', 'C.'),
        (axes[1, 1], adhp_scale_df, 'ADHP Scale', 'D.')
    ]
    
    for ax, panel_df, title, panel_letter in panels_beta:
        scales = panel_df['scale'].unique()
                
        for i, scale_name in enumerate(scales):
            subset = panel_df[panel_df['scale'] == scale_name]
            n_points = len(subset)
            offsets = np.linspace(-0.3, 0.3, n_points) if n_points > 1 else [0]
            
            for j, (_, row) in enumerate(subset.iterrows()):
                color = rater_colors.get(str(row['rater']).lower(), '#333333')
                age_val = int(row['rated_age'])
                
                ls = age_linestyles.get(str(age_val), '-')
                
                eb = ax.errorbar(row['beta_bc'], i + offsets[j], 
                            xerr=[[row['err_low']], [row['err_up']]], 
                            marker='o', linestyle='none', color=color, markersize=3.5, 
                            elinewidth=1.0, capsize=0)
                
                eb[2][0].set_linestyle(ls)
                
                
        ax.axvline(0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
        ax.set_yticks(range(len(scales)))
        ax.set_yticklabels([])
        
        if panel_letter == 'C.':
            ax.set_ylabel(title, fontsize=10)
        else:
            ax.set_ylabel(title, fontsize=10)
            
        ax.set_xlabel("Beta Coefficient (95% CI)", size=10,labelpad=10)
        
        ax.set_xlim(-0.021, 0.001) 
        ax.set_xticks([-0.020, -0.015, -0.010, -0.005, 0.000])
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(left=False, bottom=True)

        ax.text(-0.15, 1.05, panel_letter, transform=ax.transAxes, 
                size=10, weight='bold', va='bottom', ha='right')

    # --- Unified Legends ---
    custom_lines_rater = [Line2D([0], [0], color=c, lw=1.5, marker='o', markersize=3.5) for c in rater_colors.values()]
    l1 = fig.legend(custom_lines_rater, ['Mother', 'Father', 'Teacher'], title="Respondent",
                     loc='upper center', bbox_to_anchor=(0.35, 0.1), ncol=3, frameon=False,
                     fontsize=8, title_fontsize=10, handletextpad=0.4, columnspacing=1.0)
                     
    custom_lines_age = [Line2D([0], [0], color='black', lw=1.5, ls=ls) for ls in age_linestyles.values()]
    l2 = fig.legend(custom_lines_age, ['7', '10', '12'], title="Age",
                     loc='upper center', bbox_to_anchor=(0.65, 0.1), ncol=3, frameon=False,
                     fontsize=8, title_fontsize=10, handletextpad=0.4, columnspacing=1.0)

    # --- Layout & Save ---
    fig.subplots_adjust(left=0.1, right=0.95, top=0.92, bottom=0.2, wspace=0.25, hspace=0.4)

    out_path = os.path.join(output_dir, "SA_combined.pdf")
    plt.savefig(out_path, format='pdf', transparent=True)
    plt.close(fig)


if __name__ == "__main__":
    
    # ------------------ SETUP ------------------
    item_data_dir = "/Users/sezgi/Library/Mobile Documents/com~apple~CloudDocs/ADHD_paper/analyses/item_level_analyses/data_IL"
    trend_res_path ="/Users/sezgi/Library/Mobile Documents/com~apple~CloudDocs/ADHD_paper/analyses/item_level_analyses/model_results_IL/trend_plot_df.xlsx"
    output_dir = "/Users/sezgi/Library/Mobile Documents/com~apple~CloudDocs/ADHD_paper/analyses/item_level_analyses"
    
    #plot_item_level_trends(item_data_dir, output_dir)
    #plot_beta_distributions(trend_res_path, output_dir)
    
    scale_trend_path = "/Users/sezgi/Library/Mobile Documents/com~apple~CloudDocs/ADHD_paper/analyses/sum_score_analyses/model_results_SA/openmx_parameters_SA.xlsx"
    scale_data_dir = "/Users/sezgi/Library/Mobile Documents/com~apple~CloudDocs/ADHD_paper/analyses/sum_score_analyses/data_SA"
    scale_out_dir = "/Users/sezgi/Library/Mobile Documents/com~apple~CloudDocs/ADHD_paper/analyses/sum_score_analyses/descriptives_SA"
    scale_combined_figure(scale_data_dir, scale_trend_path, scale_out_dir)


    