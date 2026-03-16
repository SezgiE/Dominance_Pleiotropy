import pandas as pd


def get_sumstat():



def get_LD_info():


def sumstats_QC():


def get_independent_SNPs():


def main():


def define_gwas_loci(sumstats, ld_pairs, sig_p=(5e-8)/1060, nom_p=0.05, r2_stage1=0.6, r2_stage2=0.1, merge_kb=250):

    """
    Performs two-stage clumping and defines genomic loci from GWAS summary statistics.
    
    Expected Inputs:
    - sumstats: DataFrame with columns ['SNP', 'CHR', 'POS', 'P']
    - ld_pairs: DataFrame with columns ['SNP_A', 'SNP_B', 'r2'] containing all pairs with r2 > 0.1
    """
    
    # 1. Helper function for fast LD lookup
    def get_ld(snp1, snp2):
        if snp1 == snp2: return 1.0
        match = ld_pairs[((ld_pairs['SNP_A'] == snp1) & (ld_pairs['SNP_B'] == snp2)) | 
                         ((ld_pairs['SNP_A'] == snp2) & (ld_pairs['SNP_B'] == snp1))]
        return match['r2'].values[0] if not match.empty else 0.0


    # 2. Stage 1: Define Independent Significant SNPs (r2 < 0.6)
    sig_df = sumstats[sumstats['P'] < sig_p].sort_values('P').copy()
    indep_sig_snps = []
    stage1_clumped = set()

    for _, row in sig_df.iterrows():
        snp = row['SNP']
        if snp in stage1_clumped:
            continue
            
        indep_sig_snps.append(snp)
        
        # Find all other significant SNPs in LD (r2 >= 0.6) to clump them
        for other_snp in sig_df['SNP']:
            if other_snp not in stage1_clumped and get_ld(snp, other_snp) >= r2_stage1:
                stage1_clumped.add(other_snp)


    # 3. Stage 2: Define Lead SNPs (r2 < 0.1)
    # Subset to the independent significant SNPs and apply the stricter threshold
    stage2_df = sumstats[sumstats['SNP'].isin(indep_sig_snps)].sort_values('P')
    lead_snps = []
    stage2_clumped = set()

    for _, row in stage2_df.iterrows():
        snp = row['SNP']
        if snp in stage2_clumped:
            continue
            
        lead_snps.append(snp)
        
        for other_snp in stage2_df['SNP']:
            if other_snp not in stage2_clumped and get_ld(snp, other_snp) >= r2_stage2:
                stage2_clumped.add(other_snp)


    # 4. Define Locus Boundaries (Merging blocks < 250kb)
    loci = []
    for snp in indep_sig_snps:
        snp_chr = sumstats.loc[sumstats['SNP'] == snp, 'CHR'].values[0]
        
        # Find all SNPs with P < 0.05 in LD >= 0.6 with the independent significant SNP
        nom_snps = sumstats[(sumstats['P'] < nom_p) & (sumstats['CHR'] == snp_chr)]
        ld_snps = [s for s in nom_snps['SNP'] if get_ld(snp, s) >= r2_stage1]
        
        ld_snps_pos = sumstats[sumstats['SNP'].isin(ld_snps)]['POS']
        
        if not ld_snps_pos.empty:
            loci.append({
                'indep_snp': snp,
                'chr': snp_chr,
                'start': ld_snps_pos.min(),
                'end': ld_snps_pos.max(),
                'p_val': sumstats.loc[sumstats['SNP'] == snp, 'P'].values[0]
            })

    # Sort loci by genomic position for sequential merging
    loci_df = pd.DataFrame(loci).sort_values(['chr', 'start'])

    # 5. Merge loci within 250kb
    merged_loci = []
    merge_dist_bp = merge_kb * 1000

    if not loci_df.empty:
        current_locus = loci_df.iloc[0].to_dict()
        current_locus['all_indep_snps'] = [current_locus['indep_snp']]

        for i in range(1, len(loci_df)):
            next_locus = loci_df.iloc[i]
            
            # Check if on same chromosome and within merge distance
            if (current_locus['chr'] == next_locus['chr']) and \
               (next_locus['start'] - current_locus['end'] <= merge_dist_bp):
                
                # Expand boundaries and store the newly encompassed independent SNP
                current_locus['end'] = max(current_locus['end'], next_locus['end'])
                current_locus['all_indep_snps'].append(next_locus['indep_snp'])
                # The top SNP represents the locus (minimum P value)
                if next_locus['p_val'] < current_locus['p_val']:
                    current_locus['top_snp'] = next_locus['indep_snp']
                    current_locus['p_val'] = next_locus['p_val']
            else:
                if 'top_snp' not in current_locus:
                    current_locus['top_snp'] = current_locus['indep_snp']
                merged_loci.append(current_locus)
                current_locus = next_locus.to_dict()
                current_locus['all_indep_snps'] = [current_locus['indep_snp']]
        
        if 'top_snp' not in current_locus:
            current_locus['top_snp'] = current_locus['indep_snp']
        merged_loci.append(current_locus)

    final_loci_df = pd.DataFrame(merged_loci)
    
    return {
        'lead_snps': sumstats[sumstats['SNP'].isin(lead_snps)],
        'genomic_loci': final_loci_df[['chr', 'start', 'end', 'top_snp', 'p_val', 'all_indep_snps']]
    }