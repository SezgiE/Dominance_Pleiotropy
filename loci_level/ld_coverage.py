import pandas as pd
import glob
import sys
import os

def check_batch_ld_coverage(gwas_dir, ld_dir, output_summary):
    """
    Scans a directory of GWAS sumstats and calculates LD coverage for each one.
    """
    # ---------------------------------------------------------
    # STEP 1: Build the Master LD Dictionary (Do this ONLY ONCE)
    # ---------------------------------------------------------
    ld_files = glob.glob(os.path.join(ld_dir, "*.gz"))
    ld_positions = set()
    
    print(f"Scanning {len(ld_files)} downloaded LD manifest files...")
    for f in ld_files:
        basename = os.path.basename(f)
        try:
            chrom = int(basename.split('_')[0].replace('chr', ''))
            manifest = pd.read_csv(f, compression='gzip', sep='\t', usecols=['position'])
            
            for pos in manifest['position'].values:
                ld_positions.add((chrom, pos))
        except Exception as e:
            print(f"Could not read {basename}: {e}")
            
    print(f"Master LD Reference built: {len(ld_positions):,} unique variants available.\n")

    # ---------------------------------------------------------
    # STEP 2: Stream the 80+ GWAS Files 
    # ---------------------------------------------------------
    # Adjust the "*.tsv.gz" if your files have a different naming convention
    gwas_files = glob.glob(os.path.join(gwas_dir, "*.tsv.bgz"))
    print(f"Found {len(gwas_files)} GWAS summary statistics files. Starting batch check...")
    
    results = []
    
    for i, gwas_file in enumerate(gwas_files, 1):
        trait_name = os.path.basename(gwas_file).replace('.tsv.gz', '')
        print(f"[{i}/{len(gwas_files)}] Checking coverage for: {trait_name}")
        
        try:
            gwas_df = pd.read_csv(gwas_file, compression='gzip', sep='\t', usecols=['chr', 'pos', "dom_sig"])
            gwas_df = gwas_df[gwas_df["dom_sig"] > 0]
            gwas_positions = set(zip(gwas_df['chr'], gwas_df['pos']))
            total_gwas = len(gwas_positions)
            
            if total_gwas == 0:
                continue
                
            # --- UPDATED LOGIC ---
            # Isolate the exact missing tuples using set difference
            covered_snps = gwas_positions.intersection(ld_positions)
            missing_snps = gwas_positions - ld_positions
            
            coverage_pct = (len(covered_snps) / total_gwas) * 100 if total_gwas > 0 else 0
            
            # Extract unique chromosomes from the missing set, sort them, and join as a string
            missing_chrs = sorted(list(set([chrom for chrom, pos in missing_snps])))
            missing_chrs_str = ", ".join(map(str, missing_chrs)) if missing_chrs else "None"
            
            results.append({
                "Trait": trait_name,
                "Total_SNPs": total_gwas,
                "Covered_SNPs": len(covered_snps),
                "Missing_SNPs": len(missing_snps),
                "Coverage_Pct": round(coverage_pct, 2),
                "Where missings?": missing_chrs_str
            })
            
        except Exception as e:
            print(f"Error processing {trait_name}: {e}")

    # ---------------------------------------------------------
    # STEP 3: Generate the Final QC Report
    # ---------------------------------------------------------
    results_df = pd.DataFrame(results)
    
    # Sort by lowest coverage first so you can easily spot problematic traits
    results_df = results_df.sort_values(by="Coverage_Pct")
    
    results_path = f"{output_summary}/coverage_phen.tsv"
    results_df.to_csv(results_path, sep='\t', index=False)
    print("\n" + "="*50)
    print(f"Batch coverage check complete! Summary saved to: {results_path}")
    print("="*50)


def check_pleiotropic_coverage(all_sig_path, ld_dir, output_summary):
    """
    Filters for variants with dom_sig_total > 1 and calculates LD coverage 
    broken down by chromosome.
    """
    print("Loading and filtering significant SNPs...")
    
    # 1. Load the data and filter for dom_sig_total > 1
    # Translating your R fread logic to pandas
    df = pd.read_csv(
        all_sig_path, 
        compression='gzip', 
        sep='\t', 
        usecols=['variant', 'chr', 'pos', 'dom_sig_total']
    )
    
    pleio_df = df[df['dom_sig_total'] > 1].copy()
    
    if pleio_df.empty:
        print("No variants found with dom_sig_total > 1.")
        return pd.DataFrame()
        
    print(f"Found {len(pleio_df):,} pleiotropic variants across {pleio_df['chr'].nunique()} chromosomes.")

    # 2. Identify exactly which chromosomes we need to check
    required_chrs = sorted(pleio_df['chr'].unique())
    ld_positions = set()
    
    print("Building LD reference dictionary for required chromosomes...")
    
    # 3. Load LD manifests ONLY for the required chromosomes
    for chrom in required_chrs:
        # We only search for files matching this specific chromosome
        ld_files = glob.glob(os.path.join(ld_dir, f"chr{chrom}_*.gz"))
        
        for f in ld_files:
            try:
                manifest = pd.read_csv(f, compression='gzip', sep='\t', usecols=['position'])
                for pos in manifest['position'].values:
                    ld_positions.add((chrom, pos))
            except Exception as e:
                print(f"Could not read {os.path.basename(f)}: {e}")

    # 4. Calculate coverage grouped by chromosome
    results = []
    print("Calculating per-chromosome coverage...")
    
    for chrom in required_chrs:
        chrom_df = pleio_df[pleio_df['chr'] == chrom]
        total_snps = len(chrom_df)
        
        # Create a set of the pleiotropic positions for this specific chromosome
        gwas_positions = set(zip(chrom_df['chr'], chrom_df['pos']))
        
        # Calculate overlap
        covered_snps = gwas_positions.intersection(ld_positions)
        missing_snps = gwas_positions - ld_positions
        coverage_pct = (len(covered_snps) / total_snps) * 100 if total_snps > 0 else 0
        
        results.append({
            "Chromosome": chrom,
            "Pleiotropic_SNPs": total_snps,
            "Covered_SNPs": len(covered_snps),
            "Missing_SNPs": len(missing_snps),
            "Coverage_Pct": round(coverage_pct, 2)
        })

    # 5. Format the final output
    results_df = pd.DataFrame(results)
    
    # Add a "Total" row at the very bottom for a quick summary
    total_pleio = results_df["Pleiotropic_SNPs"].sum()
    total_covered = results_df["Covered_SNPs"].sum()
    total_missing = results_df["Missing_SNPs"].sum()
    total_pct = (total_covered / total_pleio) * 100 if total_pleio > 0 else 0
    
    total_row = pd.DataFrame([{
        "Chromosome": "ALL",
        "Pleiotropic_SNPs": total_pleio,
        "Covered_SNPs": total_covered,
        "Missing_SNPs": total_missing,
        "Coverage_Pct": round(total_pct, 2)
    }])
    
    results_df = pd.concat([results_df, total_row], ignore_index=True)

    results_path = f"{output_summary}/coverage_pleio.tsv"
    results_df.to_csv(results_path, sep='\t', index=False)
    
    print("\n" + "="*60)
    print("PLEIOTROPIC SNP COVERAGE REPORT (dom_sig_total > 1)")
    print("="*60)
    print(results_df.to_string(index=False))
    print("="*60)


# --- RUN THE SCRIPT ---
# Update these to your actual SLURM directories

if __name__ == "__main__":
    
    # Check if a task ID was passed from the terminal
    if len(sys.argv) < 2:
        print("Error: missing arguments")
        sys.exit(1)

    gwas_directory = sys.argv[1]
    ld_directory = sys.argv[2]
    output_summary = sys.argv[3]
    all_sig_path = sys.argv[4]
    
    check_batch_ld_coverage(gwas_directory, ld_directory, output_summary)
    check_pleiotropic_coverage(all_sig_path, ld_directory, output_summary)