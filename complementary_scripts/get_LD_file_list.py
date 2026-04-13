from pathlib import Path
import pandas as pd
import sys
import s3fs
import os

def download_ukb_ld_chunks(snp_file_path, ld_file_list):
    """
    Automates downloading the exact Price Lab LD chunks needed for a list of SNPs.
    """
    # Connect to the public AWS bucket
    print("Connecting to AWS and fetching the master file list...")
    fs = s3fs.S3FileSystem(anon=True)
    all_files = fs.glob('broad-alkesgroup-ukbb-ld/UKBB_LD/*.npz')
    
    # Build a catalog of available genomic chunks
    chunk_data = []
    for f in all_files:
        # Extract the file name (e.g., chr1_1000000_3000000)
        basename = f.split('/')[-1].replace('.npz', '')
        parts = basename.split('_')
        
        if len(parts) >= 3 and parts[0].startswith('chr'):
            chunk_data.append({
                'file_base': basename,
                'chr': int(parts[0].replace('chr', '')),
                'start': int(parts[1]),
                'end': int(parts[2])
            })
            
    chunk_df = pd.DataFrame(chunk_data)
    
    # Load the variants of interest
    my_snps = pd.read_csv(snp_file_path, compression="gzip", sep='\t',
                          usecols=["variant","rsid","chr","pos",
                                   "dom_sig_total"],
                          dtype={"variant": str,"rsid": str,"chr": int,"pos":int,
                                 "dom_sig_total": int
                                 }
                            )
    # Filter the variants showing significant dominance 
    my_snps = my_snps[my_snps["dom_sig_total"] > 1]
    
    chrom_col = 'chr' 
    pos_col = 'pos'
    
    # Figure out which files you actually need
    print(f"Mapping {len(my_snps)} variants to AWS chunks...")
    files_to_download = set()
    no_data_snps = []
    
    for _, row in my_snps.iterrows():
        # Find chunks that encompass this SNP's position
        match = chunk_df[(chunk_df['chr'] == row[chrom_col]) & 
                         (chunk_df['start'] <= row[pos_col]) & 
                         (chunk_df['end'] >= row[pos_col])]
        
        if not match.empty:
            for _, match_row in match.iterrows():
                files_to_download.add(match_row['file_base'])

        else:
            no_data_snps.append({"chr": row[chrom_col], "pos": row[pos_col]})

    print(f"Found {len(files_to_download)} unique LD chunks required.")
    
    if no_data_snps:
        missing_df = pd.DataFrame(no_data_snps)
        
        summary_df = missing_df.groupby("chr")["pos"].agg(["min", "max"]).reset_index()
        
        print("LD data is not available for the following regions:")
        
        for _, row in summary_df.iterrows():
            print(f" - Chromosome {row['chr']}: between position {row['min']} and {row['max']}")
        
    else:
        print("LD data is available for all required regions.")

    with open(ld_file_list, 'w') as f:
        for base in files_to_download:
            f.write(f"{base}\n")


variants = "/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/significant_SNPs/all_sig_SNPs.tsv.gz"
ld_file_list = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/required_LD_chunks.txt"

download_ukb_ld_chunks(variants, ld_file_list)
