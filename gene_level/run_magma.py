import subprocess
import bisect

def run_magma(magma_directory):
    command = [
        "./magma",
        "--annotate", "window=1",
        "--snp-loc", "snp_loc38.txt",
        "--gene-loc", "genes_loc38.loc",
        "--out", "magma_pleio_mapping"
    ]

    try:
        print("Starting MAGMA annotation...")
        result = subprocess.run(command, 
                                cwd=magma_directory,
                                check=True, 
                                capture_output=True, 
                                text=True)
        
        print("Process finished successfully!")
    except subprocess.CalledProcessError as e:
        print(f"Error: MAGMA process failed with exit code {e.returncode}")
        print(e.stderr)
    except FileNotFoundError:
        print("Error: The './magma' executable could not be found.")


def get_unmapped_snps(magma_directory):
    print("\nExtracting unmapped SNPs...")
    try:
        all_snps = set()
        with open(f"{magma_directory}/snp_loc38.txt", 'r') as f:
            for line in f:
                parts = line.strip().split()
                if parts: 
                    all_snps.add(parts[0])

        mapped_snps = set()
        with open(f"{magma_directory}/magma_pleio_mapping.genes.annot", 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split()
                if len(parts) > 2:
                    mapped_snps.update(parts[2:])

        unmapped_snps = all_snps - mapped_snps
        print(f"Found {len(unmapped_snps)} unmapped SNPs.")
        
        return unmapped_snps

    except FileNotFoundError as e:
        print(f"File Error during comparison: {e}")
        return set()


def map_nearest_gene(magma_directory, unmapped_set):
    if not unmapped_set:
        print("No unmapped SNPs to process.")
        return

    print("\nMapping unmapped SNPs and formatting as MAGMA .annot...")
    
    genes_by_chr = {}
    genes_info = {} 
    with open(f"{magma_directory}/genes_loc38.loc", 'r') as f:
        
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4:
                gene, chrom = parts[0], parts[1]
                try:
                    start, stop = int(parts[2]), int(parts[3])
                    if chrom not in genes_by_chr:
                        genes_by_chr[chrom] = []
                    genes_by_chr[chrom].append((start, stop, gene))
                    
                    genes_info[gene] = (chrom, start, stop)
                except ValueError:
                    continue


    for chrom in genes_by_chr:
        genes_by_chr[chrom].sort(key=lambda x: x[0])


    snp_locations = {}
    with open(f"{magma_directory}/snp_loc38.txt", 'r') as f:
        
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                snp = parts[0]
                if snp in unmapped_set:
                    chrom = parts[1]
                    try:
                        bp = int(parts[2])
                        snp_locations[snp] = (chrom, bp)
                    except ValueError:
                        continue
    

    rescued_mapping = {}
    distance_log = [] 
    for snp, (chrom, bp) in snp_locations.items():
       
        if chrom not in genes_by_chr:
            continue
        
        chrom_genes = genes_by_chr[chrom]
        starts = [g[0] for g in chrom_genes] 
        
        idx = bisect.bisect_left(starts, bp)
        
        min_dist = float('inf')
        closest_gene = None
        
        candidates = []
        if idx < len(chrom_genes): candidates.append(chrom_genes[idx])
        if idx > 0: candidates.append(chrom_genes[idx-1])
            
        for start, stop, gene in candidates:
            if start <= bp <= stop:
                dist = 0 
            else:
                dist = min(abs(bp - start), abs(bp - stop))
            
            if dist < min_dist:
                min_dist = dist
                closest_gene = gene
                
        if closest_gene:
            if closest_gene not in rescued_mapping:
                rescued_mapping[closest_gene] = []
            rescued_mapping[closest_gene].append(snp)
            distance_log.append((snp, chrom, bp, closest_gene, min_dist))
    
    output_file = f"{magma_directory}/manuel_pleio_mapping.genes.annot"
    
    with open(output_file, 'w') as f:
        # Add the exact same headers MAGMA outputs
        f.write("# window_up = 1000\n")
        f.write("# window_down = 1000\n")
        
        for gene, snps in rescued_mapping.items():
            chrom, start, stop = genes_info[gene]
            snp_string = "\t".join(snps)
            
            # Formatted like  MAGMA: GENE \t CHR:START:STOP \t SNP1 \t SNP2
            f.write(f"{gene}\t{chrom}:{start}:{stop}\t{snp_string}\n")

    
    log_file = f"{magma_directory}/manuel_pleio_mapping.genes_log.txt"
    with open(log_file, 'w') as f:
        # Write a clean header for easy importing into R/Pandas
        f.write("SNP\tCHR\tBP\tMAPPED_GENE\tMIN_DISTANCE\n")
        # Sort the log by Chromosome and then BP for easier reading
        distance_log.sort(key=lambda x: (int(x[1]) if x[1].isdigit() else x[1], x[2]))
        for item in distance_log:
            f.write(f"{item[0]}\t{item[1]}\t{item[2]}\t{item[3]}\t{item[4]}\n")
            
    total_rescued_snps = sum(len(snps) for snps in rescued_mapping.values())
    print(f"Successfully mapped {total_rescued_snps} SNPs to {len(rescued_mapping)} genes.")
    print(f"Results saved: {output_file}")


if __name__ == "__main__":

    magma_directory = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/magma/magma_v1"

    run_magma(magma_directory)
    unmapped_set = get_unmapped_snps(magma_directory)
    map_nearest_gene(magma_directory,unmapped_set)