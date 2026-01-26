import argparse
import pandas as pd
import subprocess
import os
import sys
import re
import shutil
import glob
import gzip

# Terminal Styling
class Style:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    RESET = ENDC

def print_header(msg):
    print(f"{Style.BOLD}{Style.HEADER}{msg}{Style.ENDC}")

def print_info(msg):
    print(f"{Style.BLUE}[INFO] {msg}{Style.ENDC}")

def print_success(msg):
    print(f"{Style.GREEN}[SUCCESS] {msg}{Style.ENDC}")

def print_warning(msg):
    print(f"{Style.WARNING}[WARNING] {msg}{Style.ENDC}")

def print_error(msg):
    print(f"{Style.FAIL}[ERROR] {msg}{Style.ENDC}")

def check_dependencies(tools):
    missing = []
    for tool in tools:
        if shutil.which(tool) is None:
            missing.append(tool)
    if missing:
        print_error(f"Missing required tools: {', '.join(missing)}")
        sys.exit(1)

def run_cmd(cmd, verbose=False):
    if verbose:
        print(f"[CMD] {cmd}")
    try:
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    except subprocess.CalledProcessError:
        print_error(f"Command failed: {cmd}")
        sys.exit(1)


def parse_ranges(range_str):
    parts = range_str.split(',')
    valid_regions = []
    
    # Valid Chromosomes: 1-22, X, Y, MT, M
    valid_chrs = set([str(i) for i in range(1, 23)] + ['X', 'Y', 'M', 'MT'])
    
    for part in parts:
        part = part.strip()
        if not part: continue
        
        # Regex to match chr:start-end
        match = re.match(r'^((?:chr)?([0-9A-Za-z]+)):(\d+)-(\d+)$', part)
        if match:
            full_chr = match.group(1)
            raw_chr = match.group(2).upper().replace('CHR', '')
            start = match.group(3)
            end = match.group(4)
            
            if raw_chr in valid_chrs:
                # Store tuple for formatted output
                valid_regions.append((full_chr, start, end))
    
    return valid_regions

def get_gene_range_from_refflat(refflat_file, gene):
    # Retrieve union of all transcripts for the gene
    # refFlat: geneName name chrom strand txStart txEnd ...
    if not os.path.exists(refflat_file):
        return None
    
    ranges = []
    try:
        with gzip.open(refflat_file, 'rt') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6: continue
                if parts[0] == gene:
                    chrom = parts[2]
                    start = parts[4]
                    end = parts[5]
                    # RVTest often uses the format chr:start-end
                    # Ensure chrom has 'chr' prefix if needed? 
                    # RefFlat usually doesn't have 'chr' in some versions, check input
                    # Our refFlat seems to not have 'chr' based on filename 'nochr'
                    # But VCF usually has non-chr or chr depending on build.
                    # The parts[2] from user `cat` was '3' (no chr)
                    # We will format as matches the parse_ranges expectation
                    ranges.append(f"{chrom}:{start}-{end}")
                    
        if not ranges:
            return None
        return ",".join(ranges)
        
    except Exception as e:
        print(f"[ERROR] Reading refFlat: {e}")
        return None

def get_rvtest_info(assoc_file, gene):
    # Columns: Gene RANGE N_INFORMATIVE NumVar NumPolyVar Q rho Pvalue
    try:
        # Check delimiter, usually tab or whitespace
        df = pd.read_csv(assoc_file, sep=r'\s+', engine='python')
        row = df[df['Gene'] == gene]
        if row.empty:
            return None, None
        
        # RANGE is comma separated list of regions
        gene_range = row.iloc[0]['RANGE']
        num_var = row.iloc[0]['NumVar']
        return gene_range, num_var
    except Exception as e:
        print(f"[ERROR] Reading assoc file: {e}")
        sys.exit(1)

def get_extended_stats(file_path, gene, fields):
    if not file_path or not os.path.exists(file_path):
        return {f: "NA" for f in fields}
    try:
        df = pd.read_csv(file_path, sep=r'\s+', engine='python')
        row = df[df['Gene'] == gene]
        if row.empty:
            return {f: "NA" for f in fields}
        res = {}
        for f in fields:
            if f in row.columns:
                res[f] = str(row.iloc[0][f])
            else:
                res[f] = "NA"
        return res
    except Exception as e:
        print(f"[WARNING] Could not read stats from {file_path}: {e}")
        return {f: "NA" for f in fields}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", required=True)
    parser.add_argument("--assoc-file", required=True)
    parser.add_argument("--burden-file", required=False, help="Path to Burden FDR file")
    parser.add_argument("--skato-file", required=False, help="Path to SKAT-O FDR file")
    parser.add_argument("--vcf-file", required=True)
    parser.add_argument("--plink-prefix", required=True)
    parser.add_argument("--tommo-vcf", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--out-log", required=True)
    parser.add_argument("--group-name", required=False, default="", help="Analysis group suffix for file naming")
    parser.add_argument("--plink2-path", required=True)
    parser.add_argument("--pheno-file", required=True)
    parser.add_argument("--refflat-file", required=False)
    
    # Sample Detail Options
    parser.add_argument("--sample-group", default="both", choices=["case", "control", "both"], help="Sample group to analyze for details (default: both)")
    parser.add_argument("--min-variant-count", type=int, default=1, help="Minimum variants to trigger alert (default: 1)")
    parser.add_argument("--no-sample-details", action='store_true', help="Disable sample level detail calculation")

    args = parser.parse_args()
    
    # 0. Robustness Checks
    check_dependencies(['bcftools'])
    # plink2 path is passed as arg, check existence
    if not os.path.exists(args.plink2_path):
        print_error(f"Plink2 executable not found at: {args.plink2_path}")
        sys.exit(1)
    
    # 1. Get Gene Info
    gene_range_str = None
    rvtest_numvar = "NA"
    
    # Try getting info from Refflat first if provided
    if args.refflat_file:
        print_info(f"Reading RefFlat: {args.refflat_file}")
        gene_range_str = get_gene_range_from_refflat(args.refflat_file, args.gene)
        
    # Get NumVar from Assoc file (and Range if refFlat failed or not provided)
    assoc_range, num_var = get_rvtest_info(args.assoc_file, args.gene)
    if num_var is not None:
         rvtest_numvar = num_var
         
    if not gene_range_str:
        if assoc_range:
            print_info("Using range from Assoc File.")
            gene_range_str = assoc_range
        else:
             print_error(f"Gene {args.gene} not found in assoc file and refFlat invalid.")
             sys.exit(1)
    else:
         print_info("Using range from RefFlat.")

    # 2. Parse Ranges
    valid_ranges = parse_ranges(gene_range_str)
    if not valid_ranges:
        print_error("No valid ranges found for gene.")
        sys.exit(1)
    
    # 3. Extract IDs from VCF using regions
    print_info("Extracting variants from VCF...")
    
    # Construct region string for -r (format: chr:start-end,chr:start-end)
    region_strs = []
    for (c, s, e) in valid_ranges:
        region_strs.append(f"{c}:{s}-{e}")
    region_arg = ",".join(region_strs)
    
    vcf_ids_file = os.path.join(args.out_dir, "vcf_extracted_ids.txt")
    vcf_info_file = os.path.join(args.out_dir, "vcf_extracted_info.txt")
    
    # Run bcftools query to get IDs
    cmd_query = (f"bcftools query -f '%ID\\n' "
                 f"-r {region_arg} "
                 f"{args.vcf_file} > {vcf_ids_file}")
                 
    run_cmd(cmd_query, verbose=True)

    # Load VCF Info into Dictionary
    vcf_info_dict = {}

    # Helper function to try extraction of ONE tag from a list of candidates
    def extract_field_robust(candidates, out_str_key, idx_in_dict):
        # candidates: list of strings like ["%INFO/IMPACT", "%INFO/impact"]
        # idx_in_dict: 0 for Impact, 1 for Effect
        
        found = False
        temp_out = os.path.join(args.out_dir, f"vcf_extract_{out_str_key}.txt")
        
        for cand in candidates:
            # Try query
            cmd = (f"bcftools query -f '%ID\\t{cand}\\n' "
                   f"-r {region_arg} "
                   f"{args.vcf_file} > {temp_out}")
            
            # print(f"[DEBUG] Trying {cand}")
            try:
                subprocess.check_call(cmd, shell=True, executable='/bin/bash', stderr=subprocess.DEVNULL)
                # If we are here, it worked.
                found = True
                break
            except subprocess.CalledProcessError:
                continue
        
        if found and os.path.exists(temp_out):
            with open(temp_out, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if not parts: continue

                    vid = parts[0]
                    if not vid: continue
                    
                    # Logic: If value is missing, empty, or ".", treat as "."
                    val = "."
                    if len(parts) >= 2:
                        raw_val = parts[1].strip()
                        if raw_val and raw_val != ".":
                            val = raw_val
                        
                    if vid not in vcf_info_dict: vcf_info_dict[vid] = [".", "."]
                    vcf_info_dict[vid][idx_in_dict] = val
        else:
             print_warning(f"Could not extract {out_str_key} (tried: {candidates})")

    # 1. Try IMPACT (upper or lower)
    extract_field_robust(["%INFO/IMPACT", "%INFO/impact"], "impact", 0)

    # 2. Try EFFECT (upper or lower)
    extract_field_robust(["%INFO/EFFECT", "%INFO/effect"], "effect", 1)

    # 4. Prepare Case/Control Lists from Pheno File
    print_info("Preparing sample lists...")
    try:
        # Expected cols: fid, iid, ..., pheno1
        # Try tab first
        pheno_df = pd.read_csv(args.pheno_file, sep='\t')
        if 'pheno1' not in pheno_df.columns:
             # Try whitespace
             pheno_df = pd.read_csv(args.pheno_file, sep=r'\s+', engine='python')
        
        if 'pheno1' not in pheno_df.columns:
            print_error("Could not find 'pheno1' in pheno file.")
            sys.exit(1)
            
        # [Strict] Filter against FAM file to ensure sample existence
        fam_file = args.plink_prefix + ".fam"
        if os.path.exists(fam_file):
            # FAM format: FID IID FAT MOT SEX PHENO
            fam_df = pd.read_csv(fam_file, sep=r'\s+', header=None, usecols=[0, 1], names=['fid', 'iid'], dtype=str)
            valid_samples = set(zip(fam_df['fid'], fam_df['iid']))
            
            # Ensure pheno cols are str
            pheno_df['fid'] = pheno_df['fid'].astype(str)
            pheno_df['iid'] = pheno_df['iid'].astype(str)
            
            # Filter
            mask = pheno_df.apply(lambda r: (str(r['fid']), str(r['iid'])) in valid_samples, axis=1)
            pheno_df = pheno_df[mask]
            
            print_info(f"Pheno file filtered by FAM. Active Samples: {len(pheno_df)}")
        else:
             print_warning(f"FAM file {fam_file} not found. Using raw pheno list.")
            
        # Case = 2, Control = 1
        case_file = os.path.join(args.out_dir, "cases.txt")
        ctrl_file = os.path.join(args.out_dir, "controls.txt")
        
        # Write FID IID for PLINK --keep
        pheno_df[pheno_df['pheno1'] == 2][['fid', 'iid']].to_csv(case_file, sep='\t', index=False, header=False)
        pheno_df[pheno_df['pheno1'] == 1][['fid', 'iid']].to_csv(ctrl_file, sep='\t', index=False, header=False)
        
        n_case = len(pheno_df[pheno_df['pheno1'] == 2])
        n_ctrl = len(pheno_df[pheno_df['pheno1'] == 1])
        print_info(f"Analysis Groups -> Cases: {n_case}, Controls: {n_ctrl}")
        
    except Exception as e:
        print_error(f"Processing pheno file: {e}")
        sys.exit(1)
    
    # Run Plink for Cases
    print_info("Running PLINK for Cases...")
    prefix_case = os.path.join(args.out_dir, "stats_case")
    # Output: .afreq (AAF), .gcount (Genocounts), .vmiss (Missing)
    # Removed 'mac' from cols, PLINK2 calculates freq/obs_ct
    cmd_mk_case = (f"{args.plink2_path} --bfile {args.plink_prefix} "
                   f"--extract {vcf_ids_file} "
                   f"--keep {case_file} "
                   f"--freq cols=chrom,pos,ref,alt,alt1,nobs,altfreq "
                   f"--geno-counts "
                   f"--missing "
                   f"--out {prefix_case} --threads 4 > /dev/null")
    run_cmd(cmd_mk_case)
    
    # Run Plink for Ctrls
    print_info("Running PLINK for Controls...")
    prefix_ctrl = os.path.join(args.out_dir, "stats_ctrl")
    cmd_mk_ctrl = (f"{args.plink2_path} --bfile {args.plink_prefix} "
                   f"--extract {vcf_ids_file} "
                   f"--keep {ctrl_file} "
                   f"--freq cols=chrom,pos,ref,alt,alt1,nobs,altfreq "
                   f"--geno-counts "
                   f"--missing "
                   f"--out {prefix_ctrl} --threads 4 > /dev/null")
    run_cmd(cmd_mk_ctrl)
    
    # 5. Extract Tommo Frequencies
    print_info("Extracting Tommo frequencies...")
    tommo_out = os.path.join(args.out_dir, "tommo_af.txt")
    tommo_dict = {}
    
    if args.tommo_vcf and os.path.exists(args.tommo_vcf):
        # Build Regions from extracted IDs (Found in VCF)
        # ToMMo requires 'chr' prefix usually, whereas our VCF might use '3'
        # Our IDs are typically formatted as chr:pos:ref:alt (from pipeline)
        
        tommo_regions_file = os.path.join(args.out_dir, "tommo_regions.txt")
        unique_regions = set()
        
        if os.path.exists(vcf_ids_file):
            with open(vcf_ids_file, 'r') as f:
                for line in f:
                    vid = line.strip()
                    if not vid: continue
                    # Parse ID: chr:pos:ref:alt
                    parts = vid.split(':')
                    if len(parts) >= 2:
                        chrom = parts[0]
                        pos = parts[1]
                        
                        # Normalize for ToMMo (Needs chr prefix usually)
                        if not chrom.startswith('chr'):
                            chrom = 'chr' + chrom
                        
                        unique_regions.add(f"{chrom}\t{pos}")
        
        if unique_regions:
            with open(tommo_regions_file, 'w') as f:
                for reg in unique_regions:
                    f.write(reg + "\n")
            
            # Tommo VCF query using specific regions
            # Added %ID to get rsID
            cmd_tommo = (f"bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF\\t%ID\\n' "
                         f"-R {tommo_regions_file} "
                         f"{args.tommo_vcf} > {tommo_out}")
            run_cmd(cmd_tommo, verbose=True)
        else:
             print("[WARNING] No variant IDs found to query Tommo.")

        # Read Tommo into dict
        if os.path.exists(tommo_out):
            try:
                # CHROM POS REF ALT AF ID
                # Check if file is empty
                if os.stat(tommo_out).st_size > 0:
                    t_df = pd.read_csv(tommo_out, sep='\t', header=None, names=['CHROM', 'POS', 'REF', 'ALT', 'AF', 'ID'])
                    for _, row in t_df.iterrows():
                        # Normalize CHROM key to strip 'chr' to match Plink output '3'
                        c_key = str(row['CHROM']).replace('chr', '')
                        key = f"{c_key}:{row['POS']}:{row['REF']}:{row['ALT']}"
                        
                        # Handle rsID
                        rsid_val = str(row['ID'])
                        if not rsid_val or rsid_val.lower() == 'nan': rsid_val = "."
                        
                        tommo_dict[key] = (row['AF'], rsid_val)
            except Exception as e:
                print_warning(f"Failed to read Tommo output: {e}")
    else:
        print_warning("Tommo VCF not found or not provided.")

    # 6. Aggregate Data and Write Log
    
    try:
        # Load Frequency / MAC
        # Plink2 .afreq format: #CHROM POS ID REF ALT ALT1 OBS_CT ALT_FREQS (or similar)
        case_freq = pd.read_csv(prefix_case + ".afreq", sep='\t')
        ctrl_freq = pd.read_csv(prefix_ctrl + ".afreq", sep='\t')
        
        # Load Geno Counts
        # Plink2 .gcount format: #CHROM POS ID REF ALT HOM_REF HET HOM_ALT MISSING
        case_gcount = pd.read_csv(prefix_case + ".gcount", sep='\t')
        ctrl_gcount = pd.read_csv(prefix_ctrl + ".gcount", sep='\t')
        
        # Load Missing
        # Plink2 .vmiss format: #CHROM POS ID REF ALT F_MISS
        case_miss = pd.read_csv(prefix_case + ".vmiss", sep='\t')
        ctrl_miss = pd.read_csv(prefix_ctrl + ".vmiss", sep='\t')

        # Clean Column Names (Strip whitespace)
        for df in [case_freq, ctrl_freq, case_gcount, ctrl_gcount, case_miss, ctrl_miss]:
            df.columns = df.columns.str.strip()
            
    except Exception as e:
        print_error(f"Failed to load PLINK stats: {e}")
        # If files empty (no variants found in plink), handle gracefully
        case_freq = pd.DataFrame()
        ctrl_freq = pd.DataFrame()
        case_gcount = pd.DataFrame()
        ctrl_gcount = pd.DataFrame()
        case_miss = pd.DataFrame()
        ctrl_miss = pd.DataFrame()
    
    # Identify Freq Column Name (ALT_FREQS or ALT_FREQ)
    freq_col_name = 'ALT_FREQS'
    if not case_freq.empty:
        found = next((c for c in case_freq.columns if 'FREQ' in c), None)
        if found: freq_col_name = found

    # Prepare List of Records
    variant_records = []
    target_allele_dict = {}
    
    total_mac_case = 0
    total_mac_ctrl = 0
    total_alleles_case = 0
    total_alleles_ctrl = 0
    
    # Helper to get value
    def get_val(df, id_val, col):
        if df.empty: return 0
        if 'ID' not in df.columns: return 0
        row = df[df['ID'] == id_val]
        if row.empty: return 0
        try:
             return row.iloc[0][col]
        except:
             return 0

    # Get set of all IDs from freq files (variants present in Plink)
    all_plink_ids = set()
    if not case_freq.empty and 'ID' in case_freq.columns:
        all_plink_ids = all_plink_ids.union(set(case_freq['ID']))
    if not ctrl_freq.empty and 'ID' in ctrl_freq.columns:
        all_plink_ids = all_plink_ids.union(set(ctrl_freq['ID']))
    
    for vid in all_plink_ids:
        # Basic Info
        row_c = pd.DataFrame()
        if not case_freq.empty and 'ID' in case_freq.columns:
             row_c = case_freq[case_freq['ID'] == vid]
        
        if row_c.empty and not ctrl_freq.empty and 'ID' in ctrl_freq.columns:
            row_c = ctrl_freq[ctrl_freq['ID'] == vid]
            
        if row_c.empty: continue

        chrom = str(row_c.iloc[0]['#CHROM'])
        pos = str(row_c.iloc[0]['POS'])
        ref = str(row_c.iloc[0]['REF'])
        alt_col = next((c for c in row_c.columns if c.startswith('ALT')), 'ALT')
        alt = str(row_c.iloc[0][alt_col]) 
        
        # MAC Calculation
        obs_case = get_val(case_freq, vid, 'OBS_CT')
        obs_ctrl = get_val(ctrl_freq, vid, 'OBS_CT')
        
        af_case = get_val(case_freq, vid, freq_col_name)
        af_ctrl = get_val(ctrl_freq, vid, freq_col_name)
        
        mac_case = int(round(obs_case * af_case))
        mac_ctrl = int(round(obs_ctrl * af_ctrl))

        # Guarantee Minor Allele Calculation
        # Check if ALT is major (Pooled AF > 0.5)
        comb_mac = mac_case + mac_ctrl
        comb_obs = obs_case + obs_ctrl
        
        is_alt_major = False
        if comb_obs > 0:
            if (comb_mac / comb_obs) > 0.5:
                is_alt_major = True
        
        if is_alt_major:
            # Flip to REF counting
            row_mac_case = obs_case - mac_case
            row_mac_ctrl = obs_ctrl - mac_ctrl
            target_allele_dict[vid] = ref
            
            # MAF is Ref Freq (1 - AltFreq)
            maf_case_val = 1.0 - af_case
            maf_ctrl_val = 1.0 - af_ctrl
            minor_allele = ref
        else:
            row_mac_case = mac_case
            row_mac_ctrl = mac_ctrl
            target_allele_dict[vid] = alt
            
            # MAF is Alt Freq
            maf_case_val = af_case
            maf_ctrl_val = af_ctrl
            minor_allele = alt
        
        total_mac_case += row_mac_case
        total_mac_ctrl += row_mac_ctrl
        total_alleles_case += obs_case
        total_alleles_ctrl += obs_ctrl
        
        # Geno Counts
        # PLINK2 default cols: HOM_REF_CT HET_REF_ALT_CTS TWO_ALT_GENO_CTS MISSING_CT
        g_homref_case = get_val(case_gcount, vid, 'HOM_REF_CT')
        g_het_case = get_val(case_gcount, vid, 'HET_REF_ALT_CTS')
        g_homalt_case = get_val(case_gcount, vid, 'TWO_ALT_GENO_CTS')
        g_miss_case = get_val(case_gcount, vid, 'MISSING_CT')

        g_homref_ctrl = get_val(ctrl_gcount, vid, 'HOM_REF_CT')
        g_het_ctrl = get_val(ctrl_gcount, vid, 'HET_REF_ALT_CTS')
        g_homalt_ctrl = get_val(ctrl_gcount, vid, 'TWO_ALT_GENO_CTS')
        g_miss_ctrl = get_val(ctrl_gcount, vid, 'MISSING_CT')
        
        # Miss Rate
        miss_rate_case = get_val(case_miss, vid, 'F_MISS')
        miss_rate_ctrl = get_val(ctrl_miss, vid, 'F_MISS')
        
        # Tommo
        tommo_key = f"{chrom}:{pos}:{ref}:{alt}"
        val = tommo_dict.get(tommo_key, ("NA", "."))
        
        val_af = "NA"
        val_rsid = "."
        
        if isinstance(val, tuple):
             val_af = val[0]
             val_rsid = val[1]
        else:
             val_af = val
             
        if val_af != "NA":
             tommo_af = f"{float(val_af):.6f}"
        else:
             tommo_af = "NA"
             
        # VCF Info
        (imp, eff) = vcf_info_dict.get(vid, (".", "."))
        
        rec = {
            "SNPID": vid,
            "rsID": val_rsid,
            "Is_Alt_Minor": "Yes" if not is_alt_major else "No",
            "MinorAllele": minor_allele,
            "Impact": imp,
            "Effect": eff,
            "Case_Geno": f"{int(g_homref_case)}/{int(g_het_case)}/{int(g_homalt_case)}/{int(g_miss_case)}",
            "Ctrl_Geno": f"{int(g_homref_ctrl)}/{int(g_het_ctrl)}/{int(g_homalt_ctrl)}/{int(g_miss_ctrl)}",
            "Case_MAF": f"{maf_case_val:.6f}",
            "Ctrl_MAF": f"{maf_ctrl_val:.6f}",
            "Case_AAF": f"{af_case:.6f}",
            "Ctrl_AAF": f"{af_ctrl:.6f}",
            "Case_MissRate": f"{miss_rate_case:.6f}",
            "Ctrl_MissRate": f"{miss_rate_ctrl:.6f}",
            "ToMMo_AAF": tommo_af
        }
        variant_records.append(rec)


    # Sort variant_records
    def variant_sort_key(record):
        # Format: chr:pos:ref:alt (e.g., chr17:80343161:G:A or 17:80343161:G:A)
        # Note: Some VCFs use 'chr1', others '1'. Logic should handle both.
        snpid = record["SNPID"]
        parts = snpid.split(':')
        if len(parts) < 2:
            return (99, 0) # Fallback
            
        chrom_str = parts[0].replace('chr', '')
        pos_str = parts[1]
        
        # Chromosome Order
        if chrom_str.isdigit():
            c_val = int(chrom_str)
        elif chrom_str == 'X':
            c_val = 23
        elif chrom_str == 'Y':
            c_val = 24
        elif chrom_str == 'M' or chrom_str == 'MT':
            c_val = 25
        else:
            c_val = 99
            
        # Position
        try:
            p_val = int(pos_str)
        except:
            p_val = 0
            
        return (c_val, p_val)

    variant_records.sort(key=variant_sort_key)
    
    # 7. Write Log File
    num_vcf_hits = len(variant_records)
    with open(args.out_log, 'w') as f:
        # Get Additional Stats
        burden_stats = get_extended_stats(args.burden_file, args.gene, ['Pvalue', 'FDR'])
        skato_stats = get_extended_stats(args.skato_file, args.gene, ['Pvalue', 'FDR', 'rho'])

        # Part 1: Gene Stats
        f.write(f"=== Gene Summary: {args.gene} ===\n")
        f.write(f"Assoc_File: {args.assoc_file}\n")
        f.write(f"VCF_File: {args.vcf_file}\n")
        f.write(f"Plink_Prefix: {args.plink_prefix}\n")
        f.write(f"Burden_File: {args.burden_file if args.burden_file else 'NA'}\n")
        f.write(f"SKAT-O_File: {args.skato_file if args.skato_file else 'NA'}\n")
        f.write(f"RVTest_NumVar: {rvtest_numvar}\n")
        f.write(f"VCF_Hit_NumVar: {num_vcf_hits}\n")
        f.write(f"Burden_Pvalue: {burden_stats['Pvalue']}\n")
        f.write(f"Burden_FDR: {burden_stats['FDR']}\n")
        f.write(f"SKAT-O_Pvalue: {skato_stats['Pvalue']}\n")
        f.write(f"SKAT-O_FDR: {skato_stats['FDR']}\n")
        f.write(f"SKAT-O_Rho: {skato_stats['rho']}\n")
        f.write(f"Total_MAC_Gene: {int(total_mac_case + total_mac_ctrl)}\n")
        
        # Burden Frequency Calculation (Moved up for Terminal Display)
        total_grp_alleles_case = 2 * n_case
        total_grp_alleles_ctrl = 2 * n_ctrl
        ratio_case = total_mac_case / total_grp_alleles_case if total_grp_alleles_case > 0 else 0
        ratio_ctrl = total_mac_ctrl / total_grp_alleles_ctrl if total_grp_alleles_ctrl > 0 else 0

        # Check if Is_Alt_Minor is always "Yes" (Needed for Terminal & Log)
        all_alt_minor = all(r['Is_Alt_Minor'] == "Yes" for r in variant_records)

        # Terminal Output for Gene Stats (Stylish & Complete)
        print(f"\n{Style.HEADER}┌── Gene Analysis: {Style.BOLD}{args.gene}{Style.RESET} {Style.HEADER}─────────────────{Style.RESET}")
        print(f"{Style.HEADER}│{Style.RESET}  RVTest NumVar  : {rvtest_numvar}")
        print(f"{Style.HEADER}│{Style.RESET}  VCF Hit NumVar : {num_vcf_hits}")
        print(f"{Style.HEADER}│{Style.RESET}  Burden P-value : {burden_stats['Pvalue']} (FDR: {burden_stats['FDR']})")
        print(f"{Style.HEADER}│{Style.RESET}  SKAT-O P-value : {skato_stats['Pvalue']} (FDR: {skato_stats['FDR']})")
        print(f"{Style.HEADER}│{Style.RESET}  SKAT-O Rho     : {skato_stats['rho']}")
        print(f"{Style.HEADER}│{Style.RESET}  Total MAC      : {int(total_mac_case + total_mac_ctrl)} (Case: {int(total_mac_case)}, Ctrl: {int(total_mac_ctrl)})")
        print(f"{Style.HEADER}│{Style.RESET}  Ratio Case     : {ratio_case:.6f}")
        print(f"{Style.HEADER}│{Style.RESET}  Ratio Ctrl     : {ratio_ctrl:.6f}")
        print(f"{Style.HEADER}└────────────────────────────────────────{Style.RESET}")

        # Terminal Output for Variant Details (Restored & Professionalized & Full Info)
        if not args.no_sample_details: # Using this flag to control verbosity generally, or just print it.
            print(f"\n{Style.BOLD}Detailed Variant List ({num_vcf_hits}):{Style.RESET}")
            
            # 1. Define Headers (Same as Log)
            header_cols = ["SNPID", "rsID", "Is_Alt_Minor", "Impact", "Effect", "Case_Geno(Ref/Het/Alt/Miss)", "Ctrl_Geno(Ref/Het/Alt/Miss)"]
            if not all_alt_minor:
                 header_cols.extend(["Case_MAF", "Ctrl_MAF"])
            header_cols.extend(["Case_AAF", "Ctrl_AAF", "ToMMo_AAF", "Case_MissRate", "Ctrl_MissRate"])
            
            table_data = [header_cols]
            
            for rec in variant_records:
                # Colorize Impact
                imp = rec["Impact"]
                if "HIGH" in imp: imp = f"{Style.FAIL}{imp}{Style.RESET}"
                elif "MODERATE" in imp: imp = f"{Style.WARNING}{imp}{Style.RESET}"
                
                # Colors for MAF (if exists)
                c_maf = rec.get("Case_MAF", "NA")
                if c_maf != "NA":
                    try:
                        if float(c_maf) > 0.01: c_maf = f"{Style.FAIL}{c_maf}{Style.RESET}"
                    except: pass

                # Genotypes: Full Format with Color
                def fmt_geno_full(g_str):
                    parts = g_str.split('/')
                    if len(parts) >= 3:
                        # Logic to handle varying parts length if Miss is present or not
                        # Assuming Ref/Het/Alt/Miss from log logic
                        # But wait, rec["Case_Geno"] holds the string we put in log logic.
                        # Let's check log logic loop above. 
                        # Log logic just writes rec["Case_Geno"].
                        # We need to see how rec["Case_Geno"] was constructed.
                        # It was constructed lines 573-590 approx.
                        # Let's assume standard logic: Ref/Het/Alt/Miss
                        r_val, h_val, a_val = parts[0], parts[1], parts[2]
                        
                        # Highlighting
                        if int(h_val) > 0: h_val = f"{Style.WARNING}{h_val}{Style.RESET}"
                        if int(a_val) > 0: a_val = f"{Style.FAIL}{a_val}{Style.RESET}"
                        
                        ret = f"{r_val}/{h_val}/{a_val}"
                        if len(parts) > 3:
                            ret += f"/{parts[3]}"
                        return ret
                    return g_str

                c_geno = fmt_geno_full(rec["Case_Geno"])
                n_geno = fmt_geno_full(rec["Ctrl_Geno"])
                
                # Build Row
                row = [
                    rec["SNPID"],
                    rec["rsID"],
                    rec["Is_Alt_Minor"],
                    imp,
                    rec["Effect"],
                    c_geno,
                    n_geno
                ]
                
                if not all_alt_minor:
                    row.extend([c_maf, rec["Ctrl_MAF"]])
                
                row.extend([
                    rec["Case_AAF"],
                    rec["Ctrl_AAF"],
                    rec["ToMMo_AAF"],
                    rec["Case_MissRate"],
                    rec["Ctrl_MissRate"]
                ])
                
                table_data.append(row)
            
            # 2. Calculate Widths (handling ANSI)
            # We need regex to strip ansi
            import re
            def get_visible_len(s):
                ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
                return len(ansi_escape.sub('', str(s)))

            col_widths = [0] * len(header_cols)
            for row in table_data:
                for i, val in enumerate(row):
                    w = get_visible_len(val)
                    if w > col_widths[i]: col_widths[i] = w
            
            # 3. Print Table
            for i, row in enumerate(table_data):
                line_parts = []
                for j, val in enumerate(row):
                    # Left Align
                    v_len = get_visible_len(val)
                    padding = " " * (col_widths[j] - v_len + 2) # +2 padding
                    line_parts.append(str(val) + padding)
                
                line_str = "".join(line_parts)
                if i == 0:
                    print(f"  {Style.UNDERLINE}{line_str}{Style.RESET}")
                else:
                    print(f"  {line_str}")



        f.write(f"Total_MAC_Case: {int(total_mac_case)}\n")
        f.write(f"Total_MAC_Ctrl: {int(total_mac_ctrl)}\n")
        f.write(f"Ratio_MAC_Case (MAC/(2*N)): {ratio_case:.6f}\n")
        f.write(f"Ratio_MAC_Ctrl (MAC/(2*N)): {ratio_ctrl:.6f}\n")
        f.write("\n")

        
        # Part 2: Variant Details
        # Check if Is_Alt_Minor is always "Yes"
        all_alt_minor = all(r['Is_Alt_Minor'] == "Yes" for r in variant_records)

        f.write("=== Variant Details ===\n")
        header = ["SNPID", "rsID", "Is_Alt_Minor", "Impact", "Effect", "Case_Geno(Ref/Het/Alt/Miss)", "Ctrl_Geno(Ref/Het/Alt/Miss)"]
        
        if not all_alt_minor:
             header.extend(["Case_MAF", "Ctrl_MAF"])
             
        header.extend(["Case_AAF", "Ctrl_AAF", "ToMMo_AAF", "Case_MissRate", "Ctrl_MissRate"])
        f.write("\t".join(header) + "\n")
        
        for rec in variant_records:
            row = [
                rec["SNPID"],
                rec["rsID"],
                rec["Is_Alt_Minor"],
                rec["Impact"],
                rec["Effect"],
                rec["Case_Geno"],
                rec["Ctrl_Geno"]
            ]
            
            if not all_alt_minor:
                row.extend([rec["Case_MAF"], rec["Ctrl_MAF"]])
            
            row.extend([
                rec["Case_AAF"],
                rec["Ctrl_AAF"],
                rec["ToMMo_AAF"],
                rec["Case_MissRate"],
                rec["Ctrl_MissRate"]
            ])
            
            f.write("\t".join(row) + "\n")

    # 8. Sample Details Module
    if not args.no_sample_details:
        print(f"\n{Style.BOLD}================ Sample Level Analysis ================{Style.RESET}")
        
        # Determine target file and count
        targets = []
        if args.sample_group == 'both':
            targets.append(('case', case_file, n_case))
            targets.append(('control', ctrl_file, n_ctrl))
        elif args.sample_group == 'case':
            targets.append(('case', case_file, n_case))
        elif args.sample_group == 'control':
            targets.append(('control', ctrl_file, n_ctrl))
        
        for (g_name, g_file, g_n) in targets:
            # Output prefix
            sample_out_prefix = os.path.join(args.out_dir, f"sample_stats_{g_name}")
            
            # Run Plink Export A
            try:
                cmd_export = (f"{args.plink2_path} --bfile {args.plink_prefix} "
                              f"--extract {vcf_ids_file} "
                              f"--keep {g_file} "
                              f"--export A "
                              f"--out {sample_out_prefix} --threads 4 > /dev/null")
                subprocess.check_call(cmd_export, shell=True, executable='/bin/bash', stderr=subprocess.DEVNULL)
            except Exception:
                 print_warning(f"Sample export failed for {g_name} (possibly no variants). Skipping.")
                 continue
            
            raw_file = sample_out_prefix + ".raw"
            if os.path.exists(raw_file):
                try:
                    df_raw = pd.read_csv(raw_file, sep='\t')
                    
                    # Columns: FID IID PAT MAT SEX PHENO ... SNPs...
                    if len(df_raw.columns) > 6:
                        snp_cols = df_raw.columns[6:]
                        
                        results = []
                        
                        for _, row in df_raw.iterrows():
                            fid = str(row['FID'])
                            iid = str(row['IID'])
                            
                            carrier_list = []
                            total_mac = 0
                            
                            for col in snp_cols:
                                val = row[col]
                                if pd.isna(val): continue
                                
                                # Parse Variant ID from column name
                                # Plink2 raw headers are typically ID_ALLELE
                                split_parts = col.rsplit('_', 1)
                                if len(split_parts) < 2: 
                                    # Fallback if no underscore (unlikely for --export A)
                                    var_id_clean = col
                                    counted_allele = None
                                else:
                                    var_id_clean = split_parts[0]
                                    counted_allele = split_parts[1]
                                
                                # Determine effective count of Minor Allele
                                target = target_allele_dict.get(var_id_clean)
                                final_val = 0
                                
                                # Logic to handle Major/Minor Inversion
                                if target and counted_allele:
                                    if counted_allele == target:
                                        # Counting the minor allele
                                        final_val = val
                                    else:
                                        # Counting the major allele -> Invert (0->2, 1->1, 2->0)
                                        final_val = 2 - val 
                                else:
                                    # Fallback: Assume the output counts the ALT/Minor allele
                                    final_val = val

                                # Determine Label (Standard 0/0, 0/1, 1/1) based on ALT Count (val)
                                gt_label = "./."
                                if not pd.isna(val):
                                    ival = int(val)
                                    if ival == 0: gt_label = "0/0"
                                    elif ival == 1: gt_label = "0/1"
                                    elif ival == 2: gt_label = "1/1"
                                    else: gt_label = f"?({ival})"

                                if final_val > 0:
                                    carrier_list.append(f"{var_id_clean}({gt_label})")
                                    total_mac += final_val
                            
                            results.append({
                                "FID": fid,
                                "IID": iid,
                                "Variant_Count": len(carrier_list),
                                "Total_MAC": int(total_mac),
                                "Variants": ";".join(carrier_list)
                            })
                        
                        df_res = pd.DataFrame(results)
                        
                        # Sort: Variant_Count Descending, Total_MAC Descending
                        df_res = df_res.sort_values(by=["Variant_Count", "Total_MAC"], ascending=[False, False])

                        # Save
                        if args.group_name:
                            fn_middle = f"{args.group_name}.{g_name}"
                        else:
                            fn_middle = g_name
                            
                        out_details = os.path.join(os.path.dirname(args.out_log), f"{args.gene}.{fn_middle}.sample_details.tsv")
                        df_res.to_csv(out_details, sep='\t', index=False)
                        print_success(f"Sample details saved to: {out_details}")
                        
                        # Terminal Output
                        df_high = df_res[df_res['Variant_Count'] >= args.min_variant_count]
                        
                        n_selected = len(df_high)
                        pct = (n_selected / g_n * 100) if g_n > 0 else 0.0
                        
                        # Stylish Summary
                        print(f"\n{Style.HEADER}┌── Summary: {g_name.capitalize()} {Style.RESET}")
                        print(f"{Style.HEADER}│{Style.RESET}  {'Total Samples':<25}: {Style.BOLD}{g_n}{Style.RESET}")
                        pct_color = Style.GREEN if pct > 0 else Style.WARNING
                        print(f"{Style.HEADER}│{Style.RESET}  {'Samples (>= ' + str(args.min_variant_count) + ' vars)':<25}: {pct_color}{n_selected} ({pct:.2f}%){Style.RESET}")
                        print(f"{Style.HEADER}└────────────────────────────────{Style.RESET}")

                        if not df_high.empty:
                            print(f"\n  {Style.BOLD}Detailed Sample List:{Style.RESET}")
                            
                            # Manual Printing for Professional Coloring
                            h_fid = "FID"
                            h_iid = "IID"
                            h_nv = "N_Vars"
                            h_mac = "MAC"
                            h_var = "Variants"
                            
                            # Header
                            print(f"  {Style.UNDERLINE}{h_fid:<15} {h_iid:<15} {h_nv:<8} {h_mac:<6} {h_var}{Style.RESET}")
                            
                            for _, row in df_high.iterrows():
                                fid = str(row['FID'])
                                iid = str(row['IID'])
                                n_vars = str(row['Variant_Count'])
                                mac = str(row['Total_MAC'])
                                variants_str = str(row['Variants'])
                                
                                # Colorize Genotypes in Variants String
                                # Format: rsID(GT) -> 0/1 (Yellow/Warning), 1/1 (Red/Fail)
                                colored_vars = []
                                if variants_str and variants_str != "nan":
                                    raw_vars = variants_str.split(';')
                                    for v in raw_vars:
                                        if '(' in v and ')' in v:
                                            try:
                                                # Split rs123(0/1)
                                                base, gt_part = v.split('(')
                                                gt = gt_part.rstrip(')')
                                                
                                                if gt == "1/1":
                                                    gt_colored = f"{Style.FAIL}{gt}{Style.RESET}"
                                                elif gt == "0/1":
                                                    gt_colored = f"{Style.WARNING}{gt}{Style.RESET}"
                                                else:
                                                    gt_colored = gt
                                                
                                                colored_vars.append(f"{base}({gt_colored})")
                                            except:
                                                colored_vars.append(v)
                                        else:
                                            colored_vars.append(v)
                                    
                                    final_var_str = "; ".join(colored_vars)
                                else:
                                    final_var_str = ""
                                
                                print(f"  {fid:<15} {iid:<15} {n_vars:<8} {mac:<6} {final_var_str}")
                            
                            print(f"{Style.CYAN}{'-' * 80}{Style.RESET}")
                        else:
                            print(f"  {Style.WARNING}No samples met the criteria.{Style.RESET}")

                    else:
                        print_warning(f"No variants found in exported sample data for {g_name}.")
                        
                except Exception as e:
                    print_error(f"Processing sample details for {g_name}: {e}")

if __name__ == "__main__":
    main()
