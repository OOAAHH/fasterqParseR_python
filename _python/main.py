import os
import pandas as pd
import numpy as np
import subprocess
from multiprocessing import Pool, cpu_count

# Define a function to correct names
def correct_names(df):
    if df is None or df.empty:
        raise ValueError("Please supply the output of GetSRAreads")
    elif df['orig_names'].str.contains("fastq.gz").any():
        new_names = df['orig_names'].str.replace(r'SRR.*', '', regex=True)
        df['new_names'] = new_names + df['SRR_ID'].str.replace(r'\_.*', '', regex=True) + df['assigned_read'].str.replace(r'R|I', '_') + ".fastq.gz"
        df['cellranger_names'] = new_names + df['SRR_ID'].str.replace(r'\_.*', '', regex=True) + "_S1_L001_" + df['assigned_read'] + "_001.fastq.gz"
    elif df['orig_names'].str.contains("fastq").any():
        new_names = df['orig_names'].str.replace(r'SRR.*', '', regex=True)
        df['new_names'] = new_names + df['SRR_ID'].str.replace(r'\_.*', '', regex=True) + df['assigned_read'].str.replace(r'R|I', '_') + ".fastq"
        df['cellranger_names'] = new_names + df['SRR_ID'].str.replace(r'\_.*', '', regex=True) + "_S1_L001_" + df['assigned_read'] + "_001.fastq"
    else:
        raise ValueError("Oops... It appears you have no assigned fastq paths")
    return df

# Define the function to assign SRA reads
def assign_sra_reads(working_dir=None, input_dir=None, outdir=None, parallel=False):
    if working_dir is None:
        raise ValueError("Please provide a working directory")
    if input_dir is None:
        raise ValueError("Please provide input directory containing all fasterq-dump outputs")

    os.chdir(working_dir)

    if outdir is None:
        outdir = os.path.join(working_dir, "read_lengths")
        if os.path.exists(outdir):
            subprocess.run(["rm", "-r", outdir])
        os.makedirs(outdir)

    print("First n=250 read lengths will be output into outdir")

    num_cores = cpu_count() - 1

    # Replace with actual whitelist data
    tenXv1, tenXv2, tenXv3 = set(), set(), set()  # Dummy sets for illustration
    whitelists = {"10xv1": tenXv1, "10xv2": tenXv2, "10xv3": tenXv3}

    fasterq_list = [os.path.join(dp, f) for dp, dn, filenames in os.walk(input_dir) for f in filenames if f.endswith((".fastq.gz", ".fastq"))]
    orig_names = [os.path.basename(f).split('.')[0] for f in fasterq_list]
    to_parse = pd.DataFrame({"fasterq_list": fasterq_list, "orig_names": orig_names})
    to_parse['cat'] = np.where(to_parse['fasterq_list'].str.contains('.gz'), 'zcat', 'cat')

    def process_file(x):
        cat_cmd = to_parse.loc[to_parse['fasterq_list'] == x, 'cat'].values[0]
        orig_name = to_parse.loc[to_parse['fasterq_list'] == x, 'orig_names'].values[0]
        read_length_file = os.path.join(outdir, f"{orig_name}.readslength.txt")
        
        subprocess.run(f"{cat_cmd} {x} | head -1000 | awk '{{if(NR%4==2) print length($1)}}' > {read_length_file}", shell=True)
        mean_length = np.mean(np.loadtxt(read_length_file, dtype=int))
        
        assigned_read, chemistry = None, None
        
        if 5 <= mean_length <= 10:
            assigned_read = "I3"
        elif mean_length > 10:
            seqs_file = os.path.join(outdir, f"{orig_name}.seqs.txt")
            subprocess.run(f"{cat_cmd} {x} | head -40000 | awk '{{if(NR%4==2) print /^@/ ? $1 : substr($0,1,16)}}' > {seqs_file}", shell=True)
            seqs = np.loadtxt(seqs_file, dtype=str)
            
            whitelist_counts = {name: np.sum(np.isin(seqs, wl)) for name, wl in whitelists.items()}
            
            if sum(whitelist_counts.values()) > 1000:
                assigned_read = "R1"
                chemistry = max(whitelist_counts, key=whitelist_counts.get)
            else:
                assigned_read = "R2"
        
        return {"SRR_ID": orig_name, "mean_length": mean_length, "assigned_read": assigned_read, "chemistry": chemistry}

    if parallel:
        with Pool(num_cores) as p:
            assigned_sra = p.map(process_file, to_parse['fasterq_list'].tolist())
    else:
        assigned_sra = [process_file(x) for x in to_parse['fasterq_list'].tolist()]

    bound = pd.DataFrame(assigned_sra)
    bound['orig_names'] = fasterq_list

    bound = correct_names(bound)

    bound.to_csv(os.path.join(outdir, "assigned_SRAreads.csv"), index=False)
    return bound

# Define the function to rename all fastas
def rename_all(assigned_sra=None, input_dir=None, format=None):
    if assigned_sra is None or assigned_sra.empty:
        raise ValueError("Please input the output of assignSRAreads")
    if input_dir is None:
        raise ValueError("Please fill input_dir with same input file used for assignSRAreads")
    if format not in ['read_correct', 'cellranger']:
        raise ValueError("Please decide on a renaming format from 'read_correct' or 'cellranger'")

    if format == "read_correct":
        assigned_sra['parent_dir'] = assigned_sra['orig_names'].str.extract(r'(.*/)')
        assigned_sra['dummy_name'] = assigned_sra['parent_dir'] + "dummy" + assigned_sra['new_names'].str.extract(r'(SRR.*)')[0]

        for i, row in assigned_sra.iterrows():
            subprocess.run(f"mv {row['orig_names']} {row['dummy_name']}", shell=True)
        for i, row in assigned_sra.iterrows():
            subprocess.run(f"mv {row['dummy_name']} {row['new_names']}", shell=True)
    elif format == "cellranger":
        for i, row in assigned_sra.iterrows():
            subprocess.run(f"mv {row['orig_names']} {row['cellranger_names']}", shell=True)

    name_check = [os.path.join(dp, f) for dp, dn, filenames in os.walk(input_dir) for f in filenames if f.endswith((".fastq.gz", ".fastq"))]
    assigned_sra['name_check'] = assigned_sra['new_names'].isin(name_check)
    
    outdir = os.path.join(input_dir, "read_lengths")
    assigned_sra.to_csv(os.path.join(outdir, "assigned_SRAreads_final.csv"), index=False)
    print("All done! Please check in assigned_SRAreads.csv for new column name_check and make sure you are happy!")

# Example usage:
# result = assign_sra_reads('/path/to/working_dir', '/path/to/input_dir', '/path/to/outdir', parallel=True)
# rename_all(result, '/path/to/input_dir', format='read_correct')
