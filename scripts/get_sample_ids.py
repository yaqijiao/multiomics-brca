# get_sample_ids.py

"""
Multi-omics Sample Intersection Extraction Script
Purpose:
1. Read sample IDs from import data files;
2. Find the intersection of sample IDs between the two omics datasets
3. Save the intersected sample IDs to a file for downstream analysis

Usage:
python get_sample_ids.py

"""

import os
import pandas as pd
 
def read_sample_ids(data_paths):
    # data_paths = {"expr": expr_path, "meth": meth_path, "mirna": mir_path}
    dfs = {}
    sample_sets = {}

    for omic_type, path in data_paths.items():
        df = pd.read_csv(path, sep="\t", index_col=0)
        dfs[omic_type] = df
        sample_ids = set("-".join(col.split("-")[:3]) for col in df.columns)
        sample_sets[omic_type] = sample_ids
        print(f"{omic_type}: {len(sample_ids)} samples")

    common_samples = set.intersection(*sample_sets.values())   
    return dfs, common_samples


def subset_dfs_by_sample_ids(dfs, common_samples):
    subset_dfs = {}

    for omic_type, df in dfs.items():
        matched_cols = [col for col in df.columns if "-".join(col.split("-")[:3]) in common_samples]
        subset_dfs[omic_type] = df[matched_cols]
        print(f"{omic_type}: {len(matched_cols)} samples after subsetting")

    return subset_dfs

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
    
if __name__ == "__main__":
    os.chdir("E:\VSCfiles\multiomics-brca")

    output_dir = "results"
    output_file = os.path.join(output_dir, "common_sample_ids.txt")
    subset_output_dir= os.path.join(output_dir, "subsetted")
    

    data_paths = {
        "expr": r"data\expression\TCGA.BRCA.sampleMap_HiSeqV2\HiSeqV2",
        "meth": r"data\methylation\TCGA.BRCA.sampleMap_HumanMethylation450\HumanMethylation450",
        "mirna": r"data\mirna\TCGA.BRCA.sampleMap_miRNA_HiSeq_gene\miRNA_HiSeq_gene"
    }

    dfs, common_samples = read_sample_ids(data_paths)
    ensure_dir(output_dir)

    with open(output_file, "w") as f:
        for sid in sorted(common_samples):
            f.write(sid + "\n")
    print(f"Common sample IDs written to {output_file}")

    subsetted_dfs = subset_dfs_by_sample_ids(dfs, common_samples)
    ensure_dir(subset_output_dir)
    
    for omic_type, df in subsetted_dfs.items():
        out_path = os.path.join(subset_output_dir, f"{omic_type}_subset.tsv")
        df.to_csv(out_path, sep="\t")
        print(f"{omic_type} subset saved to {out_path}")