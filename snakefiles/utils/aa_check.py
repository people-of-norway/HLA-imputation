# Simple script for checking if bmarker files are consistent with the amino acid sequence strings

import pysam
import pandas as pd

model = "cook"
vcf_file = f"/home/oystein/test/bmarker_cookhla_new.bcf"
aa_file = "/home/oystein/hla_imputation_pipeout/2026.02.24/aa_dict"
alleles_file = "/home/oystein/hla_imputation_pipeout/2026.02.24/merged_alleles"
samples_to_test = 10



aa_df = pd.read_csv(aa_file, sep="\t")
alleles_df = pd.read_csv(alleles_file, sep="\t")

merged_df = alleles_df.merge(aa_df[["hla", "allele", "aa"]].rename(columns={"allele": f"{model}1", "aa": f"aa_{model}1"}), on=["hla", f"{model}1"], how="left")
merged_df = merged_df.merge(aa_df[["hla", "allele", "aa"]].rename(columns={"allele": f"{model}2", "aa": f"aa_{model}2"}), on=["hla", f"{model}2"], how="left")

vcf_in = pysam.VariantFile(vcf_file)

def find_cutoff(hla):
    if hla == "A" or hla == "B" or hla == "C":
        return 24
    elif hla == "DQB1":
        return 32
    elif hla == "DRB1":
        return 29
    elif hla == "DQA1":
        return 23
    elif hla == "DPB1":
        return 29
    
def find_aa(a, s, i):
    j = i + s
    if i < 0:
        return a[j]
    if i == 0:
        return None
    else:
        return a[j-1]
    
samples_tested = 0
for record in vcf_in:
    print(f"Processing record: {record.id}")
    if record.id.startswith(f"AA_"):
        record_list = record.id.split("_")
        test_hla = record_list[1]
        s = find_cutoff(test_hla)
        aa_pos = int(record_list[2])
        samples_tested = 0
        for sample_name in record.samples:
            samples_tested += 1
            sample_row = merged_df[(merged_df["iid"] == sample_name) & (merged_df["hla"] == test_hla)]
            sample_alleles = [str(sample_row[f"{model}1"].values[0]), str(sample_row[f"{model}2"].values[0])]
            sample_aa = [str(sample_row[f"aa_{model}1"].values[0]).replace("z", ""), str(sample_row[f"aa_{model}2"].values[0]).replace("z", "")]
            genotype = record.samples.get(sample_name, {}).get("GT")
            bf_aa = [None, None]
            if genotype is not None:
                if genotype[0] == 1:
                    if record.ref != "a":
                        bf_aa[0] = record.alts[0]
                    else:
                        bf_aa[0] = record.id[-1]
                elif genotype[0] == 0:
                    if record.ref != "a":
                        bf_aa[0] = record.ref

                if genotype[1] == 1:
                    if record.ref != "a":
                        bf_aa[1] = record.alts[0]
                    else:
                        bf_aa[1] = record.id[-1]
                elif genotype[1] == 0:
                    if record.ref != "a":
                        bf_aa[1] = record.ref
            exp_aa = [find_aa(sample_aa[0], s, aa_pos), find_aa(sample_aa[1], s, aa_pos)]
            bf_aa_not_none = sorted([aa for aa in bf_aa if aa is not None])
            exp_aa_sorted = sorted(exp_aa)

            print(f"exp_aa: {exp_aa_sorted}, bf_aa: {bf_aa_not_none}")
            if len(bf_aa_not_none) == 2:
                bf_aa_sorted = sorted(bf_aa)
                if (exp_aa_sorted[0] != "." and bf_aa_sorted[0] != exp_aa_sorted[0]) or (exp_aa_sorted[1] != "." and bf_aa_sorted[1] != exp_aa_sorted[1]):
                    print(f"Discrepancy found for sample {sample_name} at record {record.id}:")
                    print(f"Expected AA: {exp_aa_sorted}, BF AA: {bf_aa_not_none}")
                    print(sample_aa)
            else:
                for aa in bf_aa_not_none:
                    if aa not in exp_aa:
                        print(f"Discrepancy found for sample {sample_name} at record {record.id}:")
                        print(f"Expected AA: {exp_aa_sorted}, BF AA: {bf_aa}")
                        print(sample_aa)
                        break
            if samples_tested >= samples_to_test:
                break