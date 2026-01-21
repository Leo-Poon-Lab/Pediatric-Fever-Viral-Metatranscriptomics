import pandas as pd
from sys import argv
# # Read the CSV files
virus_summary = argv[1]
# virus_summary = "virus_summary_with_contigs.tsv"
metagenomic_samples_file = argv[2]
# metagenomic_samples_file = "Patient_Profile_sequencing.tsv"

filtered_output = argv[3]
# filtered_output = "virus_summary_with_contigs_indexhopping.tsv"
removed_output = argv[4]  
# removed_output = "virus_summary_with_contigs_indexhopping_removed.tsv"  # New output file for filtered rows

# Read data
quantification_result = pd.read_csv(virus_summary, sep="\t")
metagenomic_samples = pd.read_csv(metagenomic_samples_file, sep="\t")

# Create sample-group mapping
sample_group_dict = dict(zip(metagenomic_samples['Sample'], metagenomic_samples['Group']))


quantification_result['Group'] = quantification_result['Sample'].map(sample_group_dict)

# --- 新增部分：处理缺失值，生成最终列 ---
# 1. Final_Supporting_Reads: 优先用 Mapped_Reads，缺失时回退到 Number_of_supporting_reads
quantification_result['Final_Supporting_Reads'] = quantification_result['Mapped_Reads'].fillna(
    quantification_result['Number_of_supporting_reads']
)

# 2. Final_RPM: 优先用 RPM_mapping_contigs，缺失时回退到 RPM
quantification_result['Final_RPM'] = quantification_result['RPM_mapping_contigs'].fillna(
    quantification_result['RPM']
)

# --- 更新 index-hopping 过滤逻辑 ---
# 计算每个 Group-Virus_Species 组合的最大 reads
quantification_result['Max_in_Group_Species'] = quantification_result.groupby(
    ['Group', 'Virus_Species']
)['Final_Supporting_Reads'].transform('max')

filter_mask = (
    (quantification_result['Final_Supporting_Reads'] >= 0.001 * quantification_result['Max_in_Group_Species']) | 
    (quantification_result['Final_Supporting_Reads'].isnull()) |
    (quantification_result['Final_Supporting_Reads'] == "")
) 
# (quantification_result['Number_of_supporting_reads'].notnull()) | (quantification_result['Number_of_supporting_reads'] != "")
# Split data
# 分割数据
filtered_result = quantification_result[filter_mask].drop(columns=['Max_in_Group_Species'])
removed_rows = quantification_result[
    ~filter_mask & 
    quantification_result['Final_Supporting_Reads'].notnull() & 
    (quantification_result['Final_Supporting_Reads'] != "")
].drop(columns=['Max_in_Group_Species'])

# Save results
filtered_result.to_csv(filtered_output, index=False, sep="\t")
removed_rows.to_csv(removed_output, index=False, sep="\t")
