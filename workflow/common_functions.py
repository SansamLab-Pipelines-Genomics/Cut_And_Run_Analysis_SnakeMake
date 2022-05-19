# this imports the pandas package functionality in an object named pd
import pandas as pd
import re

# this reads the CSV file and sets an index using the values in the "sample" column.
samples_table = pd.read_csv(config["samples_csv"]).set_index("sample", drop=False)
samples_table = samples_table.applymap(str)

# fastq filename input function definition set to Python dictionary
def fq_dict_from_sample(wildcards):
  return {
    "fq1": samples_table.loc[wildcards.sample, "fastq1"],
    "fq2": samples_table.loc[wildcards.sample, "fastq2"]
  }

# this makes a sample table with "_merged" appended to the end of each merged sample name in the "merged_sample" and "Control" columns
def add_merge_suffix_to_merged_samples(samples_df):
  df=samples_df
  merged_list=df['merged_sample'].to_list()
  def unique(list1):
    list_set = set(list1)
    unique_list = (list(list_set))
    return unique_list
  unique_merged_list=unique(merged_list)
  unique_merged_list_no_nan=[x for x in unique_merged_list if x != 'nan']
  Control_list=df['Control'].to_list()
  for mrged in unique_merged_list_no_nan:
    Control_list = [re.sub(r"\b{}\b".format(mrged), mrged + "_merged", s) for s in Control_list]
    merged_list = [re.sub(r"\b{}\b".format(mrged), mrged + "_merged", s) for s in merged_list]
  df2 = df
  df2['Control'] = Control_list
  df2['merged_sample'] = merged_list
  return df2

samples_table_w_merged_suffix = add_merge_suffix_to_merged_samples(samples_table)

# sample_type input function definition set to Python dictionary
def sample_type_dict_from_sample(wildcards):
  return {
    "treatment": 'results/aligned_speciesOfInterest/' + all_treatments_table.loc[wildcards.sample, "sample"] + '.bam',
    "control": 'results/aligned_speciesOfInterest/' + all_treatments_table.loc[wildcards.sample, "Control"] + '.bam'
  }

def keywords_to_merge(table=samples_table_w_merged_suffix):
  samples_table3 = table[~table['merged_sample'].isna()]
  samples_table4 = samples_table3.drop_duplicates(subset="merged_sample",keep='first')
  lst1 = samples_table4['merged_sample'].to_list()
  lst2 = [x for x in lst1 if x != 'nan']
  return lst2

merged_keywords_lst = keywords_to_merge()

def get_bams_to_merge(smpl):
  merged_list=samples_table_w_merged_suffix['merged_sample'].to_list()
  samples_list=samples_table_w_merged_suffix['sample'].to_list()
  indices=[i for i,x in enumerate(merged_list) if x==smpl]
  samples_list2 = [samples_list[i] for i in indices]
  samples_list3 = ['results/aligned/' + s + '.bam' for s in samples_list2]
  samples_string = ' '.join([str(item) for item in samples_list3])
  return samples_string

bams_to_merge_dict = {}
for m in merged_keywords_lst:
  bams_to_merge_dict[m] = get_bams_to_merge(m)

all_samples_lst = merged_keywords_lst + samples_table['sample'].to_list()

# this makes a new sample table with only the 'treatment' sample rows
treatments_table = samples_table_w_merged_suffix.loc[samples_table_w_merged_suffix['sampleType'] == 'treatment']
# this makes a treatments table with the merged samples in the "samples" column
merged_treatments_table = treatments_table.copy()
merged_treatments_table.loc[:,'sample'] = merged_treatments_table['merged_sample'].to_list()
merged_treatments_table = merged_treatments_table.set_index('merged_sample')
all_treatments_table = pd.concat([treatments_table, merged_treatments_table])
dup_series=all_treatments_table.duplicated(subset='sample',keep='first')
all_treatments_table = all_treatments_table[~dup_series]
all_treatments_table = all_treatments_table[all_treatments_table["sample"].str.contains("nan")==False]
treatment_samples_lst = all_treatments_table['sample'].to_list()

# get sample names for those to be merged
samples_to_merge = samples_table_w_merged_suffix[samples_table_w_merged_suffix["merged_sample"].str.contains("nan")==False]
samples_to_merge_lst = samples_to_merge['sample'].to_list()
