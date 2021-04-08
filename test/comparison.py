import yaml, json
import os

# os.system('cd ..')
# os.system('pwd')

# configs

# for hg38 test samples
config_dir=os.path.join("..","results","BSA_0000_test_hg38_atac","config_files","BSA_0000_test_hg38_atac.inputs.json")
result_dir1=os.path.join("..","results","BSA_0000_test_hg38_atac","atacseq_results")
result_dir2=os.path.join(os.sep,"nobackup","lab_bsf","users","berguener","junk","BSA_0000_test_atac","atacseq_results")

# for mm10 test samples
# config_dir=os.path.join("..","results","BSA_0000_test_mm10_atac","config_files","BSA_0000_test_mm10_atac.inputs.json")
# result_dir1=os.path.join("..","results","BSA_0000_test_mm10_atac","atacseq_results")
# result_dir2=os.path.join(os.sep,"nobackup","lab_bsf","users","berguener","junk","BSA_0000_test_atac","old","atacseq_results")

# QC stat file extensions
stats = [["",".stats.tsv"],
         ["",".tss_histogram.csv"],
         ["mapped",".samtools_flagstat.log"],
         ["mapped",".txt"]
        ]


# load config
with open(config_dir, 'r') as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exception:
        sys.stderr.write(str(exception))
# print(config)

# get sample names
samples = config["atacseq.sample_list"]

# iterate over samples
for sample in samples:
    for stat in stats:
        path_1=os.path.join(result_dir1,sample,"{}".format(stat[0]),"{}{}".format(sample,stat[1]))
        path_2=os.path.join(result_dir2,sample,"{}".format(stat[0]),"{}{}".format(sample,stat[1]))
        diff = os.system("diff {} {}".format(path_1, path_2))
        if diff==0:
            print("sample {} had no differences in statistic {}".format(sample, stat[1]))
        else:
            print("sample {} had differences in statistic {}!".format(sample, stat[1]))
    #     if diff>1:
    #         print("sample {} comparison failed!".format(sample))

