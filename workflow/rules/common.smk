
def get_raw_bams(wildcards):
    return annot.loc[wildcards.sample, "bam_file"]
#     return str.split(samples[wildcards.sample]["raw_bams"])