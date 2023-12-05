
def get_raw_bams(wildcards):
    return annot.loc[wildcards.sample, "bam_file"]

def get_quantifications(wildcards):

    if wildcards.kind=="support":
        paths = expand(os.path.join(result_path, "results", "{sample}", "peaks", "{sample}_quantification_support_counts.csv"), sample=samples_quantify)

    if wildcards.kind=="consensus":
        paths = expand(os.path.join(result_path, "results", "{sample}", "mapped", "{sample}_quantification_consensus_counts.csv"), sample=samples_quantify)

    if wildcards.kind=="promoter":
        paths = expand(os.path.join(result_path, "results", "{sample}", "mapped", "{sample}_quantification_promoter_counts.csv"), sample=samples_quantify)

    return paths