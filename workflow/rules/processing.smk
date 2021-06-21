# check if splits are defined, if yes split data
if len(data_splits)>1:
    # split dataset into subsets by split_by parameter
    rule split_data:
        input:
            counts = os.path.join(config["project_path"],'all',"all_counts.csv"),
            annotation_filtered = os.path.join(config["project_path"],'all',"all_annotation.csv"),
        output:
            split_data = expand(os.path.join(config["project_path"],'{split}', '{split}_counts.csv'),split=data_splits[1:]),
            split_annot = expand(os.path.join(config["project_path"],'{split}', '{split}_annotation.csv'),split=data_splits[1:]),
        params:
            split_by = config["split_by"],
            project_path=config["project_path"],
            # cluster parameters
            partition=partition,
        threads: threads
        resources:
            mem=mem,
        conda:
            "../envs/atacseq_analysis.yaml",
        log:
            "logs/rules/split_data.log"
        script:
            "../scripts/split_data.py"

        
# filter regions
rule filter_regions:
    input:
        counts = os.path.join(config["project_path"],'{split}',"{split}_counts.csv"),
        support = os.path.join(config["project_path"],'all',"all_support.csv"),
        annot = os.path.join(config["project_path"],'{split}',"{split}_annotation.csv"),
        consensus_regions = os.path.join(config["project_path"],'all',"consensus_regions.bed"),
        regions_annot = os.path.join(config["project_path"],'all',"consensus_regions_annotation.csv"),
    output:
        filtered_data=os.path.join(config["project_path"],'{split}',"{split}_filtered.csv"),
        filtered_plots=report(os.path.join(config["project_path"],'{split}',"{split}_filtered_regions.svg"), caption="../report/filter_regions.rst", category="{split}", subcategory="filtered"),
        filtered_consensus_regions=os.path.join(config["project_path"],'{split}',"{split}_consensus_regions_filtered.bed"),
        filtered_regions_annot=os.path.join(config["project_path"],'{split}',"{split}_consensus_regions_annotation_filtered.csv"),
    params:
        peak_support_threshold = config["peak_support_threshold"],
        proportion = config["proportion"],
        min_group = config["min_group"],
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/filter_regions_{split}.log"
    script:
        "../scripts/filter_regions.py"
        
        
# normalize data via TMM (yields log2 CPM values)
rule normalize_TMM:
    input:
        filtered_data=os.path.join(config["project_path"],'{split}',"{split}_filtered.csv"),
    output:
        normTMM_data=os.path.join(config["project_path"],'{split}',"{split}_normTMM.csv"),
    params:
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/normalize_TMM_{split}.log"
    script:
        "../scripts/normalize_TMM.py"


# normalize data via CQN
rule normalize_CQN:
    input:
        filtered_data=os.path.join(config["project_path"],'{split}',"{split}_filtered.csv"),
        consensus_regions = os.path.join(config["project_path"],'all',"consensus_regions.bed"),
    output:
        normCQN_data=os.path.join(config["project_path"],'{split}',"{split}_normCQN.csv"),
    params:
        genome_fasta=config["genome_fasta"],
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/normalize_CQN_{split}.log"
    script:
        "../scripts/normalize_CQN.py"


# select HVR
rule select_HVR:
    input:
        norm_data=os.path.join(config["project_path"],'{split}',"{split}_{norm}.csv"),
    output:
        HVR_data=os.path.join(config["project_path"],'{split}',"{split}_{norm}_HVR.csv"),
        HVR_plot=report(os.path.join(config["project_path"],'{split}',"{split}_{norm}_HVR_selection.svg"), caption="../report/region_variability.rst", category="{split}", subcategory="{norm}_HVR"),
    params:
        HVR_percentage = config["HVR_percentage"],
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/select_HVR_{split}_{norm}.log"
    script:
        "../scripts/select_HVR.py"


# plot mean-variance relationship for every split and step
rule plot_mean_var:
    input:
        data=os.path.join(config["project_path"],'{split}',"{split}_{step}.csv"),
    output:
        plot=report(os.path.join(config["project_path"],'{split}',"mean_variance_analysis","mean_variance_{split}_{step}.svg"), caption="../report/meanvar.rst", category="{split}", subcategory="{step}"),
    params:
        results_dir=lambda w: os.path.join(config["project_path"],'{}'.format(w.split),"mean_variance_analysis"),
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/plot_mean_var_{split}_{step}.log"
    script:
        "../scripts/plot_mean_var.py"