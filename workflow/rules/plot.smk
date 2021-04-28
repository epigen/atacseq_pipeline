

rule plot_dimred:
    input:
        annotation = config["atacseq.annotation_metadata"],
        quant_counts=os.path.join(config["atacseq.project_path"],"quantification_counts.csv"),
    output:
        umap_plots=report(expand(os.path.join(config["atacseq.project_path"],"UMAP_"+config['atacseq.project_name']+"_{plot_by_variable}.svg"),plot_by_variable=config['atacseq.plot_by']), caption="../report/UMAP.rst", category="Unsupervised Analysis"),
        pca_plots=report(expand(os.path.join(config["atacseq.project_path"],"PCA_"+config['atacseq.project_name']+"_{plot_by_variable}.svg"),plot_by_variable=config['atacseq.plot_by']), caption="../report/PCA.rst", category="Unsupervised Analysis"),
    params:
        variables = config['atacseq.plot_by'],
        label=config['atacseq.project_name'],
        results_dir=config["atacseq.project_path"],
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/plot_dimred.log"
    script:
        "../scripts/plot_dimred.py"