

# rule plot_dimred:
#     input:
#         annotation = config["atacseq.annotation_metadata"],
#         quant_counts=os.path.join(config["atacseq.project_path"],'all',"counts_all.csv"),
#     output:
#         umap_plots=report(expand(os.path.join(config["atacseq.project_path"],'all', 'unsupervised_analysis',"UMAP_"+config['atacseq.project_name']+"_{plot_by_variable}.svg"),plot_by_variable=config['atacseq.plot_by']), caption="../report/UMAP.rst", category="Unsupervised Analysis"),
#         pca_plots=report(expand(os.path.join(config["atacseq.project_path"],'all', 'unsupervised_analysis',"PCA_"+config['atacseq.project_name']+"_{plot_by_variable}.svg"),plot_by_variable=config['atacseq.plot_by']), caption="../report/PCA.rst", category="Unsupervised Analysis"),
#     params:
#         variables = config['atacseq.plot_by'],
#         label=config['atacseq.project_name'],
#         results_dir=config["atacseq.project_path"],
#         # cluster parameters
#         partition=partition,
#     threads: threads
#     resources:
#         mem=mem,
#     conda:
#         "../envs/atacseq_analysis.yaml",
#     log:
#         "logs/rules/plot_dimred.log"
#     script:
#         "../scripts/plot_dimred.py"

# performs UMAP dimensionality reduction on the data 
rule do_UMAP:
    input:
        quant_counts=os.path.join(config["atacseq.project_path"],'all',"counts_all.csv"),
    output:
        umap_data=os.path.join(config["atacseq.project_path"],'all', 'unsupervised_analysis',"UMAP_"+config['atacseq.project_name']+".csv"),
    params:
        results_dir=os.path.join(config["atacseq.project_path"],'all', 'unsupervised_analysis'),
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/do_UMAP.log"
    script:
        "../scripts/do_UMAP.py"
        
# performs PCA on the data 
rule do_PCA:
    input:
        quant_counts=os.path.join(config["atacseq.project_path"],'all',"counts_all.csv"),
    output:
        pca_data=os.path.join(config["atacseq.project_path"],'all', 'unsupervised_analysis',"PCA_"+config['atacseq.project_name']+".csv"),
    params:
        results_dir=os.path.join(config["atacseq.project_path"],'all', 'unsupervised_analysis'),
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/do_PCA.log"
    script:
        "../scripts/do_PCA.py"

rule plot_dimred:
    input:
        annotation = config["atacseq.annotation_metadata"],
        dimred_data=os.path.join(config["atacseq.project_path"],'all', 'unsupervised_analysis',"{dimred}_"+config['atacseq.project_name']+".csv"),
    output:
        dimred_plots=report(expand(os.path.join(config["atacseq.project_path"],'all', 'unsupervised_analysis',"{{dimred}}_"+config['atacseq.project_name']+"_{plot_by_variable}.svg"),plot_by_variable=config['atacseq.plot_by']), caption="../report/dimred.rst", category="Unsupervised Analysis"),
    params:
        variables = config['atacseq.plot_by'],
        label=config['atacseq.project_name'],
        results_dir=os.path.join(config["atacseq.project_path"],'all', 'unsupervised_analysis'),
        dimred=lambda w: "{}".format(w.dimred),
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/plot_dimred_{dimred}.log"
    script:
        "../scripts/plot_dimred.py"