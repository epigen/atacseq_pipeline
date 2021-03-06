# performs UMAP dimensionality reduction on the data 
rule do_UMAP:
    input:
        quant_counts=os.path.join(config["project_path"],'{split}',"{split}_{step}.csv"),
    output:
        umap_data=os.path.join(config["project_path"],'{split}', 'unsupervised_analysis',"UMAP_{split}_{step}.csv"),
    params:
        results_dir=os.path.join(config["project_path"],'{split}', 'unsupervised_analysis'),
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/do_UMAP_{split}_{step}.log"
    script:
        "../scripts/do_UMAP.py"
        
# performs PCA on the data 
rule do_PCA:
    input:
        quant_counts=os.path.join(config["project_path"],'{split}',"{split}_{step}.csv"),
    output:
        pca_data=os.path.join(config["project_path"],'{split}', 'unsupervised_analysis',"PCA_{split}_{step}.csv"),
        pca_expl_var=os.path.join(config["project_path"],'{split}', 'unsupervised_analysis',"PCA_{split}_{step}_explained_variance.csv"),
    params:
        results_dir=os.path.join(config["project_path"],'{split}', 'unsupervised_analysis'),
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/do_PCA_{split}_{step}.log"
    script:
        "../scripts/do_PCA.py"

rule plot_dimred:
    input:
        annotation_filtered = os.path.join(config["project_path"],'{split}',"{split}_annotation.csv"),
        dimred_data=os.path.join(config["project_path"],'{split}', 'unsupervised_analysis',"{dimred}_{split}_{step}.csv"),
    output:
        dimred_plots=report(expand(os.path.join(config["project_path"],'{{split}}', 'unsupervised_analysis',"{{dimred}}_{{split}}_{{step}}_{plot_by_variable}.svg"),plot_by_variable=config['plot_by']), caption="../report/dimred.rst", category="{split}", subcategory="{step}"),
    params:
        variables = config['plot_by'],
        label=lambda w: "{}_{}".format(w.split, w.step),
        results_dir=os.path.join(config["project_path"],'{split}', 'unsupervised_analysis'),
        dimred=lambda w: "{}".format(w.dimred),
        # cluster parameters
        partition=partition,
    threads: threads
    resources:
        mem=mem,
    conda:
        "../envs/atacseq_analysis.yaml",
    log:
        "logs/rules/plot_dimred_{dimred}_{split}_{step}.log"
    script:
        "../scripts/plot_dimred.py"