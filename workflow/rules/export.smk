# export and add all used conda environment specifications (including exact versions and builds) to report 
rule env_export:
    output:
        report(os.path.join(result_path,'envs','{env}.yaml'),
                      caption="../report/software.rst", 
                      category="Software", 
                      subcategory="{}_{}".format(config["project_name"], module_name),
                       labels={
                                    "name": config["project_name"],
                                    "module": module_name,
                                   "env": "{env}",
                                }
                     ),
    conda:
        "../envs/{env}.yaml"
    resources:
        mem_mb=config.get("mem", "1000"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","env_{env}.log"),
    shell:
        """
        conda env export > {output}
        """
        
# export and add configuration file to report        
rule config_export:
    output:
        configs = report(os.path.join(result_path,'configs','{}_config.yaml'.format(config["project_name"])), 
                         caption="../report/configs.rst", 
                         category="Configuration", 
                         subcategory="{}_{}".format(config["project_name"], module_name),
                         labels={
                                    "name": config["project_name"],
                                    "module": module_name,
                                     "type": "config"
                                }
                        )
    resources:
        mem_mb=config.get("mem", "1000"),
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","config_export.log"),
    run:
        with open(output["configs"], 'w') as outfile:
            yaml.dump(config, outfile, sort_keys=False, width=1000, indent=2)

# export and add used annotation file(s) to report
rule annot_export:
    input:
        config["annotation"],
    output:
        annot = report(os.path.join(result_path,'configs','{}_annot.csv'.format(config["project_name"])), 
                         caption="../report/configs.rst", 
                         category="Configuration", 
                         subcategory="{}_{}".format(config["project_name"], module_name),
                       labels={
                                    "name": config["project_name"],
                                    "module": module_name,
                                   "type": "annotation",
                                }
                        )
    resources:
        mem_mb=1000,
    threads: config.get("threads", 1)
    log:
        os.path.join("logs","rules","annot_export.log"),
    shell:
        """
        cp {input} {output}
        """