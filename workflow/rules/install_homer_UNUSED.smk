# install homer        
rule install_homer:
    output:
        homer_script=os.path.join("resources","homer","configureHomer.pl"),
    params:
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "8000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/homer.yaml",
    log:
        "logs/rules/install_homer.log"
    shell:
        """
        mkdir -p resources/homer
        cd resources/homer
        wget http://homer.ucsd.edu/homer/configureHomer.pl
        perl configureHomer.pl -install
        perl configureHomer.pl -install hg38 mm10
        """