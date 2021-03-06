# install homer        
rule install_homer:
    output:
        homer_script=os.path.join("resources","homer","configureHomer.pl"),
    params:
        partition=partition,
    threads: threads
    resources:
        mem=mem,
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