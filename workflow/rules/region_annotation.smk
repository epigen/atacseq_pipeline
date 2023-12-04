# annotate consensus regions

# prepare configs for uropa    
rule uropa_prepare:
    input:
        consensus_regions = os.path.join(result_path,"counts","consensus_regions.bed"),
        gencode_template = os.path.join("config","uropa","gencode_config_TEMPLATE_V4.txt"),
        reg_template = os.path.join("config","uropa","regulatory_config_TEMPLATE.txt"),
    output:
        gencode_config = os.path.join(result_path,"tmp","consensus_regions_gencode.json"),
        reg_config = os.path.join(result_path,"tmp","consensus_regions_reg.json"),
    params:
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    log:
        "logs/rules/uropa_prepare.log"
    run:
        ### generate gencode config
        with open(input.gencode_template) as f:
            gencode_template=Template(f.read())

        gencode_config=gencode_template.substitute({
                'TSS_flanking':config['tss_size'],
                'TSS_proximal_upstream':config['proximal_size_up'],
                'TSS_proximal_downstream':config['proximal_size_dn'],
                'distal_distance':config['distal_size'],
                'gtf_file':'"{}"'.format(config["gencode_gtf"]),
                'bed_file':'"{}"'.format(input.consensus_regions)
            })

        with open(output.gencode_config,'w') as out:
            out.write(gencode_config)

        ### generate reg config file
        with open(input.reg_template) as f:
            reg_template=Template(f.read())  

        reg_config=reg_template.substitute({
            'gtf_file':'"{}"'.format(config["regulatory_build_gtf"]),
            'bed_file':'"{}"'.format(input.consensus_regions)
        })

        with open(output.reg_config,'w') as out:
            out.write(reg_config)

# run uropa on consensus regions for gencode
rule uropa_gencode:
    input:
        consensus_regions = os.path.join(result_path,"counts","consensus_regions.bed"),
        gencode_config = os.path.join(result_path,"tmp","consensus_regions_gencode.json"),
    output:
        gencode_results = os.path.join(result_path,"tmp","gencode_finalhits.txt"),
    params:
        # paths
        results_dir = os.path.join(result_path,"tmp"),
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: 4*config.get("threads", 2)
    conda:
        "../envs/uropa.yaml",
    log:
        "logs/rules/uropa_run_gencode.log"
    shell:
        """
        uropa -p {params.results_dir}/gencode -i {input.gencode_config} -t {threads} -l {params.results_dir}/uropa.gencode.log
        """

# run uropa on consensus regions for regulatory build
rule uropa_reg:
    input:
        consensus_regions = os.path.join(result_path,"counts","consensus_regions.bed"),
        reg_config = os.path.join(result_path,"tmp","consensus_regions_reg.json"),
    output:
        reg_results = os.path.join(result_path,"tmp","reg_finalhits.txt"),
    params:
        # paths
        results_dir = os.path.join(result_path,"tmp"),
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 8)
    conda:
        "../envs/uropa.yaml",
    log:
        "logs/rules/uropa_run_reg.log"
    shell:
        """
        uropa -p {params.results_dir}/reg -i {input.reg_config} -t {threads} -l {params.results_dir}/uropa.reg.log
        """

# peak annotation using homer
rule homer_region_annotation:
    input:
        consensus_regions = os.path.join(result_path,"counts","consensus_regions.bed"),
        homer_script = os.path.join(HOMER_path,"configureHomer.pl"),
    output:
        homer_annotations = os.path.join(result_path,"tmp","homer_annotations.tsv"),
        homer_annotations_log = os.path.join(result_path,"tmp","homer_annotations.tsv.log"),
    params:
        # paths
        homer_bin = os.path.join(HOMER_path,"bin"),
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    conda:
        "../envs/macs2_homer.yaml",
    log:
        "logs/rules/homer_region_annotation.log"
    shell:
        """
        export PATH="{params.homer_bin}:$PATH";
        
        {params.homer_bin}/annotatePeaks.pl {input.consensus_regions} {config[genome]} \
            > {output.homer_annotations} \
            2> {output.homer_annotations_log};
        """
        
# get gc content and region length
rule bedtools_annotation:
    input:
        consensus_regions = os.path.join(result_path,"counts","consensus_regions.bed"),
    output:
        bedtools_annotation = os.path.join(result_path, "tmp", "bedtools_annotation.bed"),
    params:
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    conda:
        "../envs/pybedtools.yaml",
    log:
        "logs/rules/bedtools_annotation.log"
    shell:
        """
        bedtools nuc -fi {config[genome_fasta]} -bed {input.consensus_regions} > {output.bedtools_annotation}
        """
        
# aggregate uropa and homer annotation results
rule region_annotation_aggregate:
    input:
        gencode_results = os.path.join(result_path,"tmp","gencode_finalhits.txt"),
        reg_results = os.path.join(result_path,"tmp","reg_finalhits.txt"),
        homer_annotations = os.path.join(result_path,"tmp","homer_annotations.tsv"),
        bedtools_annotation = os.path.join(result_path, "tmp", "bedtools_annotation.bed"),
    output:
        region_annotation = os.path.join(result_path,'counts',"region_annotation.csv"),
    params:
        # cluster parameters
        partition=config.get("partition"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 2)
    log:
        "logs/rules/region_annotation_aggregate.log"
    run:
        # load and format uropa gencode results
        gencode_characterization = pd.read_csv(input.gencode_results, sep='\t')
        gencode_characterization = gencode_characterization.set_index("peak_id")
        gencode_characterization.loc[gencode_characterization['feature']=='transcript','feat_type']='transcript:'+gencode_characterization.loc[gencode_characterization['feature']=='transcript','transcript_type']
        gencode_characterization.loc[gencode_characterization['feature']=='gene','feat_type']='gene:'+gencode_characterization.loc[gencode_characterization['feature']=='gene','gene_type']
        gencode_characterization['length']=gencode_characterization['peak_end']-gencode_characterization['peak_start']
        gencode_characterization=gencode_characterization[['peak_chr','peak_start','peak_end','length','feat_anchor','distance','relative_location','feat_type','gene_id','gene_name','name']]
        gencode_characterization.columns=['chr','start','end','length','feat_anchor','distance','location','feat_type','gene_id','gene_name','characterization']
        gencode_characterization.loc[gencode_characterization['characterization'].isna(),'characterization']='NONE'
        gencode_characterization=gencode_characterization.add_prefix('gencode_')

        # load and format uropa regulatory build results
        reg_characterization=pd.read_csv(input.reg_results,sep='\t')
        reg_characterization = reg_characterization.set_index('peak_id')[['feature','ID']]
        reg_characterization.columns=['reg_feature','reg_feature_id']
        reg_characterization.loc[reg_characterization['reg_feature'].isna(),'reg_feature']='reg_NONE'
        reg_characterization=reg_characterization.add_prefix('regulatoryBuild_')
        
        # load and format homer annotation results
        homer_annotation = pd.read_csv(input.homer_annotations,sep='\t', index_col=0)
        homer_annotation = homer_annotation[['Annotation','Detailed Annotation','Distance to TSS','Nearest PromoterID','Entrez ID','Nearest Unigene','Nearest Refseq','Nearest Ensembl','Gene Name','Gene Alias','Gene Description','Gene Type']]
        homer_annotation = homer_annotation.add_prefix('homer_')
        
        # load and format bedtools annotation results
        bedtools_annotation = pd.read_csv(input.bedtools_annotation, sep='\t', index_col = 3)
        bedtools_annotation = bedtools_annotation.iloc[:,3:]
        bedtools_annotation.columns = [col.split('_', 1)[-1].replace('at', 'AT').replace('gc', 'GC').replace('oth', 'otherBases') for col in bedtools_annotation.columns]
        bedtools_annotation = bedtools_annotation.add_prefix('bedtools_')

        # join results
        base_character = gencode_characterization.join(homer_annotation)
        base_character = base_character.join(reg_characterization)
        base_character = base_character.join(bedtools_annotation)
        
        # replace whiteapaces in colnames with underscore
        base_character.columns = base_character.columns.str.replace(' ', '_')
        
        # save final results
        base_character.to_csv(output.region_annotation, index_label='peak_id')
        
        # remove tmp folder
        tmp_dir = os.path.dirname(input.gencode_results)
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)