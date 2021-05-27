#!/usr/bin/env python

import pandas as pd
import os
import sys
import csv
import pybedtools as bedtools

reg_build_file=sys.argv[1]
# chrom_file=sys.argv[2]

print("Parsing file {}".format(reg_build_file))

folder=os.path.dirname(reg_build_file)
basename=os.path.basename(reg_build_file).replace(".gff.gz","")

std_chroms=['{}'.format(i) for i in range(1,23)]
std_chroms.extend(['X','Y'])

reg_build = pd.read_csv(reg_build_file,sep='\t',
                        names=['chrom','build','type','start','end','score','strand','frame','annotations'],
                        dtype={'chrom':str})

reg_build = reg_build.loc[reg_build['chrom'].isin(std_chroms)]
reg_build = reg_build.sort_values(by=['chrom','start','end'])
reg_build['ENSEMBL_ID']=reg_build['annotations'].map(
    lambda x : x.split(";")[0].split("=")[1].split(":")[1]
)

reg_build['OUT_ID']="ID \""+reg_build['ENSEMBL_ID']+"\""
reg_build[['chrom','build','type','start','end','score','strand','frame','OUT_ID']].to_csv(
    os.path.join(folder,"{}.gtf".format(basename)),
    sep='\t',
    index=False,
    header=False,
    quoting=csv.QUOTE_NONE
)


# reg_build["ID"] = reg_build[["chrom","start","end"]].apply(tuple,axis=1)\
#     .rank(method='dense').astype(int).map(lambda x: "REG{:010d}".format(x))
# reg_build=reg_build[['chrom','start','end','type','ID','ENSEMBL_ID']]
# reg_build.to_csv(os.path.join(folder,"{}.parsed.csv".format(basename)),index=False)

# reg_build = reg_build[['chrom','start','end','ID']].drop_duplicates()
# reg_build['chrom']='chr'+reg_build['chrom']
# reg_build['score']='.'
# reg_build['strand']='.'
# consensus_bed=bedtools.BedTool().from_dataframe(reg_build[['chrom','start','end','ID','score','strand']])\
#     .sort(faidx=chrom_file)\
#     .saveas(os.path.join(folder,"{}.parsed.bed".format(basename)))