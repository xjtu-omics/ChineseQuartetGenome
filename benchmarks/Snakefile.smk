# ======================================================================================================================
# Project: ChinseQuartetGenome
# Script : Snakefile.smk TODO check 
# Author : Peng Jia
# Date   :  2023/2/18
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================
configfile: "conf/config.yaml"
include: "rules/preprocessing.smk"
include: "rules/SV_truvari.smk"


targets = expand(dir_work + "results/{prefix}",prefix=benchmark_prefix_ids[0])

# targets=expand(dir_work + "benchmarks/{ref}/{ref}.{region}_{base_region}.bed",ref="GRCh38",region=["SD","RM","SR","VNTR","STR"],base_region=[])

rule target:
    input:
        targets
