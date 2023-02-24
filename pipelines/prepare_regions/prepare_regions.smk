# ======================================================================================================================
# Project: ChinseQuartetGenome
# Script : prepare_regions.smk
# Author : Peng Jia
# Date   :  2023/2/21
# Email  : pengjia@stu.xjtu.edu.cn
# Description: Prepare regions for ChineseQuartet benchmarks
# ======================================================================================================================
import logging
import os.path
import yaml

try:
    from yaml import CDumper as Dumper
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Dumper
    from yaml import Loader
import pandas as pd

configfile: "config.yaml"
bedtools = config["software"]["bedtools"] if "bedtools" in config["software"] else "bedtools"
tabix = config["software"]["tabix"] if "tabix" in config["software"] else "tabix"
gzip = config["software"]["gzip"] if "gzip" in config["software"] else "gzip"
bgzip = config["software"]["bgzip"] if "bgzip" in config["software"] else "bgzip"
dir_work = config["dir_work"] if "dir_work" in config else os.path.abspath("./benchmark_regions")
dir_work = dir_work if dir_work[-1] == "/" else dir_work + "/"
regions_info = config["regions"]

# wildcard_constraints:
#     sample="|".join(samples_ids),
#     ref="|".join(ref_ids),
#     region="|".join(regions_ids),
#     var_type="|".join(["SNV", "Indel", "DEL", "CSV", "INV", "INS"])
targets_dict = {}
targets = []
for ref, info in config["regions_pattern"].items():
    targets_dict[ref] = {}
    for pattern in info:
        if pattern in targets_dict[ref]: continue
        targets_dict[ref][pattern] = dir_work + f"{ref}/{ref}.{pattern}.bed"
        targets.append(dir_work + f"{ref}/{ref}.{pattern}.bed")

# print(targets_dict)#
#
# targets = [ for ref, info in config["regions_pattern"].items() for pattern in info]
# targets_dict = {dir_work + f"{ref}/{ref}.{pattern}.bed" for ref, info in config["regions_pattern"].items() for pattern in info}
rule all:
    input:
        targets,
        dir_work + "regions.yaml"


# expand(dir_work + "{ref}/{ref}.{pattern}.bed",ref="GRCh38",pattern=config["regions_pattern"]["GRCh38"])
rule get_final_config:
    input:
        targets
    output:
        dir_work + "regions.yaml"
    run:
        with open(f"{output}","w") as yamlfile:
            yaml.dump({"regions": targets_dict},yamlfile,)


# print(open(f"{output}").readlines())


def get_region_input_from_input_config(wildcards):
    return regions_info[wildcards.ref][wildcards.region]


rule get_regions:
    input:
        get_region_input_from_input_config
    output:
        dir_work + "{ref}/{ref}.{region}.bed"
    run:
        if wildcards.region == "STR":
            if f"{input}"[-2:] == "gz":
                shell("{bgzip} -d -c {input} >{output}_tmp")
                shell("""
                {bedtools} sort -i {output}_tmp| awk -v OFS="\t" '{{print $1,$2-5,$3+5 }}' |{bedtools} merge > {output}
                """)
                shell("rm -rf {output}_tmp ")
            else:
                # shell("cp {input} {output}")
                shell("""
                {bedtools} sort -i  {input}| awk -v OFS="\t" '{{print $1,$2-5,$3+5 }}'> {output}
                """)
        else:
            if f"{input}"[-2:] == "gz":
                shell("{bgzip} -d -c {input} >{output}_tmp")
                shell("{bedtools} sort -i {output}_tmp |{bedtools} merge > {output}")
                shell("rm -rf {output}_tmp ")
            else:
                # shell("cp {input} {output}")
                shell("{bedtools} sort -i  {input} |{bedtools} merge > {output}")

# shell("rm -rf {output}_tmp ")


# def get_base_bench(wildcards):
#     if wildcards.base_region == "ALL":
#         return config["regions"][wildcards.ref]["benchmark_regions"]
#     elif wildcards.base_region == "HighConfSV":
#         return config["regions"][wildcards.ref]["high_confidence_regions_SV"]
#     elif wildcards.base_region == "HighConfSmall":
#         return config["regions"][wildcards.ref]["high_confidence_regions_Small"]


rule get_regions_pattern_include:
    input:
        include=dir_work + "{ref}/{ref}.{region_base}.bed",
        exclude=dir_work + "{ref}/{ref}.{region}.bed"
    output:
        dir_work + "{ref}/{ref}.{region_base}_include_{region}.bed",
    run:
        shell("{bedtools} intersect -a {input.include} -b {input.exclude} |awk '{{ if ($3-$2>10) print }}' >{output}")

rule get_regions_pattern_exclude:
    input:
        include=dir_work + "{ref}/{ref}.{region_base}.bed",
        exclude=dir_work + "{ref}/{ref}.{region}.bed"
    output:
        dir_work + "{ref}/{ref}.{region_base}_exclude_{region}.bed",
    run:
        shell("{bedtools} subtract -a {input.include} -b {input.exclude} |awk '{{ if ($3-$2>10) print }}' >{output}")
