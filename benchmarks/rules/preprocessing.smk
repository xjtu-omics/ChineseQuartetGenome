# ======================================================================================================================
# Project: ChinseQuartetGenome
# Script : preprocessing.smk TODO check 
# Author : Peng Jia
# Date   :  2023/2/21
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================
import logging as pylogging
import os.path

import pysam

pylogging.basicConfig(level=pylogging.ERROR)
import pandas as pd
# logging.setLevel(logging.DEBUG)
bedtools = config["software"]["bedtools"] if "bedtools" in config["software"] else "bedtools"
tabix = config["software"]["tabix"] if "tabix" in config["software"] else "tabix"
gzip = config["software"]["gzip"] if "gzip" in config["software"] else "gzip"
bgzip = config["software"]["bgzip"] if "bgzip" in config["software"] else "bgzip"
hap_eval = config["software"]["hap_eval"] if "hap_eval" in config["software"] else "hap.py"
truvari = config["software"]["truvari"] if "truvari" in config["software"] else "truvari"
dir_work = config["dir_work"] if "dir_work" in config else os.path.abspath("./ChineseQuartet_benchmarks")
dir_work = dir_work if dir_work[-1] == "/" else dir_work + "/"
sample_info_data = pd.read_table(config["samples"])
sample_infos = {}
benchmark_prefix_ids = []
samples_ids = []
ref_ids = []
regions_ids = []
methods_ids = []
version_ids = []
comparison_tools = {
    var_type: [] for var_type in ["DEL", "INS", "Indel", "SNV", "CSV"]
}
for var_type, info in config["comparison_tools"].items():
    if var_type in ["SV"]:
        comparison_tools["DEL"].append(info)
        comparison_tools["INS"].append(info)
        continue
    if var_type in ["Small"]:
        comparison_tools["SNV"].append(info)
        comparison_tools["Indel"].append(info)
        continue
    comparison_tools[var_type].append(info)

for vcf_id, info in sample_info_data.iterrows():

    this_regions = []
    if info["reference_version"] not in config["variants"]:
        pylogging.error(f"No match variant was provided for this version of reference ({info['reference_version']}), please check your config files")
        exit()
    if info["var_type"] not in config["variants"][info["reference_version"]]:
        pylogging.error(f"No match variant was provided for this varaint type ({info['var_type']}), please check your config files")
        exit()
    if info["reference_version"] not in config["regions"]:
        pylogging.warn(f"No match benchmarking region was provided for this version of reference ({info['reference_version']}), please check your config files."
                     f"Now, all variants on the whole genomes will be evaluated!")
    else:
        if "benchmark_regions" not in info["reference_version"]:
            pylogging.warn(f"No benchmark regions was provided for this version of reference ({info['reference_version']}), please check your config files")

        for r in config["regions"][info["reference_version"]]:
            regions_ids.append(r)
            if ("SV" in r) and (info["var_type"] in ["DEL", "INS"]):
                continue
            if ("Small" in r) and (info["var_type"] in ["SNV", "Indel"]):
                continue
            this_regions.append(r)
    if info["reference_version"] not in sample_infos:
        sample_infos[info["reference_version"]] = {}
    if info["sample_id"] not in sample_infos[info["reference_version"]]:
        sample_infos[info["reference_version"]][info["sample_id"]] = {}
    if info["var_type"] not in sample_infos[info["reference_version"]][info["sample_id"]]:
        sample_infos[info["reference_version"]][info["sample_id"]][info["var_type"]] = {}

    this_version_id = f"{info['method']}"
    if this_version_id in sample_infos[info["reference_version"]][info["sample_id"]][info["var_type"]]:
        pylogging.error(f"There are same version of input file: {this_version_id}, please check you input file!")
        exit()
    else:
        sample_infos[info["reference_version"]][info["sample_id"]][info["var_type"]][this_version_id] = info["path_vcf"]
        line_id = f"{info['reference_version']}.{info['sample_id']}.{info['var_type']}.{this_version_id}"
        benchmark_prefix_ids.extend([f"{line_id}.{i}/{vv}/{line_id}.{i}" for i in this_regions for vv in
                                     comparison_tools[info["var_type"]]])
    samples_ids.append(info["sample_id"])
    ref_ids.append(info["reference_version"])
    methods_ids.append(info["method"])

wildcard_constraints:
    sample="|".join(samples_ids),
    ref="|".join(ref_ids),
    region="|".join(regions_ids),
    method="|".join(methods_ids),
    var_type="|".join(["SNV", "Indel", "DEL", "CSV", "INV", "INS"])


def get_region_input_from_input_config(wildcards):
    return config["regions"][wildcards.ref][wildcards.region]


rule get_regions:
    input:
        get_region_input_from_input_config
    output:
        dir_work + "benchmarks_file/{ref}.{region}.bed"
    threads: config["threads"]["default"]
    run:
        if f"{input}"[-2:] == "gz":
            shell("{bgzip} -d -c {input} >{output}")
        else:
            shell("cp {input} {output}")


def get_benchmark_vcf_input(wildcards):
    return config["variants"][wildcards.ref][wildcards.var_type]


rule get_vcf_in_this_regions:
    input:
        bed=dir_work + "benchmarks_file/{ref}.{region}.bed",
        vcf=get_benchmark_vcf_input
    output:
        vcf=dir_work + "benchmarks_file/{ref}.{region}.{region}.{var_type}.vcf.gz"
    threads: config["threads"]["default"]
    run:
        ""


def get_benchmark_vcf_query(wildcards):
    return config["variants"][wildcards.ref][wildcards.var_type]


rule get_query_vcf_in_this_regions:
    input:
        bed=dir_work + "benchmarks_file/{ref}.{region}.bed",
        vcf=get_benchmark_vcf_input
    output:
        vcf=dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/{ref}.{sample}.{var_type}.{method}.{region}.vcf.gz"
    threads: config["threads"]["default"]
    run:
        vcf=pysam.VariantFile()
        chrom


rule tabix:
    input:
        "{prefix}.vcf.gz"
    output:
        "{prefix}.vcf.gz.tbi"
    threads: config["threads"]["default"]
    run:
        shell("{tabix} {input}")
