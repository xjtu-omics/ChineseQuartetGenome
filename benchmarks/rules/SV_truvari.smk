# ======================================================================================================================
# Project: ChinseQuartetGenome
# Script : SV_truvari.smk
# Author : Peng Jia
# Date   :  2023/2/23
# Email  : pengjia@stu.xjtu.edu.cn
# Description: Calculate performance of SV detection by truvari
# ======================================================================================================================
import os.path
import json

rule truvari:
    input:
        bed=dir_work + "benchmarks_file/{ref}.{region}.bed",
        bench_vcf=dir_work + "benchmarks_file/{ref}.{region}.{region}.{var_type}.vcf.gz",
        bench_vcf_idx=dir_work + "benchmarks_file/{ref}.{region}.{region}.{var_type}.vcf.gz.tbi",
        query_vcf=dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/{ref}.{sample}.{var_type}.{method}.{region}.vcf.gz",
        query_vcf_idx=dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/{ref}.{sample}.{var_type}.{method}.{region}.vcf.gz.tbi",
    output:
        dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/truvari/{ref}.{sample}.{var_type}.{method}.{region}.ok"
    threads: config["threads"]["truvari"] if "truvari" in config["threads"] else 1
    wildcard_constraints:
        var_type="DEL|INS"
    run:
        out_prefix = f"{output}"[:-3]
        if os.path.exists(f"{out_prefix}"):
            shell("rm -rf {out_prefix}")
        shell("{truvari} bench -b {input.bench_vcf} -c {input.query_vcf} -p 0 --dup-to-ins -o {out_prefix}")
        shell("touch {output}")

rule truvari_metric:
    input:
        dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/truvari/{ref}.{sample}.{var_type}.{method}.{region}.ok"
    output:
        dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/truvari/{ref}.{sample}.{var_type}.{method}.{region}.metric.csv"
    run:
        out_summary = f"{input}"[:-3] + "/summary.txt"
        data_metrics = pd.DataFrame(columns=["ref", "sample", "var_type", "method", "region", "comparison_tool",
                                             "Total_bench", "Total_query", "FP", "FN", "TP",
                                             "Precision", "Recall/Sensitivity", "F1-score"])
        with open(f'{out_summary}','r') as f:
            data = json.load(f)
        data_metrics.loc[
            f"{wildcards.ref}_{wildcards.sample}_{wildcards.var_type}_{wildcards.method}_{wildcards.region}_happy"] = [
            wildcards.ref, wildcards.sample, wildcards.var_type, wildcards.method, wildcards.ref, "happy",
            data["TP-base"] + data["FN"], data["TP-call"] + data["FP"], data["FP"], data["FN"], data["TP-base"],
            data["precision"], data["recall"], data["f1"]
        ]
        data_metrics.to_csv(f"{output}")
