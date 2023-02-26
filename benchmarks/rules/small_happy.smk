# ======================================================================================================================
# Project: ChinseQuartetGenome
# Script : small_happy.smk
# Author : Peng Jia
# Date   :  2023/2/23
# Email  : pengjia@stu.xjtu.edu.cn
# Description: Calculate the performance of Small variant detection with hap.py
# ======================================================================================================================
import os.path

import pandas as pd

rule happy:
    input:
        bed=dir_work + "benchmarks_file/{ref}.{region}.bed",
        bench_vcf=dir_work + "benchmarks_file/{ref}.{region}.{region}.{var_type}.vcf.gz",
        bench_vcf_idx=dir_work + "benchmarks_file/{ref}.{region}.{region}.{var_type}.vcf.gz.tbi",
        query_vcf=dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/{ref}.{sample}.{var_type}.{method}.{region}.vcf.gz",
        query_vcf_idx=dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/{ref}.{sample}.{var_type}.{method}.{region}.vcf.gz.tbi",
    output:
        dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/happy/{ref}.{sample}.{var_type}.{method}.{region}.ok"
    threads: config["threads"]["happy"] if "happy" in config["threads"] else 1
    wildcard_constraints:
        var_type="Indel|SNV"
    run:
        ref_fa = ""
        if "ref" in config and f"{wildcards.ref}" in config["ref"] and "fa" in config["ref"][wildcards.ref]:
            ref_fa = config["ref"][wildcards.ref]["fa"]

        if len(ref_fa) > 1:
            ref_str_path = f"export HGREF={ref_fa}"
            command_ref = f" -r {ref_fa} "
        else:
            ref_str_path = ""
            command_ref = ""

        hap_eval = os.path.abspath(config["software"]["hap_eval"] if "hap_eval" in config["software"] else "hap.py")
        hap_prefix = hap_eval[:-6]
        out_prefix = f"{output}"[:-3]
        shell("""
        {ref_str_path}
        export PATH={hap_prefix}:$PATH 
        {hap_eval} {input.bench_vcf} {input.query_vcf} --threads {threads} -f {input.bed}  -o {out_prefix} {command_ref}
        """)
        shell("touch {output}")

rule happy_metric:
    input:
        dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/happy/{ref}.{sample}.{var_type}.{method}.{region}.ok"
    output:
        dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/happy/{ref}.{sample}.{var_type}.{method}.{region}.metric.csv"
    run:
        out_summary = f"{input}"[:-2] + "summary.csv"
        data_metrics = pd.DataFrame(columns=["ref", "sample", "var_type", "method", "region", "comparison_tool",
                                             "Total_bench", "Total_query", "FP", "FN", "TP",
                                             "Precision", "Recall/Sensitivity", "F1-score"])
        data = pd.read_csv(f"{out_summary}").loc[0]
        data_metrics.loc[
            f"{wildcards.ref}_{wildcards.sample}_{wildcards.var_type}_{wildcards.method}_{wildcards.region}_happy"] = [
            wildcards.ref, wildcards.sample, wildcards.var_type, wildcards.method, wildcards.ref, "happy",
            data["TRUTH.TOTAL"], data["QUERY.TOTAL"], data["QUERY.FP"], data["TRUTH.FN"], data["TRUTH.TP"],
            data["METRIC.Precision"], data["METRIC.Recall"], data["METRIC.F1_Score"]
        ]
        data_metrics.to_csv(f"{output}")
