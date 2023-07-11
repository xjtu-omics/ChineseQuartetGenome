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
        bench_vcf=dir_work + "benchmarks_file/{ref}.{region}.{var_type}.vcf.gz",
        bench_vcf_idx=dir_work + "benchmarks_file/{ref}.{region}.{var_type}.vcf.gz.tbi",
        query_vcf=dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/{ref}.{sample}.{var_type}.{method}.{region}.vcf.gz",
        query_vcf_idx=dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/{ref}.{sample}.{var_type}.{method}.{region}.vcf.gz.tbi",
    output:
        dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/happy/{ref}.{sample}.{var_type}.{method}.{region}.ok"
    threads: config["threads"]["happy"] if "happy" in config["threads"] else 1
    wildcard_constraints:
        var_type="Indel|SNV|SNVIndel"
    run:
        ref_fa = ""
        ref_sdf = ""
        if "ref" in config and f"{wildcards.ref}" in config["ref"] and "fa" in config["ref"][wildcards.ref]:
            ref_fa = config["ref"][wildcards.ref]["fa"]

        if len(ref_fa) > 1:
            ref_str_path = f"export HGREF={ref_fa}"
            command_ref = f" -r {ref_fa} "
        else:
            ref_str_path = "echo No ref"
            command_ref = ""
        if "ref" in config and f"{wildcards.ref}" in config["ref"] and "sdf" in config["ref"][wildcards.ref]:
            ref_sdf = config["ref"][wildcards.ref]["sdf"]

        if len(ref_sdf) > 1:
            command_sdf = f" --engine-vcfeval-template {ref_sdf}"
        else:
            command_sdf = ""
        print("ref",command_ref)
        print("sdf",command_sdf)

        hap_eval = os.path.abspath(config["software"]["hap_eval"] if "hap_eval" in config["software"] else "hap.py")
        rtg = os.path.abspath(config["software"]["rtg"] if "hap_eval" in config["software"] else "rtg")
        hap_prefix = hap_eval[:-6]
        out_prefix = f"{output}"[:-3]
        shell("""
        {ref_str_path}
        export PATH={hap_prefix}:$PATH 
        {hap_eval} {input.bench_vcf} {input.query_vcf} -f {input.bed} -V -X -L -D --bcftools-norm --engine vcfeval \
        --engine-vcfeval-path {rtg} --lose-match-distance 100 {command_sdf}\
        --unhappy --threads {threads} -f {input.bed}  -o {out_prefix} {command_ref}
        """)
        shell("touch {output}")

rule happy_metric:
    input:
        tag=dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/happy/{ref}.{sample}.{var_type}.{method}.{region}.ok",
        bench_vcf=dir_work + "benchmarks_file/{ref}.{region}.{var_type}.vcf.gz",
        bench_vcf_idx=dir_work + "benchmarks_file/{ref}.{region}.{var_type}.vcf.gz.tbi",
        query_vcf=dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/{ref}.{sample}.{var_type}.{method}.{region}.vcf.gz",
        query_vcf_idx=dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/{ref}.{sample}.{var_type}.{method}.{region}.vcf.gz.tbi",
        bed=dir_work + "benchmarks_file/{ref}.{region}.bed",
    output:
        dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/happy/{ref}.{sample}.{var_type}.{method}.{region}.metric.csv"
    run:

        out_summary = f"{input.tag}"[:-2] + "summary.csv"
        vcf_path = f"{input.tag}"[:-2] + "vcf.gz"
        tp_vcf_path = f"{input.tag}"[:-2] + "tp.vcf.gz"
        fp_vcf_path = f"{input.tag}"[:-2] + "fp.vcf.gz"
        fn_vcf_path = f"{input.tag}"[:-2] + "fn.vcf.gz"

        vcf_in = pysam.VariantFile(f"{vcf_path}")
        # vcf = pysam.VariantFile(f"{input.vcf}")
        chroms_lens = {i: vcf_in.header.contigs[i].length for i in vcf_in.header.contigs.keys()}
        anno_pots_bn = {i: np.zeros(j) for i, j in chroms_lens.items()}
        for line in open(f"{input.bed}"):
            chrom, start, end = line[:-1].split("\t")[:3]
            if chrom not in anno_pots_bn: continue
            start, end = int(start), int(end)
            anno_pots_bn[chrom][start:end + 1] = 1
        tp_vcf = pysam.VariantFile(f"{tp_vcf_path}","w",header=vcf_in.header)
        fp_vcf = pysam.VariantFile(f"{fp_vcf_path}","w",header=vcf_in.header)
        fn_vcf = pysam.VariantFile(f"{fn_vcf_path}","w",header=vcf_in.header)
        tp = {}
        fp = {}
        fn = {}
        for rec in vcf_in:
            indel_len = abs(len(rec.alts[0]) - len(rec.ref))
            if indel_len >= 50:
                # print(indel_len)
                continue
            if anno_pots_bn[rec.contig][rec.pos:rec.stop + 1].mean() < 1:
                continue
            if "UNK" in [rec.samples["TRUTH"]["BD"], rec.samples["QUERY"]["BD"]]:
                # print("UNK")
                continue

            if ("TP" in [rec.samples["TRUTH"]["BD"], rec.samples["QUERY"]["BD"]]) or \
                    (rec.samples["TRUTH"]["BD"] in ["TP", "FN", "N"] and \
                     rec.samples["QUERY"]["BD"] in ["TP", "FP", "N"]):
                tp[f"{rec.contig}_{rec.pos}_{rec.stop}"] = 1
                tp_vcf.write(rec)
            elif rec.samples["TRUTH"]["BD"] in ["FN"]:
                fn[f"{rec.contig}_{rec.pos}_{rec.stop}"] = 1
                fn_vcf.write(rec)
            elif rec.samples["QUERY"]["BD"] in ["FP"]:
                fp[f"{rec.contig}_{rec.pos}_{rec.stop}"] = 1
                fp_vcf.write(rec)
                # else:
                #     print(rec)

                pass
        tp_vcf.close()
        fp_vcf.close()
        fn_vcf.close()
        vcf_in.close()
        data_metrics = pd.DataFrame(columns=["ref", "sample", "var_type", "method", "region", "comparison_tool",
                                             "Total_bench", "Total_query", "FP", "FN", "TP_bench", "TP_query",
                                             "FP_Border",
                                             "Precision", "Recall/Sensitivity", "F1-score"])
        print(len(tp),len(fp),len(fn))
        tp_num, fp_num, fn_num = len(tp), len(fp), len(fn)

        total_bench = tp_num + fn_num
        total_query = tp_num + fp_num
        precision = tp_num / total_query
        recall = tp_num / total_bench
        if (precision + recall) == 0:
            f1 = 0
        else:
            f1 = 2 * (precision * recall) / (precision + recall)
        data_metrics.loc[
            f"{wildcards.ref}_{wildcards.sample}_{wildcards.var_type}_{wildcards.method}_{wildcards.region}_happy"] = [
            wildcards.ref, wildcards.sample, wildcards.var_type, wildcards.method, wildcards.region, "happy",
            total_bench, total_query, tp_num, fn_num, tp_num, tp_num, ".",
            precision, recall, f1
        ]
        print(data_metrics)
        data_metrics.to_csv(f"{output}")

# continue


# data_metrics = pd.DataFrame(columns=["ref", "sample", "var_type", "method", "region", "comparison_tool",
#                                      "Total_bench", "Total_query", "FP", "FN", "TP",
#                                      "Precision", "Recall/Sensitivity", "F1-score"])
# data = pd.read_csv(f"{out_summary}").loc[0]
# data_metrics.loc[
#     f"{wildcards.ref}_{wildcards.sample}_{wildcards.var_type}_{wildcards.method}_{wildcards.region}_happy"] = [
#     wildcards.ref, wildcards.sample, wildcards.var_type, wildcards.method, wildcards.region, "happy",
#     data["TRUTH.TOTAL"], data["QUERY.TOTAL"], data["QUERY.FP"], data["TRUTH.FN"], data["TRUTH.TP"],
#     data["METRIC.Precision"], data["METRIC.Recall"], data["METRIC.F1_Score"]
# ]
# data_metrics.to_csv(f"{output}")
