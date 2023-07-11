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
        bench_vcf=dir_work + "benchmarks_file/{ref}.{region}.{var_type}.vcf.gz",
        bench_vcf_idx=dir_work + "benchmarks_file/{ref}.{region}.{var_type}.vcf.gz.tbi",
        query_vcf=dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/{ref}.{sample}.{var_type}.{method}.{region}.vcf.gz",
        query_vcf_idx=dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/{ref}.{sample}.{var_type}.{method}.{region}.vcf.gz.tbi",
    output:
        dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/truvari/{ref}.{sample}.{var_type}.{method}.{region}.ok"
    threads: config["threads"]["truvari"] if "truvari" in config["threads"] else 1
    wildcard_constraints:
        var_type="DEL|INS"
    run:
        SV_bench_minimum_size = max(int(paras["SV_bench_minimum_size"]) - 20,0)
        SV_query_minimum_size = max(int(paras["SV_query_minimum_size"]) - 20,0)
        SV_sizemax = int(paras["SV_sizemax"])
        out_prefix = f"{output}"[:-3]
        if os.path.exists(f"{out_prefix}"):
            shell("rm -rf {out_prefix}")
        shell("{truvari} bench -b {input.bench_vcf} -c {input.query_vcf} -p 0 -P 0.1  -S {SV_bench_minimum_size} "
              "-s {SV_query_minimum_size}  --sizemax {SV_sizemax}  --multimatch --dup-to-ins -o {out_prefix}")
        shell("touch {output}")

rule truvari_metric:
    input:
        dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/truvari/{ref}.{sample}.{var_type}.{method}.{region}.ok"
    output:
        dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/truvari/{ref}.{sample}.{var_type}.{method}.{region}.metric.csv"
    run:
        import pysam

        size_tolerate_region_boundary = int(paras["size_tolerate_region_boundary_SV"])
        out_summary = f"{input}"[:-3] + "/summary.txt"
        out_prefix = f"{input}"[:-3] + "/"
        data_metrics = pd.DataFrame(columns=["ref", "sample", "var_type", "method", "region", "comparison_tool",
                                             "Total_bench", "Total_query", "FP", "FN", "TP_bench", "TP_query",
                                             "FP_Border",
                                             "Precision", "Recall/Sensitivity", "F1-score"])
        SV_bench_minimum_size = max(int(paras["SV_bench_minimum_size"]),0)
        SV_query_minimum_size = max(int(paras["SV_query_minimum_size"]),0)
        min_size = max(SV_query_minimum_size,SV_bench_minimum_size)
        SV_sizemax = int(paras["SV_sizemax"])
        # fp_num, fn_num, tp_bench, tp_query = [0] * 4
        rec_num = {i: 0 for i in ["FP", "FN", "TP_bench", "TP_query"]}
        fp_border = 0
        for name, file_name in zip(["FP", "FN", "TP_bench", "TP_query"],["fp.vcf", "fn.vcf", "tp-base.vcf",
                                                                         "tp-call.vcf"]):
            vcf = pysam.VariantFile(f"{out_prefix}{file_name}")
            vcf_new = pysam.VariantFile(f"{out_prefix}{name}_new.vcf","w",header=vcf.header)
            idxs = set()
            for rec in vcf.fetch():
                this_size = abs(rec.info["SVLEN"])
                if name == "FP" and rec.info[f"Border.{size_tolerate_region_boundary}"] == "Yes":
                    fp_border += 1
                    continue
                if min_size <= this_size <= SV_sizemax:
                    idx = f"{rec.contig}_{rec.pos}_{rec.stop}"
                    if idx in idxs:
                        continue
                    idxs.add(idx)
                    rec_num[name] += 1
                    vcf_new.write(rec)
            vcf.close()
            vcf_new.close()
        total_query = rec_num["TP_query"] + rec_num["FP"]
        total_bench = rec_num["TP_bench"] + rec_num["FN"]
        precision = rec_num["TP_query"] / total_query
        recall = rec_num["TP_bench"] / total_bench
        if (precision + recall) == 0:
            f1 = 0
        else:
            f1 = 2 * (precision * recall) / (precision + recall)
        data_metrics.loc[
            f"{wildcards.ref}_{wildcards.sample}_{wildcards.var_type}_{wildcards.method}_{wildcards.region}_truvari"] = [
            wildcards.ref, wildcards.sample, wildcards.var_type, wildcards.method, wildcards.region, "truvari",
            total_bench, total_query, rec_num["FP"], rec_num["FN"], rec_num["TP_bench"], rec_num["TP_query"], fp_border,
            precision, recall, f1
        ]
        data_metrics.to_csv(f"{output}")
