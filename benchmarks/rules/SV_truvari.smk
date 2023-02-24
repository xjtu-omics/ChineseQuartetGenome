# ======================================================================================================================
# Project: ChinseQuartetGenome
# Script : SV_truvari.smk TODO check 
# Author : Peng Jia
# Date   :  2023/2/23
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================

rule truvari:
    input:
        bed=dir_work + "benchmarks_file/{ref}.{region}.bed",
        bench_vcf=dir_work + "benchmarks_file/{ref}.{region}.{region}.{var_type}.vcf.gz",
        bench_vcf_idx=dir_work + "benchmarks_file/{ref}.{region}.{region}.{var_type}.vcf.gz",
        query_vcf=dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/{ref}.{sample}.{var_type}.{method}.{region}.vcf.gz",
        query_vcf_idx=dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/{ref}.{sample}.{var_type}.{method}.{region}.vcf.gz.tbi",
    output:
        dir_work + "results/{ref}.{sample}.{var_type}.{method}.{region}/truvari/{ref}.{sample}.{var_type}.{method}.{region}"
    threads: config["threads"]["truvari"] if "truvari" in config["threads"] else 1
    wildcard_constraints:
        var_type="DEL|INS"
    run:
        shell ("{truvari} bench -b {input.bench_vcf} -c {input.query_vcf} -p 0 --dup-to-ins {output}")

