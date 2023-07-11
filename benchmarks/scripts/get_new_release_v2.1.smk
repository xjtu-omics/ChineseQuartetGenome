# ======================================================================================================================
# Project: Project_Chinese_Quartet
# Script : get_new_release_v2.1
# Author : Peng Jia
# Date   :  2023/2/5
# Email  : pengjia@stu.xjtu.edu.cn
# Description: xxx
# ======================================================================================================================
import pysam
import numpy as np
import pandas as pd
import gzip

tabix = "/data/home/pengjia/miniconda3/envs/default/bin/tabix"

dir_work_old = "/data/home/pengjia/Project_Human_DATA/reNBT/"
dir_work = "/data/home/pengjia/Project_Human_DATA/reNBT2/"
dir_release_fa = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/assm_res/05_post_assembly/rename/"
hap1 = dir_release_fa + "ChineseQuartet.hap1.hifiasm_latest.rename.fa"
hap2 = dir_release_fa + "ChineseQuartet.hap2.hifiasm_latest.rename.fa"
path_hg38 = "/data/DATA/Reference/human/GRCh38.p13/GRCh38.p13.genome.fa"

path_t2t = "/data/DATA/Reference/human/T2T_v2.0/chm13v2.0.fa"
my_haps = {"pat": hap1, "mat": hap2, }
ref_haps = {"GRCh38": path_hg38, "T2T": path_t2t}
var_types = ["INS", "DEL", "SNV", "Indel", "CSV", "INV"]
seqtk = "/data/home/pengjia/miniconda3/envs/assm/bin/seqtk"
bedtools = "/data/home/pengjia/miniconda3/envs/default/bin/bedtools"
bcftools = "/data/home/pengjia/miniconda3/envs/default/bin/bcftools"
dir_jdlin = "/data/home/pengjia/REF/GRCh38_database/repeat/fromJdLin/"
beds = {}
# beds["RM"] = dir_jdlin + "rmsk.bed.gz"
beds["SD"] = dir_jdlin + "seg_dup.bed.gz"
beds["SR"] = dir_jdlin + "simplerepeat.bed.gz"
beds["STR"] = "/data/home/pengjia/REF/GRCh38_database/repeat/GRCh38_full_analysis_set_plus_decoy_hla.fa.str.bed.gz"
beds["VNTR"] = "/data/home/pengjia/REF/table_browser/GRCh38_VNTR.sort.bed.gz"

wildcard_constraints:
    var_type="|".join(var_types)

rule all:
    input:
        expand(dir_work + "02_variants/release/ChineseQuartet.twin.germline.GRCh38.{var_type}.freeze.v2.1.vcf.gz.tbi",
            var_type=["DEL", "INS", "SNV", "Indel","SNVIndel"]),# "DEL", "INS",



def get_vcf(wildcards):
    # dir_release_vcf = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/ref_based_analysis/comparation2/detail/final"
    dir_release_vcf = "/data/home/pengjia/Project_Human_DATA/reNBT2/02_variants/merged_final"
    # if wildcards.ref_hap in ["GRCh38"]:
    return f"{dir_release_vcf}/{wildcards.var_type}/ChineseQuartet.merged3tech.{wildcards.var_type}.sort.pass.val.RP.DPMQ.final.vcf.gz"


# return f"{dir_release_vcf}/{wildcards.var_type}/ChineseQuartet.merged3tech.{wildcards.var_type}.sort.pass.val.RP.vcf.gz"


def get_ref_fai(wildcards):
    return ref_haps[f"GRCh38"] + ".fai"


rule get_new_vcf:
    input:
        vcf=get_vcf,
        bed=dir_work + "02_variants/release/ChineseQuartet.twin.germline.GRCh38.benchmark_region.v2.1.bed",
        fai=get_ref_fai
    output:
        dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.vcf.gz"
    wildcard_constraints: var_rype="|".join(var_types)
    threads: 4
    run:
        chroms = [f"chr{i}" for i in range(1,23)] + ["chrX"] + ["chrY"]
        chroms_lens = {}
        for item in open(f"{input.fai}"):
            chrom, chrom_len = item[:-1].split("\t")[:2]
            chroms_lens[chrom] = int(chrom_len)
        anno_pots = {i: np.zeros(chroms_lens[i]) for i in chroms}

        for line in open(f"{input.bed}"):
            lineinfo = line[:-1].split("\t")
            chrom, start, end = lineinfo[:3]
            if chrom not in chroms: continue
            start, end = int(start), int(end)
            anno_pots[chrom][start:end] += 1
        total_num = 0
        check_num = 0
        vcf = pysam.VariantFile(f"{input.vcf}",threads=int(f"{threads}"))
        my_header = vcf.header
        my_header.add_line(f'##INFO=<ID=Benchmarks_Ratio,Number=1,Type=Float,Description="Ratio of this variant in benchmark regions.">')
        vcf_new = pysam.VariantFile(f"{output}",header=my_header,mode="w",threads=int(f"{threads}"))
        for rec in vcf.fetch():
            total_num += 1
            chrom = rec.contig
            start = rec.pos
            # svlen = abs(rec.info["SVLEN"])
            # svlen_new = 1000 if svlen < 1000 else svlen
            end = rec.stop + 1
            ratio = np.mean(anno_pots[chrom][start:end])
            my_rec = f"{rec}"[:-1]
            rec.info["Benchmarks_Ratio"] = ratio

            if wildcards.var_type in ["SNV", "Indel"]:
                if rec.info["Benchmarks_Ratio"] < 1:
                    continue
            elif wildcards.var_type in ["DEL", "INS", "INV"]:
                if rec.info["Benchmarks_Ratio"] < 0.8:
                    continue

            vcf_new.write(rec)

rule rephased:
    input:
        vcf=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.vcf.gz",
        bed=dir_work_old + "01_assembly/read_phasing/hap_regions/ChineseQuartet.ref.GRCh38.cov.bed",
        fai=get_ref_fai,
    output:
        dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.rephased.vcf.gz"
    wildcard_constraints: var_rype="|".join(var_types)
    threads: 4
    run:
        chroms = [f"chr{i}" for i in range(1,23)] + ["chrX"] + ["chrY"]
        chroms_lens = {}
        for item in open(f"{input.fai}"):
            chrom, chrom_len = item[:-1].split("\t")[:2]
            chroms_lens[chrom] = int(chrom_len)
        anno_pots = {i: np.zeros(chroms_lens[i]) for i in chroms}

        for line in open(f"{input.bed}"):
            lineinfo = line[:-1].split("\t")
            chrom, start, end = lineinfo[:3]
            if chrom not in chroms: continue
            start, end = int(start), int(end)
            anno_pots[chrom][start:end] += 1
        total_num = 0
        check_num = 0
        vcf = pysam.VariantFile(f"{input.vcf}",threads=int(f"{threads}"))
        my_header = vcf.header
        my_header.add_line(f'##INFO=<ID=Phased_Ratio,Number=1,Type=Float,Description="Ratio of this variant in phased regions.">')
        vcf_new = pysam.VariantFile(f"{output}",header=my_header,mode="w",threads=int(f"{threads}"))
        sample_names = vcf.header.samples
        for rec in vcf.fetch():
            total_num += 1
            chrom = rec.contig
            start = rec.pos
            # svlen = abs(rec.info["SVLEN"])
            # svlen_new = 1000 if svlen < 1000 else svlen
            end = rec.stop + 1
            ratio = np.mean(anno_pots[chrom][start:end])
            my_rec = f"{rec}"[:-1]
            rec.info["Phased_Ratio"] = ratio
            if wildcards.var_type in ["SNV", "Indel"]:
                if ratio < 1:
                    for sample_name in sample_names:
                        # print(rec.samples[sample_name].phased)
                        rec.samples[sample_name].phased = False
            elif wildcards.var_type in ["DEL", "INS", "INV"]:
                if ratio < 0.5:
                    for sample_name in sample_names:
                        # print(rec.samples[sample_name].phased)
                        rec.samples[sample_name].phased = False
            vcf_new.write(rec)

#
rule extract_SV_for_Small:
    input:
        get_vcf
    output:
        dir_work + "02_variants/confidence/Filter_regions/ChineseQuartet.twin.germline.GRCh38.{var_type}.filter_forsmall.bed"
    wildcard_constraints:
        var_type="DEL|INS"
    run:
        border_size = 0
        vcf = pysam.VariantFile(F"{input}")
        out_bed = open(f"{output}","w")
        for rec in vcf.fetch():
            if abs(rec.info["SVLEN"]) > 80:
                out_bed.write(f"{rec.contig}\t{rec.pos - border_size}\t{rec.stop + border_size}\n")
        out_bed.close()


rule merge_sv_regions:
    input:
        expand(dir_work + "02_variants/confidence/Filter_regions/ChineseQuartet.twin.germline.GRCh38.{var_type}.filter_forsmall.bed",
            var_type=["DEL", "INS", "CSV"])
    output:
        dir_work + "02_variants/confidence/Filter_regions/ChineseQuartet.twin.SV.filter.bed"
    run:
        shell("""cat {input} | cut -f1,2,3 |awk 'BEGIN {{OFS="\t"}}{{ print $1,$2-50,$3+50}}' | awk 'BEGIN {{OFS="\t"}}{{if($2<0)print $1,0,$3; else print $0}}' |{bedtools} sort |{bedtools} merge |{bedtools} sort >{output}""")

rule get_first_regions:
    input:
        # bed=dir_work + "02_variants/confidence/ChineseQuartet.twin.GRCh38.benchmark_region.diploid.bed",
        bed=dir_work + "02_variants/release/ChineseQuartet.twin.germline.GRCh38.SV_Tier1_draft.v2.1.bed",
        SV=dir_work + "02_variants/confidence/Filter_regions/ChineseQuartet.twin.SV.filter.bed",
        filter=dir_work + "02_variants/confidence/Filter_regions/GRCh38_filter.bed",
        gap=dir_work + "02_variants/confidence/Filter_regions/GRCh38_gap.bed",
    # sd=dir_work + "02_variants/confidence/Filter_regions/GRCh38.SD.bed",
    output:
        dir_work + "02_variants/confidence/ChineseQuartet.twin.GRCh38.small_high_confidence_region.first.bed"
    run:
        shell("cat {input.SV} {input.filter} {input.gap}|cut -f 1,2,3|{bedtools} sort > {output}_tmp")
        shell("{bedtools} subtract -a {input.bed} -b {output}_tmp |{bedtools} sort >{output}")
        shell("rm {output}_tmp")

rule low_confidence_Small:
    input:
        vcf=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.rephased.vcf.gz",
        vcf_idx=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.rephased.vcf.gz.tbi",
        # sv=dir_work + "02_variants/confidence/Filter_regions/ChineseQuartet.twin.SV.filter.bed",
        repeat=dir_work + "02_variants/confidence/Filter_regions/ChineseQuartet.twin.repeat.filter_forSmall2.bed",
        fai=get_ref_fai,
        bed=dir_work + "02_variants/confidence/ChineseQuartet.twin.GRCh38.small_high_confidence_region.first.bed"
    output:
        vcf=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.draft.vcf.gz",
        high_conf_bed=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.rephased.high_confidence.bed",
        low_conf_bed=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.rephased.low_confidence.bed",
    threads: 4
    wildcard_constraints:
        var_type="SNV|Indel"
    run:
        border_size = 200
        chroms = [f"chr{i}" for i in range(1,23)] + ["chrX"] + ["chrY"]
        chroms_lens = {}
        for item in open(f"{input.fai}"):
            chrom, chrom_len = item[:-1].split("\t")[:2]
            chroms_lens[chrom] = int(chrom_len)
        anno_pots_bn = {i: np.zeros(chroms_lens[i]) for i in chroms}
        # anno_pots_sv = {i: np.zeros(chroms_lens[i]) for i in chroms}
        anno_pots_repeat = {i: np.zeros(chroms_lens[i]) for i in chroms}
        for line in open(f"{input.bed}"):
            chrom, start, end = line[:-1].split("\t")[:3]
            start, end = int(start), int(end)
            anno_pots_bn[chrom][start:end + 1] = 1
        # for line in open(f"{input.sv}"):
        #     chrom, start, end = line[:-1].split("\t")[:3]
        #     start, end = int(start), int(end)
        #     anno_pots_sv[chrom][start:end + 1] = 1
        for line in open(f"{input.repeat}"):
            chrom, start, end = line[:-1].split("\t")[:3]
            if chrom not in chroms: continue
            start, end = int(start), int(end)
            anno_pots_repeat[chrom][start:end + 1] = 1
        vcf = pysam.VariantFile(f"{input.vcf}",threads=int(f"{threads}"))
        my_header = vcf.header
        my_header.add_line(f'##INFO=<ID=Qual_ANN,Number=1,Type=String,Description="High_confidence.">')
        my_header.add_line(f'##INFO=<ID=Tier,Number=1,Type=String,Description="Tier of variants, the Tier1 variants have more evidences than Tier2.">')
        my_header.add_line(f'##INFO=<ID=Border_Tier1,Number=1,Type=String,Description="Is this mutation occurring at the Tier1 region boundary?.">')
        my_header.add_line(f'##INFO=<ID=Indel_length,Number=1,Type=Integer,Description="Indel length">')

        vcf_new = pysam.VariantFile(f"{output.vcf}",header=my_header,mode="w",threads=int(f"{threads}"))
        high_conf_bed = open(f"{output.high_conf_bed}","w")
        low_conf_bed = open(f"{output.low_conf_bed}","w")
        num = 0
        total_num = 0
        num_bn = 0
        num_bn_supp2 = 0
        for rec in vcf.fetch():
            if abs(len(rec.alts[0]) - len(rec.ref)) > 80: continue
            total_num += 1
            if anno_pots_bn[rec.chrom][rec.pos:rec.stop + 1].mean() < 1:
                rec.info["Qual_ANN"] = "Region_Filter"
            else:
                if rec.info["SUPP"] == 3:
                    rec.info["Qual_ANN"] = "High_confidence"
                elif rec.info["SUPP"] == 2:
                    if anno_pots_repeat[rec.chrom][rec.pos:rec.stop + 1].sum() > 1:
                        rec.info["Qual_ANN"] = "Small_variant_in_repeat_region"
                    elif rec.info["Qual"] != "PASS":
                        rec.info["Qual_ANN"] = "NOPASS"
                    else:
                        rec.info["Qual_ANN"] = "High_confidence"
                else:
                    # if anno_pots_repeat[rec.chrom][rec.pos:rec.stop + 1].sum() > 1:
                    #     rec.info["Qual_ANN"] = "Small_variant_in_repeat_region"
                    # elif rec.info["Qual"] == "PASS":
                    #     rec.info["Qual_ANN"] = "High_confidence"
                    # else:
                    rec.info["Qual_ANN"] = "Tech_Specific"
            #
            # else:
            #     rec.info["Qual_ANN"] = "Tech_Specific"
            #
            # rec.info["Qual_ANN"] = ""
            #
            # if rec.info["Qual"] != "PASS":
            #     rec.info["Qual_ANN"] = "NOPASS"
            # elif "Yes" in [rec.info[i] for i in ["R_SD", "R_SR", "R_STR", "R_VNTR"]]:
            #     rec.info["Qual_ANN"] = "Repeat_region"
            # elif rec.info["Mappability"]>1:
            #     rec.info["Qual_ANN"] = "LowMapQ"
            # else:
            #     rec.info["Qual_ANN"] = "High_confidence"


            rec.info["Tier"] = "Tier1" if "High_confidence" in rec.info["Qual_ANN"] else "Tier2"
            vcf_new.write(rec)
            if rec.info["Benchmarks_Ratio"] < 1: continue
            # if anno_pots_bn[rec.chrom][rec.pos:rec.stop + 1].mean() < 1: continue
            if rec.info["Qual_ANN"] in ["High_confidence"]:
                high_conf_bed.write(f"{rec.chrom}\t{rec.pos - border_size}\t{rec.stop + border_size}\n")
                num_bn += 1
            else:
                low_conf_bed.write(f"{rec.chrom}\t{rec.pos - border_size}\t{rec.stop + border_size}\n")
        # print(wildcards.var_type,num,total_num,num / total_num,total_num - num)
        # print(num_bn,total_num,num_bn_supp2)
        vcf_new.close()
        high_conf_bed.close()
        low_conf_bed.close()


rule low_confidence_SV:
    input:
        vcf=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.rephased.vcf.gz",
        vcf_idx=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.rephased.vcf.gz.tbi"
    output:
        vcf=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.draft.vcf.gz",
        high_conf_bed=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.rephased.high_confidence.bed",
        low_conf_bed=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.rephased.low_confidence.bed",
    wildcard_constraints:
        var_type="DEL|INS"
    threads: 4
    run:
        board_region = 1000
        vcf = pysam.VariantFile(f"{input.vcf}",threads=int(f"{threads}"))
        my_header = vcf.header
        my_header.add_line(f'##INFO=<ID=Qual_ANN,Number=1,Type=String,Description="Quality annotation.">')
        my_header.add_line(f'##INFO=<ID=Tier,Number=1,Type=String,Description="Tier of variants, the Tier1 variants have more evidences than Tier2.">')
        my_header.add_line(f'##INFO=<ID=Border_Tier1,Number=1,Type=String,Description="Is this mutation occurring at the Tier1 region boundary?.">')

        vcf_new = pysam.VariantFile(f"{output.vcf}",header=my_header,mode="w",threads=int(f"{threads}"))
        high_conf_bed = open(f"{output.high_conf_bed}","w")
        low_conf_bed = open(f"{output.low_conf_bed}","w")
        num = 0
        total_num = 0
        for rec in vcf.fetch():
            total_num += 1
            if rec.info["Qual"] in ["PASS"]:
                if rec.info["SUPP"] == 1:
                    if "Yes" in [rec.info["R_SD"], rec.info["R_SR"], rec.info["R_VNTR"]]:
                        rec.info["Qual_ANN"] = "Call_in_SD_SR_VNTR"
                        num += 1
                    elif "ILM" in rec.info["SUPP_Set"]:
                        rec.info["Qual_ANN"] = "ILM_specific_not_val"
                        num += 1
                    else:
                        rec.info["Qual_ANN"] = "High_confidence"
                else:
                    rec.info["Qual_ANN"] = "High_confidence"

            else:
                rec.info["Qual_ANN"] = "Tech_Specific_call"
                num += 1
            rec.info["Tier"] = "Tier1" if "High_confidence" in rec.info["Qual_ANN"] else "Tier2"
            vcf_new.write(rec)
            if rec.info["Benchmarks_Ratio"] < 0.8: continue
            if rec.info["Qual_ANN"] in ["High_confidence"]:
                high_conf_bed.write(f"{rec.chrom}\t{rec.pos - board_region}\t{rec.stop + board_region}\n")
            else:
                low_conf_bed.write(f"{rec.chrom}\t{rec.pos - board_region}\t{rec.stop + board_region}\n")
        vcf_new.close()
        high_conf_bed.close()
        low_conf_bed.close()


rule process_simple_repeat:
    input:
        lambda wildcards: beds[wildcards.rp_type]
    output:
        dir_work + "02_variants/confidence/Filter_regions/GRCh38_{rp_type}_forSV.bed"
    run:
        # shell("""less {input} |awk 'BEGIN {{OFS="\t"}}{{ print $1,$2,$3}}' |head""")
        # shell("""less {input} |awk '{{ print $1,$2-1000,$3+1000}}' |head """)
        shell("""less {input}|awk 'BEGIN {{OFS="\t"}}{{ print $1,$2-200,$3+200}}' | awk 'BEGIN {{OFS="\t"}}{{if($2<0)print $1,0,$3; else print $0}}' |{bedtools} sort| \
        {bedtools} merge   |{bedtools} sort  >{output}""")


rule process_simple_repeat_forSmall:
    input:
        lambda wildcards: beds[wildcards.rp_type]
    output:
        dir_work + "02_variants/confidence/Filter_regions/GRCh38_{rp_type}_forSmall.bed"
    run:
        # shell("""less {input} |awk 'BEGIN {{OFS="\t"}}{{ print $1,$2,$3}}' |head""")
        # shell("""less {input} |awk '{{ print $1,$2-1000,$3+1000}}' |head """)
        shell("""less {input}|awk 'BEGIN {{OFS="\t"}}{{ print $1,$2-5,$3+5}}' | awk 'BEGIN {{OFS="\t"}}{{if($2<0)print $1,0,$3; else print $0}}' |{bedtools} sort| \
        {bedtools} merge   |{bedtools} sort  >{output}""")


rule process_filter_regions:
    input:
        filter="/data/home/pengjia/REF/GRCh38_database/filter/filter.sort.merge.bed.gz",
    output:
        dir_work + "02_variants/confidence/Filter_regions/GRCh38_filter.bed"
    run:
        shell("less {input} >{output}")

rule extract_gap_bed:
    input:
        ref_haps[f"GRCh38"]
    output:
        dir_work + "02_variants/confidence/Filter_regions/GRCh38_gap.bed",
    run:
        shell("{seqtk} cutN -n 1 -g {input} >{output}")

rule get_SV_low_regions_with_repeat:
    input:
        bed_del=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.DEL.rephased.low_confidence.bed",
        bed_ins=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.INS.rephased.low_confidence.bed",
        # bed_indel=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.Indel.rephased.low_confidence_forSV.bed",
        bed_sr=dir_work + "02_variants/confidence/Filter_regions/GRCh38_SR_forSV.bed",
        bed_vntr=dir_work + "02_variants/confidence/Filter_regions/GRCh38_VNTR_forSV.bed",
        filter=dir_work + "02_variants/confidence/Filter_regions/GRCh38_filter.bed",
        gap=dir_work + "02_variants/confidence/Filter_regions/GRCh38_gap.bed",

    output:
        bed_sv=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.SV.rephased.low_confidence_rp.bed",
    run:
        # shell("cat ")
        shell("cat {input.bed_sr} {input.bed_vntr}| {bedtools} sort | {bedtools} merge |{bedtools} sort  >{output}_tmp2")
        shell("{bedtools} intersect -a {output}_tmp2 -b {input.bed_del} {input.bed_ins} -wa |{bedtools} sort|uniq >{output}_tmp.bed")
        shell("cat {input.bed_del} {input.bed_ins}  {output}_tmp.bed {input.gap} |cut -f 1,2,3 |{bedtools} sort |{bedtools} merge |{bedtools} sort > {output}")
        shell("rm {output}_tmp2 {output}_tmp.bed")

rule get_small_low_qual_regions_with_repeat:
    input:
        bed_snv=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.SNV.rephased.low_confidence.bed",
        bed_indel=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.Indel.rephased.low_confidence.bed",
        bed_sr=dir_work + "02_variants/confidence/Filter_regions/ChineseQuartet.twin.repeat.filter_forSmall.bed"
    output:
        dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.Small.rephased.low_confidence_rp.bed",
    run:
        shell("{bedtools} intersect -a {input.bed_sr} -b {input.bed_snv} {input.bed_indel} -wa |{bedtools} sort|uniq >{output}_tmp.bed")
        shell("cat {input.bed_indel} {input.bed_snv}  {output}_tmp.bed  |cut -f 1,2,3 |{bedtools} sort |{bedtools} merge |{bedtools} sort > {output}")


# shell("rm {output}_tmp.bed")

# shell("cat {input.bed} {input.del_high_conf_bed} {input.ins_high_conf_bed} {input.csv_bed}  |cut -f1,2,3| sort|{bedtools} sort | {bedtools} merge |{bedtools} sort >{output}_tmp ")
# shell("{bedtools} subtract -a {output}_tmp -b {input.del_low_conf_bed} {input.ins_low_conf_bed}  | {bedtools} sort  > {output} ")
#
# # shell("cat {output}_tmp {input.del_high_conf_bed} {input.ins_high_conf_bed} {input.csv_bed}  |cut -f1,2,3| sort|{bedtools} sort | {bedtools} merge |{bedtools} sort >{output} ")
# shell("rm {output}_tmp")


rule only_keep_diploid:
    input:
        base=dir_work + "02_variants/release/ChineseQuartet.twin.germline.{ref_hap}.benchmark_region.v2.1.bed",
        # diploid=dir_work + "03_benchmark_regions/benchmark_regions_diploid/{ref_hap}.uniq.bed"
        diploid=dir_work_old + "02_variants/benchmark_regions_diploid/assm/twin_v1.0.GRCh38.uniq.bam.bed"
    output:
        dir_work + "02_variants/confidence/ChineseQuartet.{sample}.{ref_hap}.benchmark_region.diploid.bed",
    run:
        # shell("cp {input.base} {output}")

        shell("{bedtools} intersect -a {input.base} -b {input.diploid} |{bedtools} sort |awk '$3-$2>1000'> {output}")

rule get_sv_confidence_regions_release_small:
    input:
        snv_high_conf_bed=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.SNV.rephased.high_confidence.bed",
        indel_high_conf_bed=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.Indel.rephased.high_confidence.bed",
        low_conf_bed=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.Small.rephased.low_confidence_rp.bed",
        filter="/data/home/pengjia/REF/GRCh38_database/filter/filter.sort.merge.bed.gz",
        bed=dir_work + "02_variants/confidence/ChineseQuartet.twin.GRCh38.small_high_confidence_region.first.bed"
    output:
        dir_work + "02_variants/confidence/ChineseQuartet.{sample}.{ref_hap}.Small_high_confidence_region.release.bed",
    run:
        shell("cat {input.bed} |cut -f1,2,3| sort|{bedtools} sort | {bedtools} merge |{bedtools} sort >{output}_tmp ")
        shell("{bedtools} subtract -a {output}_tmp -b {input.low_conf_bed}   | {bedtools} sort  > {output} ")


rule get_sv_confidence_regions_release:
    input:
        del_high_conf_bed=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.DEL.rephased.high_confidence.bed",
        ins_high_conf_bed=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.INS.rephased.high_confidence.bed",
        sv_low_conf_bed=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.SV.rephased.low_confidence_rp.bed",
        filter="/data/home/pengjia/REF/GRCh38_database/filter/filter.sort.merge.bed.gz",
        csv_bed=dir_work_old + "02_variants/benchmark_regions/SV_regions/{sample}.{ref_hap}.CSV.extend.csv",
        bed=dir_work + "02_variants/confidence/ChineseQuartet.{sample}.{ref_hap}.benchmark_region.diploid.bed",
    output:
        dir_work + "02_variants/confidence/ChineseQuartet.{sample}.{ref_hap}.SV_high_confidence_region.release.bed",
    # wildcard_constraints:

    run:
        shell("cat {input.bed} {input.del_high_conf_bed} {input.ins_high_conf_bed} {input.csv_bed}  |cut -f1,2,3| sort|{bedtools} sort | {bedtools} merge |{bedtools} sort >{output}_tmp ")
        shell("{bedtools} subtract -a {output}_tmp -b {input.sv_low_conf_bed}   | {bedtools} sort  > {output} ")


# shell("cat {output}_tmp {input.del_high_conf_bed} {input.ins_high_conf_bed} {input.csv_bed}  |cut -f1,2,3| sort|{bedtools} sort | {bedtools} merge |{bedtools} sort >{output} ")
# shell("rm {output}_tmp")


# rule get_small_confidence_regions_first:
#     input:
#         sv=dir_work + "02_variants/confidence/Filter_regions/ChineseQuartet.twin.SV.filter.bed",
#         bed=dir_work + "02_variants/benchmark_regions/{sample}.{ref_hap}.benchmark_region.release.bed",
#         repeat=dir_work + "02_variants/confidence/Filter_regions/ChineseQuartet.twin.repeat.filter.bed"
#     output:
#         dir_work + "02_variants/confidence/{sample}.{ref_hap}.small_high_confidence_region.first.bed"
#     run:
#         shell("{bedtools} subtract -a {input.bed} -b {input.sv} {input.repeat}  | {bedtools} sort  > {output} ")
#
#
def get_conf(wildcards):
    if wildcards.var_type in ["DEL", "INS"]:
        return dir_work + f"02_variants/release/ChineseQuartet.twin.germline.{wildcards.ref_hap}.SV_Tier1.v2.1.bed"
    else:
        return dir_work + f"02_variants/release/ChineseQuartet.twin.germline.{wildcards.ref_hap}.small_Tier1.v2.0.bed"


rule get_final_benchamrk_regions:
    input:
        dir_work_old + "02_variants/benchmark_regions/twin.{ref_hap}.benchmark_region.release.bed",
    output:
        dir_work + "02_variants/release/ChineseQuartet.twin.germline.{ref_hap}.benchmark_region.v2.1.bed",
    run:
        shell("cp {input} {output}")

rule get_final_benchamrk_regions_tier1:
    input:
        tier=dir_work + "02_variants/confidence/ChineseQuartet.twin.{ref_hap}.{size}_high_confidence_region.release.bed",
        base=dir_work + "02_variants/release/ChineseQuartet.twin.germline.{ref_hap}.benchmark_region.v2.1.bed",
    output:
        dir_work + "02_variants/release/ChineseQuartet.twin.germline.{ref_hap}.{size}_Tier1_draft.v2.1.bed"
    wildcard_constraints:
        site="SV|Small"
    run:
        # shell("cp {input.tier} {output}")
        if "SV" in f"{wildcards.size}":
            shell("{bedtools} intersect -a {input.base} -b {input.tier} |bedtools sort | awk '$3-$2>1000'> {output}")
        else:
            shell("{bedtools} intersect -a {input.base} -b {input.tier} |bedtools sort | awk '$3-$2>1000'> {output}")


#
# rule get_final_benchamrk_regions_tier1:
#     input:
#         tier=dir_work + "02_variants/confidence/ChineseQuartet.twin.{ref_hap}.{size}_high_confidence_region.release.bed",
#         base=dir_work + "02_variants/release/ChineseQuartet.twin.germline.{ref_hap}.benchmark_region.v2.1.bed",
#     output:
#         dir_work + "02_variants/release/ChineseQuartet.twin.germline.{ref_hap}.{size}_Tier1_draft.v2.1.bed"
#     wildcard_constraints:
#         site="SV"
#     run:
#         # shell("cp {input.tier} {output}")
#
#         shell("{bedtools} intersect -a {input.base} -b {input.tier} |bedtools sort | awk '$3-$2>1000'> {output}")


rule check_small:
    input:
        bed=dir_work + "02_variants/release/ChineseQuartet.twin.germline.{ref_hap}.Small_Tier1_draft.v2.1.bed",
        vcf=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.draft.vcf.gz",
        vcf_idx=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.draft.vcf.gz.tbi",
        fai=get_ref_fai,
    output:
        vcf=dir_work + "02_variants/release/ChineseQuartet.twin.germline.{ref_hap}.{var_type}.freeze.v2.1.vcf.gz",
    threads: 4
    wildcard_constraints:
        var_type="Indel|SNV"
    run:
        safe_border = 5
        vcf = pysam.VariantFile(f"{input.vcf}")
        new_header_records = [r for r in vcf.header.records if
                              not (r.key == 'bcftools_concatCommand')]
        new_header = pysam.VariantHeader()
        for record in new_header_records:
            new_header.add_line(str(record))
        for i in vcf.header.samples:
            new_header.add_sample(i)
        vcf_new = pysam.VariantFile(f"{output.vcf}",header=new_header,mode="w")

        chroms = [f"chr{i}" for i in range(1,23)] + ["chrX"]
        chroms_lens = {}
        for item in open(f"{input.fai}"):
            chrom, chrom_len = item[:-1].split("\t")[:2]
            chroms_lens[chrom] = int(chrom_len)
        anno_pots_bn = {i: np.zeros(chroms_lens[i]) for i in chroms}
        for line in open(f"{input.bed}"):
            chrom, start, end = line[:-1].split("\t")[:3]
            start, end = int(start), int(end)
            anno_pots_bn[chrom][start:end + 1] = 1
        num = 0
        tier1 = 0
        tier2 = 0
        tier2_a = 0
        tier2_regions = []
        for rec in vcf.fetch():
            # rec.info[]=rec.info["SVLEN"])
            num += 1
            rec.info["Indel_length"] = len(rec.alts[0]) - len(rec.ref)
            rec.info[f"Border_Tier1"] = "No"
            # if rec.info["Benchmarks_Ratio"] < 1:
            #     print(rec.info["Benchmarks_Ratio"])
            if rec.info["Tier"] in "Tier1":
                if anno_pots_bn[rec.contig][rec.pos:rec.stop].mean() < 1:
                    rec.info["Tier"] = "Tier2"
                    tier2_a += 1
                    if anno_pots_bn[rec.contig][rec.pos - safe_border:rec.stop + safe_border].mean() > 0:
                        # print(anno_pots_bn[rec.contig][rec.pos:rec.stop].mean() )
                        rec.info[f"Border_Tier1"] = "Yes"
                else:
                    rec.info["Tier"] = "Tier1"
                    tier1 += 1
                    if anno_pots_bn[rec.contig][rec.pos - safe_border:rec.pos + safe_border].mean() < 1:
                        rec.info[f"Border_Tier1"] = "Yes"
            else:
                tier2 += 1
                rec.info["Tier"] = "Tier2"
                if anno_pots_bn[rec.contig][rec.pos - safe_border:rec.stop + safe_border].mean() > 0:
                    # print(anno_pots_bn[rec.contig][rec.pos:rec.stop].mean())
                    rec.info[f"Border_Tier1"] = "Yes"
            # if anno_pots_bn[rec.contig][rec.pos :rec.stop + 1].mean()==1:
            #     print("lllll")

            vcf_new.write(rec)
        vcf_new.close()
        vcf.close()
        print(num,tier1,tier2,tier2_a)
rule merge_small:
    input:
        vcf=expand(dir_work + "02_variants/release/ChineseQuartet.twin.germline.{{ref_hap}}.{var_type}.freeze.v2.1.vcf.gz",var_type=["SNV","Indel"]),
        vcf_idx=expand(dir_work + "02_variants/release/ChineseQuartet.twin.germline.{{ref_hap}}.{var_type}.freeze.v2.1.vcf.gz.tbi",var_type=["SNV","Indel"]),
    output:
        dir_work + "02_variants/release/ChineseQuartet.twin.germline.{ref_hap}.SNVIndel.freeze.v2.1.vcf.gz",
    run:
        shell("{bcftools} concat -Oz -a   -o  {output} {input.vcf}")


rule check_SV:
    input:
        bed=dir_work + "02_variants/release/ChineseQuartet.twin.germline.{ref_hap}.SV_Tier1_draft.v2.1.bed",
        vcf=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.draft.vcf.gz",
        vcf_idx=dir_work + "02_variants/confidence/ChineseQuartet.twin.germline.GRCh38.{var_type}.draft.vcf.gz.tbi",
        fai=get_ref_fai,
    output:
        vcf=dir_work + "02_variants/release/ChineseQuartet.twin.germline.{ref_hap}.{var_type}.freeze.v2.1.vcf.gz",
    threads: 4
    wildcard_constraints:
        var_type="DEL|INS"
    run:
        safe_border = 500
        vcf = pysam.VariantFile(f"{input.vcf}")
        new_header_records = [r for r in vcf.header.records if
                              not (r.key == 'bcftools_concatCommand')]
        new_header = pysam.VariantHeader()
        for record in new_header_records:
            new_header.add_line(str(record))
        # new_header.add_line(f'##INFO=<ID=Border.{safe_border}_Tier1,Number=1,Type=String,Description="Is this mutation occurring at the Tier1 region boundary?.">')
        for i in vcf.header.samples:
            new_header.add_sample(i)
        vcf_new = pysam.VariantFile(f"{output.vcf}",header=new_header,mode="w")

        chroms = [f"chr{i}" for i in range(1,23)] + ["chrX"]
        chroms_lens = {}
        for item in open(f"{input.fai}"):
            chrom, chrom_len = item[:-1].split("\t")[:2]
            chroms_lens[chrom] = int(chrom_len)
        anno_pots_bn = {i: np.zeros(chroms_lens[i]) for i in chroms}
        for line in open(f"{input.bed}"):
            chrom, start, end = line[:-1].split("\t")[:3]
            start, end = int(start), int(end)
            anno_pots_bn[chrom][start:end + 1] = 1
        num = 0
        tier1 = 0
        tier2 = 0
        tier2_a = 0
        tier2_regions = []
        for rec in vcf.fetch():
            # if abs(rec.info["SVLEN"]) < 50: continue
            num += 1

            rec.info[f"Border_Tier1"] = "No"
            # if rec.info["Benchmarks_Ratio"] < 1:
            #     print(rec.info["Benchmarks_Ratio"])
            if rec.info["Tier"] in "Tier1":
                if anno_pots_bn[rec.contig][rec.pos:rec.stop].mean() < 1:
                    rec.info["Tier"] = "Tier2"
                    if anno_pots_bn[rec.contig][rec.pos - safe_border:rec.stop + safe_border].mean() > 0:
                        # print(anno_pots_bn[rec.contig][rec.pos:rec.stop].mean() )
                        rec.info[f"Border_Tier1"] = "Yes"
                else:
                    rec.info["Tier"] = "Tier1"
                    tier1 += 1
                    if anno_pots_bn[rec.contig][rec.pos - safe_border:rec.pos + safe_border].mean() < 1:
                        rec.info[f"Border_Tier1"] = "Yes"
            else:
                tier2 += 1
                rec.info["Tier"] = "Tier2"
                if anno_pots_bn[rec.contig][rec.pos - safe_border:rec.stop + safe_border].mean() > 0:
                    # print(anno_pots_bn[rec.contig][rec.pos:rec.stop].mean())
                    rec.info[f"Border_Tier1"] = "Yes"
            vcf_new.write(rec)
        vcf_new.close()
        vcf.close()
        print(num,tier1,tier2,tier2_a)




rule extend_regions_CSV:
    input:
        bed=dir_work_old + "02_variants/benchmark_regions/SV_regions/twin.{ref_hap}.CSV.extend.csv"
    output:
        bed=dir_work + "02_variants/confidence/Filter_regions/ChineseQuartet.twin.germline.{ref_hap}.CSV.filter_forsmall.bed",
    run:
        out_file = open(f"{output}","w")
        for line in open(f"{input}"):
            chrom, start, end = line[:-1].split("\t")[:3]
            start, end = int(start), int(end)
            out_file.write(f"{chrom}\t{start}\t{end}\n")
        # print("lllll",line)
        out_file.close()


# def get_repeat_input(wildcards):
#     return beds[wildcards.repeat]


# rule get_repeat_bed:
#     input:
#         get_repeat_input
#     output:
#         bed=dir_work + "02_variants/confidence/Filter_regions/GRCh38.{repeat}.bed",
#     run:
#         shell(f"gzip -d -c {input} >{output}")

rule merge_repeat_regions:
    input:
        expand(dir_work + "02_variants/confidence/Filter_regions/GRCh38_{repeat}_forSmall.bed",
            repeat=["STR", "VNTR", "SR", "SD"])
    output:
        dir_work + "02_variants/confidence/Filter_regions/ChineseQuartet.twin.repeat.filter_forSmall.bed"
    run:
        shell("""cat {input} |{bedtools} sort |{bedtools} merge >{output}""")

rule merge_repeat_regions2:
    input:
        expand(dir_work + "02_variants/confidence/Filter_regions/GRCh38_{repeat}_forSmall.bed",
            repeat=["STR", "VNTR", "SD", "SR"])
    output:
        dir_work + "02_variants/confidence/Filter_regions/ChineseQuartet.twin.repeat.filter_forSmall2.bed"
    run:
        shell("""cat {input} |{bedtools} sort |{bedtools} merge >{output}""")


rule tabix:
    input:
        "{prefix}.vcf.gz"
    output:
        "{prefix}.vcf.gz.tbi"
    run:
        shell("{tabix} {input}")
