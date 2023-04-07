# ======================================================================================================================
# Project: ChinseQuartetGenome
# Script : read_phase.smk TODO check 
# Author : Peng Jia
# Date   :  2023/4/7
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================


rule all:
    input:
        expand("{sample}/{sample}.{tech}.untag.h1.bam",tech=["ONT", "HiFi"],sample=["LCL5", "LCL6"])


rule split_untag_bam:
    input:
        bam="{sample}/{sample}.{tech}.untag.bam",
        bai="{sample}/{sample}.{tech}.untag.bam.bai",
    output:
        hap1="{sample}/{sample}.{tech}.untag.h1.bam",
        hap2="{sample}/{sample}.{tech}.untag.h2.bam",
    run:
        import pysam

        samfile = pysam.AlignmentFile(str(input.bam),"rb")
        hap1_bam = pysam.AlignmentFile(str(output.hap1),"wb",template=samfile)
        hap2_bam = pysam.AlignmentFile(str(output.hap2),"wb",template=samfile)
        read_num = 0
        for read in samfile.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary: continue
            read_num += 1
            if read_num % 2 == 1:
                hap1_bam.write(read)
            else:
                hap2_bam.write(read)
        hap1_bam.close()
        hap2_bam.close()
        samfile.close()

rule bam2fq:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.fq.gz"
    threads: 24
    run:
        shell("{samtools} fastq -0 {output} -@ {threads} {input}")
