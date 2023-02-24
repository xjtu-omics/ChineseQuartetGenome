# ======================================================================================================================
# Project: ChinseQuartetGenome
# Script : SV.benchmark.smk TODO check 
# Author : Peng Jia
# Date   :  2023/2/21
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================

rule truvari_compare:
    input:
        a="a.txt",
        b=dir_work+"{}"
    output:
        c="c.txt"
    run:
        shell("touch {input}")
