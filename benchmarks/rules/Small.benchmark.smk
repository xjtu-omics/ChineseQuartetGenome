# ======================================================================================================================
# Project: ChinseQuartetGenome
# Script : Small.benchmark.smk TODO check 
# Author : Peng Jia
# Date   :  2023/2/21
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================






rule hap_eval:
    input:
        a="a.txt",
        b="b.txt"
    output:
        c="c.txt"
    run:
        shell("touch {input}")
