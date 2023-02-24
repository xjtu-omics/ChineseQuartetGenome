# Chinese Quartet Genome

The [Chinese Quartet Project](http://chinese-quartet.org/) is designed for quality control in biological analysis in China. The samples 
in this project have been approved as standard reference materials by 
the [State Administration for Market Regulation in China](https://www.samr.gov.cn/). In this repository,
we provided high-quality genome and variant benchmarks for this quartet.

## Samples 
There are four samples in Chinse Quartet, including father (LCL7), mother (LCL8), and two monozygotic 
twin daughters (LCL5 and LCL6). This quartet family is from Taizhou in China. For more information of the quartet, 
please click [here](http://chinese-quartet.org/). 



## Data information
The Quartet multi-omics reference materials and data are publicly available and accessible. 
The recipients of the Reference Materials are highly encouraged to share their data with the 
community in order for to improve the evaluation of the technologies, pipelines, batch effects, and so on.

For genomic data, we generated 50x PacBio HiFi, 100x coverage of Oxford Nanopore, 30x coverage of ultra 
long Oxford Nanopore (only for LCL5), 100x PacBio CLR, 100x 10x genomics linked reads and BioNano. 
All raw data is available from public, if you want to use genomic data in this study, please fell free to contact with us.
## Assembly 
* Details are available at [github](/docs/assm_stat.md)
* Data from two monozygotic twin daughters were merged for assembly.
* Long reads were phased into two haplotypes with heterozygous mutations.
* HiFi reads were assembled with hifiasm, hicanu, and flye.
* ONT reads were assembled with flye and shasta.
* Pipeline for assembly is available at https://github.com/PengJia6/AssmPipe.
* Pipeline for assembly evaluation is available at https://github.com/PengJia6/Postassm

## Variant Benchmarks
* Variants are discovery with Illumina reads, HiFi reads, and haplotyped-resolved assemblies.
* Varaints are filtered with read depth, allele depth, Mendelian rules.
* Methods for variant establishing is available at https://github.com/PengJia6/NGSGemlineMutPipe and https://github.com/PengJia6/TGSGermlineMutPipe

## Availability
Haplotype-resolved assemblies and variant benchmarks of Chinese Quartet are available at [OneDrive](https://stuxjtueducn-my.sharepoint.com/:f:/g/personal/pengjia_stu_xjtu_edu_cn/Eqc2HjImbKFHiJHluAbLH68Bm6wzmY25v48y3kjVg5iRvg?e=g784Aw).

## Citation
  Jia P, Dong L, Yang X, Wang B, Wang T, Lin J, Wang S, Zhao X, Xu T, Che Y, et al: Haplotype-resolved assemblies and variant benchmark of a Chinese Quartet. bioRxiv 2022:2022.2009.2008.504083. [PDF](https://www.biorxiv.org/content/10.1101/2022.09.08.504083v1.full.pdf)

## Contributions 

* Software developers are encouraged to benchmark their software using our samples and data.
* We also encourage the community to submit their assemblies and variant to improve the assemblies and variant benchmarks. 
* Please contact with Peng Jia (pengjia@stu.xjtu.edu.cn)



## Contact

* Kai Ye (kaiye@xjtu.edu.cn)
* Peng Jia (pengjia@stu.xjtu.edu.cn)
