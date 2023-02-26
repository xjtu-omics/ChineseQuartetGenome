# Pipeline for variant benchmarking of Chinese Quartet. 
We develop a pipeline for benchmarking of variants in the genomes of 
Chinese Quartet, enabling evaluation of the performance of 
different sequencing technologies, variant calling algorithms, and pipelines.

## Quick start

### Run with snakemake

#### 1. Prepare you environment

The following software/packages are required in same environment: 
  * python
  * snakemake
  * pysam 
  * numpy 
  * pandas 

You can use conda to install all of these packages, for example:   

   ```shell
    conda install package-name 
   ```
You also need to install the following software:
  * bedtools 
  * bcftools 
  * tabix 
  * bgzip
  * [hap.py](https://github.com/Illumina/hap.py) (for small variants benchmarking)
  * [truvari](https://github.com/ACEnglish/truvari) (for structural variants benchmarking)

#### 2. Prepare the benchmark regions and variants
Download the latest version of variants and benchmark regions of Chinese Quartet according to the [instruction.](https://github.com/xjtu-omics/ChineseQuartetGenome#variant-benchmarks) 

#### 3. Config you task
* Config your own config.yaml according to the [template](/conf/config.yaml).
* Config your own vcf file in a tsv (Tab-Separated-Values) file according to the [template](/conf/samples.tsv). 
  

#### 4. Run the pipeline
Run the piepline with snakemake 

```shell
snakemake -s ./Sankefile -j 40 -k --ri # on a local computer
snakemake -s ./Sankefile -j 10 -k --ri  --cluster 'qsub -l nodes=1:ppn=12 -l walltime=99:00:00' >sublog 2>&1 & # on a cluster 
```


### Run with docker

Under construction!

## Citation
  Jia P, Dong L, Yang X, Wang B, Wang T, Lin J, Wang S, Zhao X, Xu T, Che Y, et al: Haplotype-resolved assemblies and variant benchmark of a Chinese Quartet. bioRxiv 2022:2022.2009.2008.504083. [PDF](https://www.biorxiv.org/content/10.1101/2022.09.08.504083v1.full.pdf)



## Contact

* Kai Ye (kaiye@xjtu.edu.cn)
* Peng Jia (pengjia@stu.xjtu.edu.cn)


