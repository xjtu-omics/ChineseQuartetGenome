# Quick start

## Run with snakemake

### 1. Prepare you environment

The following software/packages are required in same environment: 
  * python
  * snakemake
  * pysam 
  * numpy 
  * pandas 

You can use conda to install all of these packages. 

You also need to install the following software and  
  * bedtools 
  * bcftools 
  * tabix 
  * bgzip
  * [hap.py](https://github.com/Illumina/hap.py) (for small variants benchmarking)
  * [truvari](https://github.com/ACEnglish/truvari) (for structural variants benchmarking)

### 2. Prepare the benchmark regions and variants

### 3. Config you task

### 4. Run the pipeline

### 5. Check the result



## Run with docker


# Known and Likely Limitations of this version

* For assemblies:
    * 20 original gaps in hifiasm assembly have abnormal read depths.
* For variants:
    * SV detected by only one technology need to be validated in the future.
    * Only variants at autosome are benchmarking. 


