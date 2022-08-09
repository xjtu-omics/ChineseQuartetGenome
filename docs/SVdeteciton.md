## Mutation detection of HiFi reads 

### SNV and Indel detection  


### Structural variation (SV) detection
CuteSV: 
* version: 1.0.11
* script
   ```
  cuteSV --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 -s 2 \
     --diff_ratio_merging_DEL 0.5 --max_cluster_bias_DEL 1000 --threads {threads}  \
     --sample {sample} --min_size 20 \
     {input.bam} {input.ref} {vcf} {workdir}
    ```
* filtering: 
    ```
  ```

Sniffles 
* version: 1.0.12

SVision
* version: 1.3.6

pbsv: 
* version: 2.6.2 (commit SL-release-10.1.0-175-g866ef3c)


SURVIVOR:
* version: 1.0.7
* 


## Mutation detection of Illumina reads 
### SNV and Indel detection for Illumina reads 


### Structural variation (SV) detection Illumina reads 


# SV Detection
## Illumina short reads 
  ### Pindel (v)
    ` pindel `
  ### manta  (v)
    `manta `
  ### delly (v)
    `delly `
  ### lumpy (v)
    `lumpy  `
    
    
