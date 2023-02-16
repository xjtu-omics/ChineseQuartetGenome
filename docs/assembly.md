# ChineseQuartet assembly summary 

-----


See details from our preprint [paper.](https://www.biorxiv.org/content/10.1101/2022.09.08.504083v1)
## 1. Read phasing 

First, we obtain the high quality SNVs and Indels from a public [research](https://www.ncbi.nlm.nih.gov/pubmed/34980216) 
of these samples. Next, the variants were phased according to their genetic (parent-child) information. 
The raw reads were aligned to reference genome and phased into two haplotypes by the phased variants. 
The unmapped reads and unphased reads were also assigned to the two haplotypes randomly.

## 2. Assembly

For reads from each haplotype, we assembled them with hifiasm  (v0.15.5), hicanu (v-r10117), flye (v2.8.3-b1695), and shasta (v0.7.0) .
The parameters of each assembler was as follows:
* hifiasm for HiFi reads
```shell
 hifiasm  -o {output}  -t {threads}  {input.fqs}
```

* hicanu for HiFi reads (with purge dups) 

```shell
ulimit -Su 100000 && canu useGrid=false maxThreads={threads}-p {assm_prefix} -d {assm_dir} genomeSize=3.1g -pacbio-hifi {input.fqs}
minimap2 -t {threads} -x map-hifi {assm_dir}/{assm_prefix}.contigs.fasta {input.fqs}|gzip -c - > {output.paf}
pbcstat {output.paf} && calcuts PB.stat >{output.cutoff} 
split_fa {assm_dir}/{assm_prefix}.contigs.fasta > {output.split_fa}
minimap2 -t {threads} -xasm5 -DP {output.split_fa} {output.split_fa} | gzip -c > {output.self_paf}
purge_dups -2 -T {output.cutoff} -c PB.base.cov  {input.self_paf} > {output.bed} 
get_seqs {input.bed} {assm_dir}/{assm_prefix}.contigs.fasta 
``` 

* flye for HiFi reads

```shell
flye --pacbio-hifi {input.fqs} --genome-size 3.1g --out-dir {output} --threads {threads} 
```
* flye for ONT reads

```shell
flye --nano-raw {input.fqs} --genome-size 3.1g --out-dir {output} --threads {threads} 
```
* shasta for ONT reads

```shell
shasta --input {input.fqs}  --threads {threads} --assemblyDirectory {output} --config {shasta_config} --memoryMode filesystem --memoryBacking 2M --command assemble 
```


