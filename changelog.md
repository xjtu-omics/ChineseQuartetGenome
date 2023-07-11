Please see the latest information of the project at https://github.com/xjtu-omics/ChineseQuartetGenome

## v2.1
### Release notes
* Defined a more exclusive Tier 1 benchmark regions.

* For SVs, high-confidence calls (Tier 1) and corresponding regions of v2.1 of the benchmark were obtained according to the following steps:
  * Initially, we identified candidate high-confidence SVs within benchmark regions. This included all SVs supported by two or more technologies. For technology-specific SVs, we considered the following three criteria to classify them as candidate high-confidence SVs:    
    *   They were not located within segmental duplications, VNTR (variable number tandem repeat), or simple repeat regions. 
    *	They were not reported through short-read analysis.
    *	They were reported in an orthogonal callset.
  *	Next, we refined the benchmark regions to create Tier 1 (v2.1) regions by excluding the following: 
       *  Regions spanning SVs that did not meet the criteria in step a) were excluded.
       *  	Non-diploid regions were excluded.
       *	Simple repeats that were overlap with variants not meet the criteria in step a) were excluded.
  * The candidate high-confidence SVs fully within the Tier 1 regions were defined as high-confidence SVs.

* For small variants, high-confidence calls (Tier 1) of v2.1 were obtained according to the three steps:
  * Initially, we identified candidate high-confidence small variants within benchmark regions. This included all small variants supported by all three technologies. For variants supported by two technologies, we considered those were not located within segmental duplications, VNTR, simple repeat, or short tandem repeat regions as candidate high-confidence calls. 
  * Next, we refined the benchmark regions to create Tier 1 (v2.1) regions for small variants by excluding the following:
    *	Regions spanning SVs were excluded.
    *	Non-diploid regions were excluded.
    *	Regions spanning small variants that did not meet the criteria in step a) were excluded. 
    *	Repeat regions (Simple repeats, segmental duplications, short tandem repeats, and VNTRs) that were overlap with small variants not meet the criteria in step a) were excluded.
  *	The candidate high-confidence small variants fully within the Tier 1 regions were defined as high-confidence small variants.

* In this release, to determine whether a variant belongs to a specific region,
we require that the variant is 100% within that region. 
Users have the flexibility to adjust this parameter (_SV_overlap_ratio_ and _Small_overlap_ratio_) through the configuration file. 
* In the query set, we also include variants in the bordering regions of a given region, 
such as Tier 1. This inclusion helps mitigate false negatives caused by inaccurate 
variant breakpoints. Users can manually verify the accuracy of variants 
in these border regions. When calculating false positives, variants in the 
bordering regions are excluded. Users can also adjust parameter _size_tolerate_region_boundary_SV_ and
_size_tolerate_region_boundary_Small_ in light of requirements.



## v2.0
### Release notes
* Defined benchmark regions.
* Defined phased regions. 
* Defined high-confidence regions. 
* Annotated the benchmark calls in the above regions.
* Release a pipeline for user to use variant benchmarks.

### Known and Likely Limitations of v2.0
* For assemblies:
    * 20 original gaps in hifiasm assembly have abnormal read depths.
    * 965 and 1,464 novel sequence regions in paternal and maternal, respectively, are with abnormal read depth. 
* For variants:
    * SV detected by only one technology need to be validated in the future.
    * Only variants at chr1-chr22, and chrX are benchmarking. 

    
# v1.0
### Release note 
* First release of the Chinese Quartet benchmark
* Only variants at chr1-chr22 are benchmarking. 
* First release of the Chinese Quartet benchmark