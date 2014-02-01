

Explaining the ceRNA effects from genomic variants
=============================================

1. Data
------------------------------
All data are downloaded from TCGA, snp level2 is controlled access. WXS bam can access using lab key.
![TCGA data](/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/report/data_1.png?raw=true) <imd src="/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/report/data_1.png">


- rnaSeq2:  tested; not used in model because the current methods for DEG is basically  using raw countslevel 3; rsem score for RNA seq, de nove assembly. 
  [birdsuite algorithm description][bs]; Affymetrix SNP 6.0  
  
- rnaSeq:  tested and used in edgeR for DEG
  level 3; rowcounts and rpkm value for each gene 839 - 91(duplicated)  - 
  
- snpArray:  level 2 birdseed data in hg18 annotation; for snp data usage

- somaticMutation:  level 3 called somatic variants data; used in modelling

- CNV_snparray:  level 3:  used in later modelling

- wxs:  bam file list of whole exome sequencing breast cancer in TCGA for the samples involved in the analysis(731) 

- WGS:  whole genome bams of TCGA breast cancer, only 86 samples. not used

- miRNA:  project stage 2 data; for future use.

2. data preprocessing 
---------------------------------
- differential expressed genes
   - voom transformation of raw count
   
   - EdgeR: using rowcounts as expression
      - TMM normalization
      - figures here to show the whole process
      
- snp: confidence filtering: necessary to make sure all the SNPs are real SNP < 0.01 sample size > 0.1  total(731) 78043 out

  - MAF fréquence filtration: 
  
  - MAF top small would make no sense to do the calculation  > 0.05 98 out
  
  - [Hardy–Weinberg principle][hwe]: to make sure the sample comes from a normal population, no need to use in cancer population.
  
  - Annotation use platform information(life over from hg18 to hg19), lost 5808; retain 822651
  
  - [Kruskal–Wallis test][kw]: threshold is an important thing
      #plot showing how many pairs of testing we are doing, and how many we are selecting using different threhold.
```
```
- cnv : useing the region position 

- somatic mutations: 

- indel matrix: ongoing

3. hirarchical Modeling of genomic variants
-------------------------------------------------------
- ftest:(-SNP/indel) contribution
  - simple linear model and lose threshold 
      - sig
      - no sig

- group lasso regression --<still testing other model...>
  - [cohen's kappa][ck] similarity distance 
  - svd: dimension that catch 80% variation 
  - group:
      - kmeans:
      - DBSCAN
- regression

  - logit transform
  - fitting 
  - random sampling residual 
  - permutation 
  - p value as well as coeff

- filtering:
  - p value < 0.05
  - .rda file created for all sig genes


4. bridging gene-variants association with gene-ceRNA effects
---------------------------------------------


5. population genomic variants
------------------------------------------------
- 1000 Genome SNP
- 1000 Genome Indel
- UK10K 

[ck]: http://en.wikipedia.org/wiki/Cohen's_kappa
[kw]: http://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance
[hwe]:http://en.wikipedia.org/wiki/HWE
[bs]: http://www.broadinstitute.org/science/programs/medical-and-population-genetics/birdsuite/birdsuite-analysis#birdsuite_snps
