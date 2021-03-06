

Explaining the ceRNA effects from genomic variants
=============================================
**Jing; Hua-sheng; Pavel**
1. Data
------------------------------
All data are downloaded from TCGA, snp level2 is controlled access. WXS bam can access using lab key. [data table](data_1.png) 

- rnaSeq2:  tested; not used in model because the current methods for DEG is basically  using raw countslevel 3; rsem score for RNA seq, de nove assembly. 
 
- rnaSeq:  tested and used in edgeR for DEG
  level 3; rowcounts and rpkm value for each gene 839 - 91(duplicated)  - 
  
- snpArray:  level 2 birdseed data in hg18 annotation; for snp data usage
   [birdsuite algorithm description][bs]; Affymetrix SNP 6.0  
  
- somaticMutation:  level 3 called somatic variants data; used in modelling

- CNV_snparray:  level 3:  used in later modelling

- wxs:  bam file list of whole exome sequencing breast cancer in TCGA for the samples involved in the analysis(731) 


2. Methods
---------------------------------

- [current work flow](workflow.png)
- Differential expressed genes
   - EdgeR: using rowcounts as expression
      - TMM normalization
      - figures here to show the whole process  
      
- snp: confidence filtering: necessary to make sure all the SNPs are real SNP < 0.01 sample size > 0.1  total(731)
  - MAF fréquence filtration: 
  - MAF top small would make no sense to do the calculation  > 0.05 
  - [Hardy–Weinberg principle][hwe]: to make sure the sample comes from a normal population, no need to use in cancer population.
  
- cnv : useing the region position 

- somatic mutations: 

- indel matrix: ongoing

- **[genomic variants count](data_summary.txt)**

3. hierachical Model of genomic variants's effect on gene expression 
-------------------------------------------------------
- [current model](model.png)

- [Kruskal–Wallis test][kw]:  table: [kw filtration snp count: All genes](kwtestSummary.xlsx)


           
- ftest:(-SNP/indel) contribution
- group lasso regression --<still testing other model...>
  - [cohen's kappa][ck] similarity distance 
  - svd: dimension that catch 80% variation 
  - grouping of SNPs: kmeans
- regression
  - logit transform of expression 
  - fitting 
  - random sampling residual 
  - permutation 
  - p value as well as coeff
- filtering:  p value < 0.005 for 1000 permutation

4. bridging gene-variants association with gene-ceRNA effects
---------------------------------------------
- Thoughts

5. population genomic variants
------------------------------------------------
- 1000 Genome SNP / Indel
- UK10K 

6. Test Run using GWAScatalog related genes
----------------------------------------
- [GWAS catalog genes](GWAS_catalog_brca_allGeneName.txt): 70 in total, 27 was found in the having the filtered snp around +/- 1M of the gene TSS

- [Gene-SNP pairs] : 3767 for KW test cutoff 0.01, 1555 for KW test cutoff 1e-6

- snp significantly contribution [genes](fTest_pval_all__snp_cnv_som.txt.sig): relatively lose currently pvalue 0.05
  
- Find [snps] significantly associate with gene: 3 genes
  - [CCND1](grplasso_coeff_grplasso_CCND1.txt)
  - [GDI2](grplasso_coeff_grplasso_GDI2.txt)
  - [NT5C1B](grplasso_coeff_grplasso_NT5C1B.txt)
- ceRNA effect for those genes

- some interesting SNPs
  - [RANBP9 genetic variants network](figure/RANBP9_geneticVars_net.pdf) 
      - rs16874698 [[ref 1 mooney lab]](http://mutdb.org/cgi-bin/mutdb.pl?id=CD83&geneid=9308)       [[ref 2 umbc]](http://bioinf.umbc.edu/dmdm/mut_on_prot.php?id=232223&range=174.25_184.5)
  - [RANBP9 ceRNA network] (figure/RANBP9_ceRNA_net.pdf)
    
[ck]: http://en.wikipedia.org/wiki/Cohen's_kappa
[kw]: http://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance
[hwe]:http://en.wikipedia.org/wiki/HWE
[bs]: http://www.broadinstitute.org/science/programs/medical-and-population-genetics/birdsuite/birdsuite-analysis#birdsuite_snps
