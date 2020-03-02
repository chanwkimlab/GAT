GAT (Generic Genome-Wide Association Tool) `1.0`
========================================


`GAT` is a generic tool conducting GWAS(Genome-Wide Association Tool) for phased/unphased(=dosage) bialleic/multialleic markers simultaneously.

## Getting Started

1. `cd ANY_DESIRED_PATH`
2. `git clone https://github.com/ch6845/GAT`
3. `export PATH=$PATH:$ANY_DESIRED_PATH/GAT`
    
```
GAT.py \
--assoc linear \
--out sample_output \
--bgl-phased /data/ch6845/MHC_phewas_testbench/data/genotype/4_merge/HLA_AA_SNP.bgl.phased \
--bfile /data/ch6845/MHC_phewas_testbench/data/genotype/4_merge/HLA_SNP_1000G \
--multialleic (?P<name>HLA_[0-9A-Z]*)\*(?P<allele>[0-9:]*) \
--multialleic-always (?P<name>AA_[A-Z0-9]*_[\-0-9]*_[0-9]*_exon[0-9]*)_*(?P<allele>[A-Z]*) \
--pheno /data/ch6845/MHC_phewas_testbench/data/out_pheno/ALP.phe \
--covar /data/ch6845/MHC_phewas_testbench/data/out_assoc/ALP/step_01.covar \
--condition-list /data/ch6845/MHC_phewas_testbench/data/out_assoc/ALP/step_01.cond\
```

## Why is it needed
For multialleic(categorical) marker, to get the significance of the overall effect of alleles on the phenotype, statistical test called omnibus test is conducted.(https://en.wikipedia.org/wiki/Omnibus_test, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6292650/)
```
Our univariate tests of binary SNP and SNP allele markers, and our omnibus tests of polymorphic HLA amino acid positions both highlighted HLA-DRβ1, amino acid position 11 as the MHC feature most significantly associated with UC.
(from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3341846)
```
In some cases, assoicating phased genotype is required. When we condition a marker that may have some correlation with other variants(especially markers in same gene), distinguishing each haplotype may have a more reliable, larger-power result. (https://www.nature.com/articles/s41398-017-0010-9) Two haplotypes from individuals are consided to have independent effect on phenotype. For example, it can be used to condition amino acid polymorphism ina gene (ex. HLA gene)

In reality both phased genotype file and unphased genotype file exist for one dataset. Testing the variants of such diverse cases simultaneously is very labor-intensive and time-consuming.
For example, 
peak from phased ->  binary(phased->dosage), omnibus(condition)
peak from unphased -> binary(condition=dosage), omnibus(dosage)
gene->


Plink(latest version 2), one of the most popular tool for GWAS, does not support association for multi-alleic variants, although it takes categorical covariate as input.
In some conditional analysis, if a variant having significant signal is in transcription region of a gene, all of the variants in the gene is used as covariate to robustly attribute the signal to the gene. Then, a number of covariates are addded to the association model. However, --glm in plink2.0 gives up continuing association test when covariate have some colinearity. (VIF) glm in statsmodels used by this tool partly automatically resolves the situation.

Therefore, specialized tool for this job is needed.


## Functionality
* Conduct GWAS(Genome-Wide Association Tool) for phased/unphased(=dosage) bialleic/multialleic markers simultaneously.
* Flexible in defining multialleic markers (by regular expression)
* Faster than R with same simulation settings and parameters
    * Manually convert categorial input to onehot-encoded input (https://stackoverflow.com/questions/48898712/very-slow-glm-logistic-regression-in-r)
    * Numpy use algorithm optimized for cpu (BLAS)
    * Built-in function(data.table) for loading text file in R is slow.

## Reference

https://gatkforums.broadinstitute.org/gatk/discussion/6455/biallelic-vs-multiallelic-sites
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6292650/




## Input / Output

### Input: 
* Genotype
    Either or both of the following genotype files are required.
    If a marker coexists on both of the files, priority is on --bgl.phased
    * --bfile .bed/bim/fam
    * --bgl-phased .bgl.phased (You can convert plink file to .bgl.phased using Beagle 5)
* Phenotype
    * --pheno .phe (the same as plink. values set as -9 is considered missing)
* Regular expression for identifying multialleic marker
    * --multialleic regex
    * --multialleic-always regex
* (Optional) .covar (the same as plink)
    --covar .cover
* (Optional) .condition (the same as plink --condition-list)
    --condition-list .cond
    
### Ouput:
* Association report
    p-value for each merged marker (AA_A_1_Residue->AA_A_1) (HLA_A\*--:-- ->HLA_A)  chisquare - degree freedom - p-value
    *.tsv
    for bialleic markers, report (effect size(glm regression coeffient)/z score/p value)
    for multialleic markers, report (chisquare/degree of freedom/p value) p-value is from deviance from the null model.
* Log file
    *. log
    

    
## Dependencies    
* numpy              1.17.4 
* pandas             0.25.3   
* scipy              1.3.2 
* statsmodels        0.10.1
* pyplink            1.3.5

## Citation