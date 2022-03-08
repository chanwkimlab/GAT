GAT (Generic Genome-Wide Association Tool) `1.0`
========================================


`GAT` is a generic tool for conducting GWAS(Genome-Wide Association Study) on phased/unphased(=dosage) biallelic/multiallelic markers simultaneously.

## Getting Started

If you have used `PLINK` before, you will find it easy to use `GAT`.

1. `cd ANY_DESIRED_PATH`
2. `git clone https://github.com/ch6845/GAT`
3. `export PATH=$PATH:$ANY_DESIRED_PATH/GAT`
    
```
python ANY_DESIRED_PATH/GAT/GAT.py \
python /data/ch6845/MHC_phewas_testbench/Generic_Association_Tool/GAT.py \
--assoc linear \
--out step_02.GAT \
--bgl-phased /data/ch6845/MHC_phewas_testbench/data/genotype/4_merge/KCHIP_HLA_AA_SNP.bgl.phased \
--bfile /data/ch6845/MHC_phewas_testbench/data/genotype/4_merge/1000G \
--multiallelic "(?P<name>HLA_[0-9A-Z]*)\*(?P<allele>[0-9:]*)" \
--multiallelic-always "(?P<name>AA_[A-Z0-9]*_[\-0-9]*_[0-9]*_exon[0-9]*)_*(?P<allele>[A-Z]*)" \
--skip "(?P<name>6:[0-9]*_[A-Z]*/[\<\>A-Z\:0-9]*),(?P<name>AX\-[0-9]*),(?P<name>AFFX\-SP\-[0-9]*),(?P<name>SNPS_.*),(?P<name>INS_SNPS_.*)" \
--pheno /data/ch6845/MHC_phewas_testbench/data/out_pheno/height.phe \
--covar /data/ch6845/MHC_phewas_testbench/data/out_assoc/height/covar \
--condition-list /data/ch6845/MHC_phewas_testbench/data/out_assoc/height/step_02.cond


plink2 \
--glm \
--bfile /data/ch6845/MHC_phewas_testbench/data/genotype/4_merge/KCHIP_SNP_1000G \
--pheno /data/ch6845/MHC_phewas_testbench/data/out_pheno/height.phe \
--out step_02.plink \
--covar step_02.GAT.covar_unphased.tsv \
--covar-variance-standardize \
--threads 40

```

## Motivation
For multiallelic(categorical) marker, to get the significance of the overall effect of alleles on the phenotype, statistical test called omnibus test is conducted.(https://en.wikipedia.org/wiki/Omnibus_test, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6292650/)
```
Our univariate tests of binary SNP and SNP allele markers,
and our omnibus tests of polymorphic HLA amino acid positions both highlighted HLA-DRβ1 amino acid position 11 as the MHC feature most significantly associated with UC.
(from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3341846)
```
In some cases, phased genotype is required for association study. When we condition a marker that may have some correlation with other variants(especially markers in same gene), distinguishing each haplotype may show a more reliable, larger-power result. (https://www.nature.com/articles/s41398-017-0010-9) Two haplotypes from individuals are consided to have independent effect on phenotype. For example, it can be used to condition amino acid polymorphism in a gene (ex. HLA gene)

In reality, both phased genotype file and unphased genotype file can coexist for one dataset. Testing the variants of such diverse cases(biallelic-multiallelic/phased-unphased) simultaneously is very labor-intensive and time-consuming especially when conducting stepwise conditional analysis.
For example, 
peak from phased ->  binary(phased->dosage), omnibus(condition)
peak from unphased -> binary(condition=dosage), omnibus(dosage)
gene->


Plink(latest version 2), one of the most popular tool for GWAS, does not support association for multi-allelic variants, although it takes categorical covariate as input.
In some conditional analysis, if a variant having significant signal is in transcription region of a gene, the variants in the gene is together used as covariate to robustly attribute the signal to the gene. 
```
When the top-associated variant itself was the HLA gene polymorphism or the SNV and indel in strong LD with any of the HLA gene polymorphisms (r2 ≥ 0.7), 
we additionally included all the two-digit, four- digit and six-digit alleles and the amino acid polymorphisms of the corresponding HLA gene
as covariates in the regression to robustly condition the associations attributable to the HLA gene, as previously described
Genetic and phenotypic landscape of the major histocompatibilty complex region in the Japanese population.
(Hirata, Jun, et al. "Genetic and phenotypic landscape of the major histocompatibilty complex region in the Japanese population." Nature genetics 51.3 (2019): 470-480.)
```
Then, the association model has a number of covariates. Unfortunately, --glm in plink2.0 gives up continuing association test when covariate have some colinearity. (VIF) However, GLM in `Statsmodels` package(https://www.statsmodels.org) used by this tool partly automatically resolves the situation.

Therefore, specialized tool for this job was needed and this tool was developed.


## Functionality
* Conduct GWAS(Genome-Wide Association Study) for phased/unphased(=dosage) biallelic/multiallelic markers simultaneously.
* Flexible in defining multiallelic markers (by regular expression)
* Faster than R with same simulation settings and parameters
    * Manually convert categorial input to onehot-encoded input (https://stackoverflow.com/questions/48898712/very-slow-glm-logistic-regression-in-r)
    * Numpy use algorithm optimized for Intel CPU (`BLAS`)
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
* Regular expression for identifying multiallelic marker
    * --multiallelic regex
    * --multiallelic-always regex
* (Optional) .covar (the same as plink)
    --covar .cover
* (Optional) .condition (the same as plink --condition-list)
    --condition-list .cond
* (Optional) Regular expression for identifying markers to skip
    You can exclude bi-allelic markers from analysis and run association test using plink, since glm in plink is faster
    --skip .cond    
    
### Ouput:
* Association report
    p-value for each merged marker (AA_A_1_Residue->AA_A_1) (HLA_A\*--:-- ->HLA_A)  chisquare - degree freedom - p-value
    *.tsv
    for biallelic markers, report (effect size(glm regression coeffient)/z score/p value)
    for multiallelic markers, report (chisquare/degree of freedom/p value) p-value is from deviance from the null model.
* Log file
    *. log

    
## Dependencies    
* `numpy`              1.17.4 
* `pandas`             0.25.3   
* `scipy`              1.3.2 
* `statsmodels`        0.10.1
* `pyplink`            1.3.5

## Citation
If you use any part of this code, please cite our
[paper](https://doi.org/10.1093/hmg/ddac016), in which we introduced this software package.
```
@article{kim2022phenome,
  title={Phenome-wide association study of the major histocompatibility complex region in the Korean population identifies novel association signals},
  author={Kim, Chanwoo and Kim, Young Jin and Choi, Wanson and Jang, Hye-Mi and Hwang, Mi Yeong and Jung, Sunwoo and Lim, Hyunjoon and Hong, Sang Bin and Yoon, Kyungheon and Kim, Bong-Jo and others},
  journal={Human molecular genetics},
  year={2022}
}
```

## Contact
If you have any inquiries, please feel free to contact
- [Chanwoo Kim](https://chanwoo.kim) (Seoul National University)



