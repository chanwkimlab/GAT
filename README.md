Generic Association Tool `1.0`
========================================

## Functionality
* Flexible in defining multialleic markers
    * For markers of multiple allelels, the reported p-value is from deviance from the null model.
    * It can be used to test gene amino acid polymorphism (ex. HLA gene)
* Support for phased markers. Two haplotypes from individuals are consided to have independent effect on phenotype.
* Faster than R with same simulation settings and parameters
    * Manually convert categorial input to onehot-encoded input (https://stackoverflow.com/questions/48898712/very-slow-glm-logistic-regression-in-r)
    * Numpy use algorithm optimized for cpu (BLAS)
    * Built-in function(data.table) for loading text file in R is slow.
* Extract individuals in .phe from .bgl.phased and continue association test only with individuals
* Integrity check between .fam, pedigree information from .bgl.phased



## Getting Started

1. `cd ANY_DESIRED_PATH`
2. `git clone https://github.com/ch6845/Generic_Association_Tool`
3. `export PATH=$PATH:$ANY_DESIRED_PATH/Generic_Association_Tool`


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
* Log file
    *. log
    
    
## Example
```
Python GAT.py
--assoc linear
--out data/out_assoc/ALP/step_01
--bgl_phased data/genotype/4_merge/HLA_AA_SNP.bgl.phased
--bfile data/genotype/4_merge/HLA_SNP_1000G
--multialleic (?P<name>HLA_[0-9A-Z]*)\*(?P<allele>[0-9:]*)
--multialleic_always (?P<name>AA_[A-Z0-9]*_[0-9]*_[0-9]*_exon[0-9]*)_*(?P<allele>[A-Z]*)
--pheno data/out_pheno/ALP.phe
--covar data/out_assoc/ALP/step_01.covar
--condition_list data/out_assoc/ALP/step_01.cond
```
    
## Dependencies    
* numpy              1.17.4 
* pandas             0.25.3   
* scipy              1.3.2 
* statsmodels        0.10.1
* pyplink            1.3.5

## Citation