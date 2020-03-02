Multiple Alleic Association Test Module
========================================

## Input/Output

### Input: 
* .bgl.phased (You can convert plink file to .bgl.phased using Beagle 5)  (TODO: support for plink file if conditional analysis is not done)
* .fam (the same as plink)
* .phe (the same as plink)
* regular expression for identifying multialleic marker
* (optional) .covar (the same as plink)
* (optional) .condition (the same as plink except that the )
    
### Ouput:
* p-value for each merged marker (AA_A_1_Residue->AA_A_1) (HLA_A\*--:-- ->HLA_A)  chisquare - degree freedom - p-value
    
    
## Usage

`
example
omnibus_test.py    --assoc linear 
    --out data/out_assoc/hba1c/step_02.omnibus 
    --pheno data/out_assoc/hba1c/phenotype.phe 
    --fam data/genotype/4_merge/???.fam 
    --covar data/out_assoc/hba1c/step_02.omnibus.covar.temp 
    --bglphased data/out_assoc/hba1c/step_02.bgl.phased
    --condition-list data/out_assoc/hba1c/step_02.omnibus.cond\
`
    

## Functionality
* can define any marker as multialleic and do omnibus test
    * support HLA gene multialleic test by omnibus test
* faster than R with same simulation settings and parameters
    * manually convert categorial input to onehot-encoded input (https://stackoverflow.com/questions/48898712/very-slow-glm-logistic-regression-in-r)
    * numpy use modules optimized for cpu(blas)
    * data.table module in R is slow
* extract individuals in .phe from .bgl.phased and continue association test only with individuals
* Integrity check between .fam, pedigree information from .bgl.phased

## Dependencies    
* numpy              1.17.4 
* pandas             0.25.3   
* scipy              1.3.2 
* statsmodels        0.10.1
* pyplink            1.3.5

## Citation