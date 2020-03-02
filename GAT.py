#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import re
import sys

import argparse
import logging
import time
import socket


from pyplink import PyPlink
import pandas as pd
import numpy as np
from scipy.stats import chi2

import statsmodels.api as sm

#jupyter nbconvert GAT.ipynb --to script


# In[32]:


def dir_path(path):
    if os.path.dirname(path)=='' or os.path.exists(os.path.dirname(path)):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")

def file_path(path):
    if os.path.isfile(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")

def bfile_path(path):
    if os.path.isfile(path+'.fam') and os.path.isfile(path+'.bed') and os.path.isfile(path+'.bim'):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")
        
parser = argparse.ArgumentParser(description='GAT (Generic Genome-Wide Association Tool) For more information, please see https://github.com/ch6845/GAT')

#required mode
parser.add_argument('--assoc',choices=['linear','logistic'],required=True)

#required output file
parser.add_argument('--out', type=dir_path,required=True,help='output file prefix. (prefix.log, prefix.assoc will be generated)')

#required input files
parser.add_argument('--bgl-phased', type=file_path,help='bgl-phased (See Beagle 5.1 documentation)')
parser.add_argument('--bfile', type=bfile_path,help='plink binary format')
parser.add_argument('--multialleic',type=str,help='regular expression for specifying multiple alleic marker (comma delimiter)')
parser.add_argument('--multialleic-always',type=str,help='regular expression for specifying multiple alleic marker (comma delimiter)')

parser.add_argument('--pheno', type=file_path,required=True,help='format is the same as plink. Tab-delimited file without header, of which the first and second columns is family and within-family IDs respectively, and the third column is pheotype')

#optional
parser.add_argument('--covar', type=file_path,help='format is the same as plink')
parser.add_argument('--condition-list',type=file_path,help='format is the same as plink')


# In[33]:


debug=False
#debug=True

if debug:
    arg_split='--assoc linear --out sample_output --bgl-phased /data/ch6845/MHC_phewas_testbench/data/genotype/4_merge/KCHIP_HLA_AA_SNP.bgl.phased --bfile /data/ch6845/MHC_phewas_testbench/data/genotype/4_merge/KCHIP_HLA_SNP_1000G --multialleic (?P<name>HLA_[0-9A-Z]*)\*(?P<allele>[0-9:]*) --multialleic-always (?P<name>AA_[A-Z0-9]*_[\-0-9]*_[0-9]*_exon[0-9]*)_*(?P<allele>[A-Z]*) --pheno /data/ch6845/MHC_phewas_testbench/data/out_pheno/FEV_predicted.phe --covar /data/ch6845/MHC_phewas_testbench/data/out_assoc/FEV_predicted/step_01.plink.covar --condition-list /data/ch6845/MHC_phewas_testbench/data/out_assoc/FEV_predicted/step_01.plink.cond'.split(' ')
    args=parser.parse_args(arg_split)
else:
    args=parser.parse_args()
    
if args.bfile is None and args.bgl_phased is None:
    raise argparse.ArgumentTypeError("either --bfile or --bgl-phased parameter is needed")    


# In[7]:


log = logging.getLogger('logger')
log.setLevel(logging.DEBUG)

log_file_path=args.out+'.log'
fileHandler = logging.FileHandler(log_file_path,'w')
streamHandler = logging.StreamHandler()

formatter = logging.Formatter('%(message)s')
fileHandler.setFormatter(formatter)
streamHandler.setFormatter(formatter)

log.addHandler(fileHandler)
log.addHandler(streamHandler)

log.info_head=lambda x: log.info('\n'+'*'*int((100-len(x))/2)+x+'*'*int((100-len(x))/2)+'\n')


# In[8]:


log.info_head("*********************************")
log.info("* GAT (Generic Genome-Wide Association Tool)")
log.info("* Description: Generic module for associating bialleic/multialleic phased/unphased markers")
log.info("* version 1.0")
log.info("* (C) 2020-, Seoul National University")
log.info("* Please report bugs to: Chanwoo Kim <ch6845@snu.ac.kr>")
log.info("* https://github.com/ch6845/GAT")
log.info_head("*********************************")

log.info("Start time: "+time.strftime('%c', time.localtime(time.time())))

log.info('Working directory: '+os.getcwd())
log.info('Hostname: '+socket.gethostname())

log.info('Parameters\n'+'\n'.join(['--{} {}'.format(key,value) for key,value in vars(args).items()]))
#print('Parameters\n'+'\n'.join(['--{} {}'.format(key,value) for key,value in vars(args).items()]))


# In[9]:


assoc=args.assoc
out=args.out


# In[10]:


log.info_head("Data Loading")


# # parse input files

# In[11]:


plink=None
plink_bim=None
plink_fam=None

if args.bfile is not None:
    plink=PyPlink(args.bfile)
    plink_bim=plink.get_bim()
    plink_fam=plink.get_fam().astype({'fid':str,'iid':str}).rename(columns={'fid':'FID','iid':'IID','father':'fID', 'mother':'mID','gender':'sex'})
    
    log.info("{} samples ({} males, {} females) loaded from {}".format(plink_fam.shape[0],(plink_fam['sex']==1).sum(),(plink_fam['sex']==2).sum(),args.bfile))
    log.info("{} unphased variants loaded from {}".format(plink_bim.shape[0],args.bfile))


# In[12]:


phased_FID_list=None
phased_IID_list=None
phased_fID_list=None
phased_mID_list=None
phased_sex_list=None

phased_marker_name_list=None
phased_marker_data_list=None

if args.bgl_phased is not None:
    log.info("Loading bgl phased")

    with open(args.bgl_phased,'r') as f:
        line_cnt=0
        while True:
            
            line=f.readline()            
            
            if not line or line_cnt%1000==5:
                sys.stdout.write('\r read %5d markers' % (line_cnt-5))
                sys.stdout.flush()             
                if not line:
                    break
                    
            line_cnt+=1
            line_split=line.strip().split(' ')
            line_type,line_id,line_data=line_split[0],line_split[1],line_split[2:]
            if line_type=='P':
                phased_FID_list1=np.array([line_data[i+0] for i in range(0,len(line_data),2)])
                phased_FID_list2=np.array([line_data[i+1] for i in range(0,len(line_data),2)])
                if np.all(phased_FID_list1==phased_FID_list2):
                    phased_FID_list=phased_FID_list1
                else:
                    raise
            elif line_type=='fID':
                phased_fID_list1=np.array([line_data[i+0] for i in range(0,len(line_data),2)])
                phased_fID_list2=np.array([line_data[i+1] for i in range(0,len(line_data),2)])
                if np.all(phased_fID_list1==phased_fID_list2):
                    phased_fID_list=phased_fID_list1
                else:
                    raise
            elif line_type=='mID':
                phased_mID_list1=np.array([line_data[i+0] for i in range(0,len(line_data),2)])
                phased_mID_list2=np.array([line_data[i+1] for i in range(0,len(line_data),2)])
                if np.all(phased_mID_list1==phased_mID_list2):
                    phased_mID_list=phased_mID_list1
                else:
                    raise      
            elif line_type=='I':        
                phased_IID_list1=np.array([line_data[i+0] for i in range(0,len(line_data),2)])
                phased_IID_list2=np.array([line_data[i+1] for i in range(0,len(line_data),2)])
                if np.all(phased_IID_list1==phased_IID_list2):
                    phased_IID_list=phased_IID_list1
                else:
                    raise   
            elif line_type=='C':
                phased_sex_list1=np.array([line_data[i+0] for i in range(0,len(line_data),2)])
                phased_sex_list2=np.array([line_data[i+1] for i in range(0,len(line_data),2)])
                if np.all(phased_sex_list1==phased_sex_list2):
                    phased_sex_list=np.array(phased_sex_list1).astype(int)
                else:
                    raise  
            elif line_type=='M':
                if phased_marker_name_list is None:
                    phased_marker_name_list=[]
                if phased_marker_data_list is None:
                    phased_marker_data_list=[]                    
                phased_marker_name_list.append(line_id)
                line_data=np.array(line_data)
                phased_marker_data_list.append(line_data)
            else:
                print(line_type)
                raise 
                
    assert phased_FID_list is not None
    assert phased_IID_list is not None
    assert phased_fID_list is not None
    assert phased_mID_list is not None
    assert phased_sex_list is not None
    assert len(phased_marker_name_list)!=0
    assert len(phased_marker_data_list)!=0
    
    log.info("{} phsaed variants loaded from {}".format(len(phased_marker_name_list),args.bgl_phased))
    log.info("{} samples ({} males, {} females) loaded from {}".format(len(phased_IID_list),(np.array(phased_sex_list).astype(int)==1).sum(),(np.array(phased_sex_list).astype(int)==2).sum(),args.bgl_phased))


# In[13]:


pheno=pd.read_csv(args.pheno,header=None,sep='\t',names=['FID','IID','pheno'])
pheno['pheno']=pheno['pheno'].replace(-9,np.nan)

if args.assoc=='linear':
    assert len(pheno['pheno'].unique())>2
else:
    assert np.all(np.isnan(pheno['pheno'])|(pheno['pheno']==1)|(pheno['pheno']==2))
    pheno['pheno']=pheno['pheno']-1
    
log.info("{} pheotype loaded from {}".format(pheno.shape[0],args.pheno))
log.info("Among them, valid: {}, missing: {}".format((~pheno['pheno'].isnull()).sum(),pheno['pheno'].isnull().sum()))
if assoc=='linear':
    log.info("mean={:.3f} std={:.3f} median={:.3f} min={:.3f} max={:.3f}".format(pheno['pheno'].mean(),pheno['pheno'].std(),pheno['pheno'].median(),pheno['pheno'].min(),pheno['pheno'].max()))
else:
    log.info("case: {} / control: {}".format((pheno['pheno']==1).sum(),(pheno['pheno']==0).sum()))    


# # parse multialleic regular exp

# In[14]:


log.info_head("Multialleic expression parsing")

plink_multialleic_dict={}
plink_multialleic_always_dict={}

phased_multialleic_dict={}
phased_multialleic_always_dict={}

for expression in args.multialleic.split(','):
    re_exp=re.compile(expression)
    if plink is not None:
        for marker in plink_bim.index:
            name,allele=(re_exp.search(marker).group('name'),re_exp.search(marker).group('allele')) if re_exp.search(marker) is not None else (None,None)
            if name is not None:
                plink_multialleic_dict[marker]=name
    if phased_marker_name_list is not None:
        for marker in phased_marker_name_list:
            name,allele=(re_exp.search(marker).group('name'),re_exp.search(marker).group('allele')) if re_exp.search(marker) is not None else (None,None)
            if name is not None:
                phased_multialleic_dict[marker]=name  


for expression in args.multialleic_always.split(','):
    re_exp=re.compile(expression)
    if plink is not None:
        for marker in plink_bim.index:
            name,allele=(re_exp.search(marker).group('name'),re_exp.search(marker).group('allele')) if re_exp.search(marker) is not None else (None,None)
            if name is not None:
                plink_multialleic_always_dict[marker]=name
                
    if phased_marker_name_list is not None:
        for marker in phased_marker_name_list:
            name,allele=(re_exp.search(marker).group('name'),re_exp.search(marker).group('allele')) if re_exp.search(marker) is not None else (None,None)
            if name is not None:
                phased_multialleic_always_dict[marker]=name
                
plink_multialleic_df=pd.DataFrame(list(zip(plink_multialleic_dict.keys(),plink_multialleic_dict.values())),columns=['marker','name'])
plink_multialleic_df['from']='plink'
plink_multialleic_df['always']=False
plink_multialleic_always_df=pd.DataFrame(list(zip(plink_multialleic_always_dict.keys(),plink_multialleic_always_dict.values())),columns=['marker','name'])
plink_multialleic_always_df['from']='plink'
plink_multialleic_always_df['always']=True

phased_multialleic_df=pd.DataFrame(list(zip(phased_multialleic_dict.keys(),phased_multialleic_dict.values())),columns=['marker','name'])
phased_multialleic_df['from']='phased'
phased_multialleic_df['always']=False
phased_multialleic_always_df=pd.DataFrame(list(zip(phased_multialleic_always_dict.keys(),phased_multialleic_always_dict.values())),columns=['marker','name'])
phased_multialleic_always_df['from']='phased'
phased_multialleic_always_df['always']=True

multialleic_df_concat=pd.concat([plink_multialleic_df,plink_multialleic_always_df,phased_multialleic_df,phased_multialleic_always_df],sort=False)

multialleic_collapse=set(multialleic_df_concat[multialleic_df_concat['always']==True]['name']).intersection(set(multialleic_df_concat[multialleic_df_concat['always']==False]['name']))

assert len(multialleic_collapse)==0
# Note: duplicated mulltialleic marker in --bfile and --bgl-phased is available. In that case, priority is on --bgl-phased.    
    
log.info("plink, multialleic: {}".format(','.join(multialleic_df_concat[(multialleic_df_concat['from']=='plink')&(multialleic_df_concat['always']==False)]['name'].unique())))
log.info("plink, multialleic always: {}".format(','.join(multialleic_df_concat[(multialleic_df_concat['from']=='plink')&(multialleic_df_concat['always']==True)]['name'].unique())))
log.info("phased, multialleic: {}".format(','.join(multialleic_df_concat[(multialleic_df_concat['from']=='phased')&(multialleic_df_concat['always']==False)]['name'].unique())))
log.info("phased, multialleic always: {}".format(','.join(multialleic_df_concat[(multialleic_df_concat['from']=='phased')&(multialleic_df_concat['always']==True)]['name'].unique())))    


# # parse optional input files

# In[15]:


if args.covar is None:
    covar=fam.iloc[:,:2]
else:
    covar=pd.read_csv(args.covar,sep='\t')
    covar.columns=['FID','IID']+covar.columns[2:].tolist()
    covar=covar.astype({'FID':str,'IID':str})
    
    covar.iloc[:,2:]=covar.iloc[:,2:].astype(float)
    covar.iloc[:,2:]=covar.iloc[:,2:].replace(-9,np.nan)
    
    log.info("{} covariates loaded from {}".format(len(covar.columns[2:]),args.covar))


# In[16]:


if args.condition_list is None:
    condition_list=[]
else:
    with open(args.condition_list,'r') as f:
        condition_list=f.read().strip().split('\n')
        if condition_list[0]=='':
            condition_list=[]
            log.warning("Empty --condition-list {}".format(args.condition_list))
        else:
            #condition_list.append()
            log.info("{} conditions loaded from --condition-list {}".format(len(condition_list),args.condition_list))
            if len(np.unique(condition_list))!=len(condition_list):
                condition_list=np.unique(condition_list).tolist()
                log.info("After removing duplicated conditions, {} conditions remains".format(len(condition_list)))
            for condition1 in condition_list:
                if condition1 not in multialleic_df_concat['name'].values:
                    for condition2 in condition_list:
                        if condition1 in multialleic_df_concat[multialleic_df_concat['name']==condition2]['marker'].values:
                            condition_list.remove(condition1)
                            log.info("Removed bialleic condition({}) with correponding multialleic condition({})".format(condition1,condition2))
            log.info("Finally {} conditions remains".format(len(condition_list)))
            log.info('*********\n '+', '.join(condition_list)+'\n*********')


# # check idx integrity

# In[17]:


log.info_head("Input integrity check")
if plink_fam is not None and phased_FID_list is not None:
    assert np.all(plink_fam['FID']==phased_FID_list)
    assert np.all(plink_fam['IID']==phased_IID_list)
    assert np.all(plink_fam['fID']==phased_fID_list)
    assert np.all(plink_fam['mID']==phased_mID_list)
    assert np.all(plink_fam['sex']==phased_sex_list)
    log.info("Passed individual integrity check (Individuals from --bfile is the same as individuals from --bgl-phased)")

assert np.all(covar['FID']==(plink_fam['FID'] if plink_fam is not None else phased_FID_list))
assert np.all(covar['IID']==(plink_fam['IID'] if plink_fam is not None else phased_IID_list))
log.info("Passed individual integrity check (Individuals from --bfile or --bgl-phased is the same as individuals from --covar)")

diff=set(condition_list)
if phased_marker_name_list is not None:
    diff=diff.difference(phased_marker_name_list)
if plink_bim is not None:
    diff=diff.difference(plink_bim.index)
diff=diff.difference(multialleic_df_concat['name'])
assert len(diff)==0
log.info("Passed condition integrity check (All variants in --condition-list are identified from loaded variants)")


# In[18]:


log.info_head("Converting condtion to covariate")


# In[19]:


def plink_get_dosage(marker,keep_allele_order=True,repeat=1):
    dosage=plink.get_geno_marker(marker).astype(float)
    dosage[dosage==-1]=np.nan
    if keep_allele_order or ((dosage==0).sum()>(dosage==2).sum()):
        a1=plink_bim.loc[marker]['a1']
    else:
        a1=plink_bim.loc[marker]['a2']
        dosage=2-dosage
    return a1,np.repeat(dosage,repeat)

def phased_get_dosage(marker,a1=None):
    phased_marker_idx=phased_marker_name_list.index(marker)
    phased_marker_data=phased_marker_data_list[phased_marker_idx]
    phased_marker_data_unique=np.unique(phased_marker_data)
    if len(phased_marker_data_unique)>2:
        raise NotImplementedError
    if a1 is not None:
        if a1==phased_marker_data_unique[1]:
            a2=phased_marker_data_unique[0]     
        elif a1==phased_marker_data_unique[0]:
            a2=phased_marker_data_unique[1]       
        else:
            raise NotImplementedError
    elif (phased_marker_data==phased_marker_data_unique[0]).sum()>(phased_marker_data==phased_marker_data_unique[1]).sum():
        a1=phased_marker_data_unique[1]
        a2=phased_marker_data_unique[0]
    else:
        a1=phased_marker_data_unique[0]
        a2=phased_marker_data_unique[1]
    phased_marker_data=np.where(phased_marker_data==a1, 1, phased_marker_data)
    phased_marker_data=np.where(phased_marker_data==a2, 0, phased_marker_data)
    return a1, phased_marker_data.astype(float)


# In[20]:


def find_trivial_index(array2d):
    array2d_sumcol=array2d.sum(axis=1)
    array2d_sumrow=array2d.sum(axis=0)   
    
    array2d_sumrow_argmax=np.argmax(array2d_sumrow)
    
    if array2d_sumrow.shape[0]>1:
        return array2d_sumrow_argmax
    else:
        return None


# In[21]:


covar_phased=covar.loc[covar.index.repeat(2)].reset_index().drop(columns='index')


# In[22]:


pheno_phased=pheno.loc[pheno.index.repeat(2)].reset_index().drop(columns='index')#['pheno']


# In[23]:


for condition in condition_list:
    if condition in multialleic_df_concat['name'].values:
        if condition in multialleic_df_concat[(multialleic_df_concat['from']=='phased')]['name'].values:
            bialleic_marker_list=multialleic_df_concat[(multialleic_df_concat['from']=='phased')&(multialleic_df_concat['name']==condition)]['marker'].values
            bialleic_marker_info_list=[phased_get_dosage(bialleic_marker) for bialleic_marker in bialleic_marker_list]
        elif condition in multialleic_df_concat[(multialleic_df_concat['from']=='plink')]['name'].values:
            bialleic_marker_list=multialleic_df_concat[(multialleic_df_concat['from']=='plink')&(multialleic_df_concat['name']==condition)]['marker'].values
            bialleic_marker_info_list=[plink_get_dosage(bialleic_marker,repeat=2) for bialleic_marker in bialleic_marker_list]  
        else:
            raise
        if len(np.unique(bialleic_marker_list))!=len(bialleic_marker_list):
            raise
            
        bialleic_marker_dosage=np.array([dosage for allele,dosage in bialleic_marker_info_list]).transpose()
        trivial_index=find_trivial_index(bialleic_marker_dosage)
        bialleic_marker_list_cut=bialleic_marker_list if trivial_index is None else np.delete(bialleic_marker_list, trivial_index)
        bialleic_marker_dosage_cut=bialleic_marker_dosage if trivial_index is None else np.delete(bialleic_marker_dosage, trivial_index,axis=1)        

        for bialleic_marker_idx,bialleic_marker in enumerate(bialleic_marker_list_cut):
            covar_phased[bialleic_marker]=bialleic_marker_dosage_cut[:,bialleic_marker_idx]
            
        log.info("{} bialleic marker(s) from mulitalleic marker specifier({}) were added.".format(len(bialleic_marker_list_cut),condition))
        
        if trivial_index is not None:
            log.info("==> To avoid coliearity, {} removed from {}".format(bialleic_marker_list[trivial_index] ,', '.join(bialleic_marker_list)))
    elif phased_marker_name_list is not None and condition in phased_marker_name_list:
        allele,dosage=phased_get_dosage(condition)
        covar_phased[condition]=dosage
        log.info("1 bialleic marker {} was added from --bgl-phased".format(condition))
    elif plink_bim is not None and condition in plink_bim.index:
        allele,dosage=plink_get_dosage(condition,repeat=2)
        covar_phased[condition]=dosage        
        log.info("1 bialleic marker {} was added from --bfile".format(condition))
    else:
        raise NotImplementedError


# # Run regression

# In[24]:


log.info_head("Regression")
log.info("Start time: "+time.strftime('%c', time.localtime(time.time())))


# In[25]:


test_marker_list=sorted(np.unique(plink_bim.index.tolist()+phased_marker_name_list+multialleic_df_concat['name'].tolist()).tolist())

test_marker_list=pd.Index(test_marker_list)
test_marker_list=test_marker_list.difference(multialleic_df_concat[multialleic_df_concat['always']==True]['marker'])

test_marker_list=list(test_marker_list)#np.random.shuffle(test_marker_list_temp)


# In[26]:


x_data_intercept=np.array([np.ones(2*plink_fam.shape[0])]).transpose()    
x_data_covariate=covar_phased.iloc[:,2:].values

x_data_null=np.concatenate([x_data_intercept,x_data_covariate],axis=1)
x_data_null_names=['const']+covar_phased.columns[2:].tolist()

y_data=pheno_phased['pheno'].values


# In[87]:


reduce_1d=lambda x: np.mean(y_data.reshape(-1,2),axis=1)
reduce_2d=lambda x: np.mean(x.reshape(int(x.shape[0]/2),-1,x.shape[1]),axis=1)


# In[88]:


assoc_result_list=[]
assoc_result_list_keys=['marker_name','P','nobs','coef','std','Z','chisq','df','term','A1','multi_allele','note']
def assoc_result_record(marker_name='',P=np.nan,nobs=np.nan,coef=np.nan,std=np.nan,Z=np.nan,chisq=np.nan,df=np.nan,term=np.nan,A1=np.nan,multi_allele=np.nan,note=''):
    assoc_result_list.append({'marker_name':marker_name,'P':P,'nobs':nobs,'coef':coef,'std':std,'Z':Z,'chisq':chisq,'df':df,'term':term,'A1':A1,'multi_allele':multi_allele,'note':note})


# In[99]:


family=(sm.families.Gaussian() if assoc=='linear' else sm.families.Binomial())

for marker_idx,marker in enumerate(test_marker_list):
    if marker_idx%5000==0:
        log.info("Time: {} - {:.3f} %".format(time.strftime('%c', time.localtime(time.time())),marker_idx/len(test_marker_list)*100))
    try:
        if phased_marker_name_list is not None and marker in phased_marker_name_list:
            allele,dosage=phased_get_dosage(marker)        
            dosage=np.expand_dims(dosage,axis=1)
            x_data_full=np.concatenate([x_data_null,dosage],axis=1)
            x_data_full_names=x_data_null_names+['THIS']


            model=sm.GLM(y_data, x_data_full, family=family,missing='drop')
            model_result=model.fit()

            for model_result_idx in range(len(model_result.params)):
                assoc_result_record(marker_name=marker,
                                    P=model_result.pvalues[model_result_idx],
                                    coef=model_result.params[model_result_idx],
                                    std=model_result.bse[model_result_idx],
                                    Z=model_result.tvalues[model_result_idx],
                                    term=x_data_full_names[model_result_idx],
                                    A1=allele if x_data_full_names[model_result_idx]=='THIS' else np.nan,
                                    nobs=model_result.nobs,
                                    note='phased bialleic')        

        elif marker in multialleic_df_concat[multialleic_df_concat['from']=='phased']['name'].values:
            bialleic_marker_list=multialleic_df_concat[(multialleic_df_concat['from']=='phased')&(multialleic_df_concat['name']==marker)]['marker'].values
            if len(np.unique(bialleic_marker_list))!=len(bialleic_marker_list):
                raise
            bialleic_marker_info_list=[phased_get_dosage(bialleic_marker) for bialleic_marker in bialleic_marker_list]

            bialleic_marker_dosage=np.array([dosage for allele,dosage in bialleic_marker_info_list]).transpose()
            trivial_index=find_trivial_index(bialleic_marker_dosage)
            bialleic_marker_list_cut=bialleic_marker_list if trivial_index is None else np.delete(bialleic_marker_list, trivial_index)
            bialleic_marker_dosage_cut=bialleic_marker_dosage if trivial_index is None else np.delete(bialleic_marker_dosage, trivial_index,axis=1)

            x_data_full=np.concatenate([x_data_null,bialleic_marker_dosage_cut],axis=1)       


            common_idx=(~np.isnan(y_data))&(np.isnan(x_data_null).sum(axis=1)==0)&(np.isnan(x_data_full).sum(axis=1)==0)

            model_null=sm.GLM(y_data[common_idx], x_data_null[common_idx], family=family,missing='raise')
            model_null_result=model_null.fit()

            model_full=sm.GLM(y_data[common_idx], x_data_full[common_idx], family=family,missing='raise')
            model_full_result=model_full.fit()        

            chisq_diff=2*(model_full_result.llf-model_null_result.llf)
            df_diff=model_full_result.df_model-model_null_result.df_model
            p_value=chi2.sf(chisq_diff,df_diff)          

            assoc_result_record(marker_name=marker,
                                P=p_value,
                                chisq=chisq_diff,
                                df=df_diff,
                                multi_allele=','.join(bialleic_marker_list),
                                nobs=np.sum(common_idx),
                                note='phased multialleic')          

        elif plink_bim is not None and marker in plink_bim.index:
            allele,dosage=plink_get_dosage(marker,repeat=2)
            dosage=np.expand_dims(dosage,axis=1)
            x_data_full=np.concatenate([x_data_null,dosage],axis=1)
            x_data_full_names=x_data_null_names+['THIS']
            
            y_data_reduce=reduce_1d(y_data)
            x_data_full_reduce=reduce_2d(x_data_full)
            
            model=sm.GLM(y_data_reduce, x_data_full_reduce, family=family,missing='drop')
            model_result=model.fit()

            for model_result_idx in range(len(model_result.params)):
                assoc_result_record(marker_name=marker,
                                    P=model_result.pvalues[model_result_idx],
                                    coef=model_result.params[model_result_idx],
                                    std=model_result.bse[model_result_idx],
                                    Z=model_result.tvalues[model_result_idx],
                                    term=x_data_full_names[model_result_idx],
                                    A1=allele if x_data_full_names[model_result_idx]=='THIS' else np.nan,
                                    nobs=model_result.nobs,
                                    note='plink bialleic')

        elif marker in multialleic_df_concat[multialleic_df_concat['from']=='plink']['name'].values:
            bialleic_marker_list=multialleic_df_concat[(multialleic_df_concat['from']=='plink')&(multialleic_df_concat['name']==marker)]['marker'].values
            if len(np.unique(bialleic_marker_list))!=len(bialleic_marker_list):
                raise
            bialleic_marker_info_list=[plink_get_dosage(bialleic_marker,repeat=2) for bialleic_marker in bialleic_marker_list]

            bialleic_marker_dosage=np.array([dosage for allele,dosage in bialleic_marker_info_list]).transpose()
            trivial_index=find_trivial_index(bialleic_marker_dosage)
            bialleic_marker_list_cut=bialleic_marker_list if trivial_index is None else np.delete(bialleic_marker_list, trivial_index)
            bialleic_marker_dosage_cut=bialleic_marker_dosage if trivial_index is None else np.delete(bialleic_marker_dosage, trivial_index,axis=1)

            x_data_full=np.concatenate([x_data_null,bialleic_marker_dosage_cut],axis=1)       
            
            x_data_full_reduce=reduce_2d(x_data_full)
            x_data_null_reduce=reduce_2d(x_data_null)
            
            y_data_reduce=reduce_1d(y_data)
            
            common_idx=(~np.isnan(y_data_reduce))&(np.isnan(x_data_null_reduce).sum(axis=1)==0)&(np.isnan(x_data_full_reduce).sum(axis=1)==0)

            model_null=sm.GLM(y_data_reduce[common_idx], x_data_null_reduce[common_idx], family=family,missing='raise')
            model_null_result=model_null.fit()

            model_full=sm.GLM(y_data_reduce[common_idx], x_data_full_reduce[common_idx], family=family,missing='raise')
            model_full_result=model_full.fit()        

            chisq_diff=2*(model_full_result.llf-model_null_result.llf)
            df_diff=model_full_result.df_model-model_null_result.df_model
            p_value=chi2.sf(chisq_diff,df_diff)          

            assoc_result_record(marker_name=marker,
                                P=p_value,
                                chisq=chisq_diff,
                                df=df_diff,
                                multi_allele=','.join(bialleic_marker_list),
                                nobs=np.sum(common_idx),
                                note='plink multialleic') 
    
    except sm.tools.sm_exceptions.PerfectSeparationError as e:
        log.waring("{} PerfectSeparationError".format(marker))
        assoc_result_record(marker_name=marker,note='PerfectSeparationError')
    else:
        pass


# In[100]:


assoc_result_df_verbose=pd.DataFrame(assoc_result_list)[assoc_result_list_keys]
assoc_result_df_concise=assoc_result_df_verbose[(assoc_result_df_verbose['term'].isnull())|(assoc_result_df_verbose['term']=='THIS')]


# In[101]:


assoc_result_df_verbose.to_csv(out+'_verbose.tsv',sep='\t',index=None)
assoc_result_df_concise.to_csv(out+'.tsv',sep='\t',index=None)


# In[102]:


log.info("End time: "+time.strftime('%c', time.localtime(time.time())))


# In[103]:


assoc_result_df_concise.head()

