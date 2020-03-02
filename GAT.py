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


# In[2]:


def dir_path(path):
    if os.path.exists(os.path.dirname(path)):
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
        
parser = argparse.ArgumentParser(description='Multialleic association test')

#required mode
parser.add_argument('--assoc',choices=['linear','logistic'],required=True)

#required output file
parser.add_argument('--out', type=dir_path,required=True,help='output file prefix. (prefix.log, prefix.assoc will be generated)')

#required input files
parser.add_argument('--bgl-phased', type=file_path,help='bgl-phased (See Beagle 5.1 documentation)')
parser.add_argument('--bfile', type=bfile_path,help='plink binary format')
parser.add_argument('--multialleic',type=str,help='regular expression for specifying multiple alleic marker (comma delimiter)')
parser.add_argument('--multialleic-always',type=str,help='regular expression for specifying multiple alleic marker (comma delimiter)')

parser.add_argument('--pheno', type=file_path,required=True,help='format is the same as plink. Tab-delimited file without header of which the first and second columns is family and within-family IDs respectively and the third columns are pheotype')

#optional
parser.add_argument('--covar', type=file_path,help='format is the same as plink')
parser.add_argument('--condition-list',type=file_path,help='format is the same as plink')


# In[3]:


debug=True
#debug=True

if debug:
    arg_split='--assoc linear --out data/out_assoc/ALP/step_01 --bgl-phased data/genotype/4_merge/KCHIP_HLA_AA_SNP.bgl.phased --bfile data/genotype/4_merge/KCHIP_HLA_SNP_1000G --multialleic (?P<name>HLA_[0-9A-Z]*)\*(?P<allele>[0-9:]*) --multialleic-always (?P<name>AA_[A-Z0-9]*_[0-9]*_[0-9]*_exon[0-9]*)_*(?P<allele>[A-Z]*) --pheno data/out_pheno/ALP.phe --covar data/out_assoc/ALP/step_01.covar --condition-list data/out_assoc/ALP/step_01.cond'.split(' ')
    args=parser.parse_args(arg_split)
else:
    args=parser.parse_args()
    
if args.bfile is None and args.bgl_phased is None:
    raise argparse.ArgumentTypeError("either --bfile or --bgl-phased parameter is needed")    


# In[4]:


log = logging.getLogger('logger')
log.setLevel(logging.DEBUG)

log_file_path=args.out+'.log'
fileHandler = logging.FileHandler(log_file_path)
streamHandler = logging.StreamHandler()

formatter = logging.Formatter('%(message)s')
fileHandler.setFormatter(formatter)
streamHandler.setFormatter(formatter)

log.addHandler(fileHandler)
log.addHandler(streamHandler)


# In[5]:


log.info_head=lambda x: log.info('*'*int((100-len(x))/2)+x+'*'*int((100-len(x))/2))


# In[6]:


log.info_head("*********************************")
log.info("* Generic Association Tool")
log.info("* Description: Generic module for testing bialleic/multialleic phased/unphased markers")
log.info("* version 1.0")
log.info("* (C) 2020-, Seoul National University")
log.info("* Please report bugs to: Chanwoo Kim <ch6845@snu.ac.kr>")
log.info("* https://github.com/ch6845/Generic_Association_Tool")
log.info_head("*********************************")


# In[7]:


log.info("Start time: "+time.strftime('%c', time.localtime(time.time())))


# In[8]:


log.info('Working directory: '+os.getcwd())
log.info('Hostname: '+socket.gethostname())


# In[9]:


log.info('Parameters\n'+'\n'.join(['--{} {}'.format(key,value) for key,value in vars(args).items()]))


# In[10]:


assoc=args.assoc
out=args.out


# In[11]:


log.info_head("Data Loading")


# # parse input files

# In[12]:


if args.bfile is not None:
    plink=PyPlink(args.bfile)
    plink_bim=plink.get_bim()
    plink_fam=plink.get_fam().astype({'fid':str,'iid':str}).rename(columns={'fid':'FID','iid':'IID','father':'fID', 'mother':'mID','gender':'sex'})
    
    log.info("{} samples ({} males, {} females) loaded from {}".format(plink_fam.shape[0],(plink_fam['sex']==1).sum(),(plink_fam['sex']==2).sum(),args.bfile))
    log.info("{} unphased variants loaded from {}".format(plink_bim.shape[0],args.bfile))
else:
    plink=None
    plink_bim=None
    plink_fam=None


# In[13]:


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
    #aa_marker_name_list_aaonly=pd.Series(aa_marker_name_list)[pd.Series(aa_marker_name_list).str.slice(stop=3)=='AA_'].values            


# In[14]:


pheno=pd.read_csv(args.pheno,header=None,sep='\t',names=['FID','IID','pheno'])
pheno['pheno']=pheno['pheno'].replace(-9,np.nan)


# In[15]:


if args.assoc=='linear':
    assert len(pheno['pheno'].unique())>2
else:
    assert np.all(np.isnan(pheno['pheno'])|(pheno['pheno']==1)|(pheno['pheno']==2))
    pheno['pheno']=pheno['pheno']-1


# In[16]:


log.info("{} pheotype loaded from {}".format(pheno.shape[0],args.pheno))
log.info("Among them, valid: {}, missing: {}".format((~pheno['pheno'].isnull()).sum(),pheno['pheno'].isnull().sum()))
if assoc=='linear':
    log.info("mean={:.3f} std={:.3f} median={:.3f} min={:.3f} max={:.3f}".format(pheno['pheno'].mean(),pheno['pheno'].std(),pheno['pheno'].median(),pheno['pheno'].min(),pheno['pheno'].max()))
else:
    log.info("case: {} / control: {}".format((pheno['pheno']==1).sum(),(pheno['pheno']==0).sum()))


# # parse multialleic regular exp

# In[17]:


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
plink_multialleic_df['from']='plink_multialleic'
plink_multialleic_always_df=pd.DataFrame(list(zip(plink_multialleic_always_dict.keys(),plink_multialleic_always_dict.values())),columns=['marker','name'])
plink_multialleic_always_df['from']='plink_multialleic_always'

phased_multialleic_df=pd.DataFrame(list(zip(phased_multialleic_dict.keys(),phased_multialleic_dict.values())),columns=['marker','name'])
phased_multialleic_df['from']='phased_multialleic'
phased_multialleic_always_df=pd.DataFrame(list(zip(phased_multialleic_always_dict.keys(),phased_multialleic_always_dict.values())),columns=['marker','name'])
phased_multialleic_always_df['from']='phased_multialleic_always'


multialleic_df_concat=pd.concat([plink_multialleic_df,plink_multialleic_always_df,phased_multialleic_df,phased_multialleic_always_df],sort=False)


multialleic_df_concat_notalways=multialleic_df_concat[(multialleic_df_concat['from']=='plink_multialleic') | (multialleic_df_concat['from']=='phased_multialleic')]
multialleic_df_concat_always=multialleic_df_concat[(multialleic_df_concat['from']=='plink_multialleic_always') | (multialleic_df_concat['from']=='phased_multialleic_always')]
multialleic_collapse=set(multialleic_df_concat_always['name']).intersection(set(multialleic_df_concat_notalways['name']))

if len(multialleic_collapse)!=0:
    raise


# In[18]:


log.info("plink, multialleic: {}".format(','.join(plink_multialleic_df['name'].unique())))
log.info("plink, multialleic always: {}".format(','.join(plink_multialleic_always_df['name'].unique())))
log.info("phased, multialleic: {}".format(','.join(phased_multialleic_df['name'].unique())))
log.info("phased, multialleic always: {}".format(','.join(phased_multialleic_always_df['name'].unique())))

multialleic_df_concat['from']


# In[ ]:





# # parse optional input files

# In[19]:


if args.covar is None:
    covar=fam.iloc[:,:2]
else:
    covar=pd.read_csv(args.covar,sep='\t')
    covar.columns=['FID','IID']+covar.columns[2:].tolist()
    covar=covar.astype({'FID':str,'IID':str})
    
    covar.iloc[:,2:]=covar.iloc[:,2:].astype(float)
    covar.iloc[:,2:]=covar.iloc[:,2:].replace(-9,np.nan)
    
    log.info("{} covariates loaded from {}".format(len(covar.columns[2:]),args.covar))


# In[20]:


if args.condition_list is None:
    condition_list=[]
else:
    with open(args.condition_list,'r') as f:
        condition_list=f.read().strip().split('\n')
        if condition_list[0]=='':
            condition_list=[]
            log.warning("Empty --condition-list {}".format(args.condition_list))
        else:
            condition_list.append('AA_A_9_30018537_exon2_F')
            log.info("{} conditions loaded from --condition-list {}".format(len(condition_list),args.condition_list))
            if len(np.unique(condition_list))!=len(condition_list):
                condition_list=np.unique(condition_list).tolist()
                log.info("After removing duplicated conditions, {} conditions remains".format(len(condition_list)))
            #log.info(', '.join(condition_list))
            for condition1 in condition_list:
                if condition1 not in multialleic_df_concat['name'].values:
                    for condition2 in condition_list:
                        #print(condition2,multialleic_df_concat[multialleic_df_concat['name']==condition2]['marker'].values)
                        if condition1 in multialleic_df_concat[multialleic_df_concat['name']==condition2]['marker'].values:
                            condition_list.remove(condition1)
                            log.info("Removed bialleic condition({}) with correponding multialleic condition({})".format(condition1,condition2))
            log.info("Finally {} conditions remains".format(len(condition_list)))
            log.info('*********\n '+', '.join(condition_list)+'\n*********')


# # check idx integrity

# In[21]:


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


# # Run regression

# In[22]:


log.info_head("Converting condtion to covariate")


# def marker_data_to_onehot(aa_marker_data):
#     
#     aa_marker_data_unique=np.unique(aa_marker_data)
#     aa_marker_data_unique_nonan=aa_marker_data_unique[aa_marker_data_unique!='nan'].tolist()
#     aa_marker_data_unique_nan=aa_marker_data_unique_nonan+['nan']
#     
#     aa_marker_data_int_nan=list(map(lambda x: aa_marker_data_unique_nan.index(x),aa_marker_data))
#     
#     aa_marker_data_onehot_nan=np.zeros((len(aa_marker_data),len(aa_marker_data_unique_nan)))
#     aa_marker_data_onehot_nan[np.arange(len(aa_marker_data)),aa_marker_data_int_nan]=1
#     aa_marker_data_onehot_nonan=aa_marker_data_onehot_nan[:,:-1]
#     
#     return aa_marker_data_unique_nonan,aa_marker_data_onehot_nonan
# 
# def prepare_onehot(aa_marker_data_unique_nonan,aa_marker_data_onehot_nonan,set_nan=True,cut_mostfrequent=True):
#     assert len(aa_marker_data_unique_nonan)==aa_marker_data_onehot_nonan.shape[1]
#     assert np.isnan(aa_marker_data_onehot_nonan).sum()==0
#     assert (np.nan not in aa_marker_data_unique_nonan) and ('nan' not in aa_marker_data_unique_nonan)
#     
#     aa_marker_data_onehot_nonan_sumrow=np.sum(aa_marker_data_onehot_nonan,axis=1)
#     aa_marker_data_onehot_nonan_sumcol=np.sum(aa_marker_data_onehot_nonan,axis=0)
#     #print(aa_marker_data_onehot_nonan_sumcol,np.argmax(aa_marker_data_onehot_nonan_sumcol))
#     #print(aa_marker_data_onehot_nonan_sumrow.shape,aa_marker_data_onehot_nonan_sumcol.shape)
#     #print(np.argmax(aa_marker_data_onehot_nonan_sumcol))
#     if set_nan:
#         aa_marker_data_onehot_nonan[aa_marker_data_onehot_nonan_sumrow!=1,:]=np.nan
#     if cut_mostfrequent:
#         aa_marker_data_unique_nonan=np.delete(aa_marker_data_unique_nonan, np.argmax(aa_marker_data_onehot_nonan_sumcol))
#         aa_marker_data_onehot_nonan=np.delete(aa_marker_data_onehot_nonan, np.argmax(aa_marker_data_onehot_nonan_sumcol), axis=1)
#     return aa_marker_data_unique_nonan,aa_marker_data_onehot_nonan

# In[23]:


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


# In[34]:


covar_phased=covar.loc[covar.index.repeat(2)].reset_index().drop(columns='index')


# In[35]:


for condition in condition_list:
    if condition in multialleic_df_concat['name'].values:
        if condition in phased_multialleic_df['name'].values:
            bialleic_marker_list=np.array(phased_multialleic_df['marker'][phased_multialleic_df['name']==condition].values)
            bialleic_marker_info_list=[phased_get_dosage(bialleic_marker) for bialleic_marker in bialleic_marker_list]
        elif condition in phased_multialleic_always_df['name'].values:
            bialleic_marker_list=np.array(phased_multialleic_always_df['marker'][phased_multialleic_always_df['name']==condition].values)
            bialleic_marker_info_list=[phased_get_dosage(bialleic_marker) for bialleic_marker in bialleic_marker_list]        
            
        elif condition in plink_multialleic_df['name'].values:   
            bialleic_marker_list=np.array(plink_multialleic_df['marker'][plink_multialleic_df['name']==condition].values)
            bialleic_marker_info_list=[plink_get_dosage(bialleic_marker,repeat=2) for bialleic_marker in bialleic_marker_list]   
        elif condition in plink_multialleic_always_df['name'].values:   
            bialleic_marker_list=np.array(plink_multialleic_always_df['marker'][plink_multialleic_always_df['name']==condition].values)
            bialleic_marker_info_list=[plink_get_dosage(bialleic_marker,repeat=2) for bialleic_marker in bialleic_marker_list]               
        if len(np.unique(bialleic_marker_list))!=len(bialleic_marker_list):
            raise        
        bialleic_marker_allele=[allele for allele,dosage in bialleic_marker_info_list]
        bialleic_marker_dosage=np.array([dosage for allele,dosage in bialleic_marker_info_list]).transpose()
        bialleic_marker_dosage_sumcol=bialleic_marker_dosage.sum(axis=1)
        bialleic_marker_dosage_sumrow=bialleic_marker_dosage.sum(axis=0)
        
        if bialleic_marker_dosage_sumrow.shape[0]>1:
            bialleic_marker_list_cut=np.delete(bialleic_marker_list, np.argmax(bialleic_marker_dosage_sumrow))
            bialleic_marker_dosage_cut=np.delete(bialleic_marker_dosage, np.argmax(bialleic_marker_dosage_sumrow),axis=1)
        else:
            bialleic_marker_list_cut=bialleic_marker_list
            bialleic_marker_dosage_cut=bialleic_marker_dosage
            
        for bialleic_marker_idx,bialleic_marker in enumerate(bialleic_marker_list_cut):
            covar_phased[bialleic_marker]=bialleic_marker_dosage_cut[:,bialleic_marker_idx]
            
        log.info("{} bialleic markers from mulitalleic marker({}) were added.".format(len(bialleic_marker_list_cut),condition))
        
        if bialleic_marker_dosage_sumrow.shape[0]>1:
            log.info("==> To avoid coliearity, {} removed from {}".format(bialleic_marker_list[np.argmax(bialleic_marker_dosage_sumrow)] ,', '.join(bialleic_marker_list)))
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


# In[36]:


test_marker_list=np.unique(plink_bim.index.tolist()+phased_marker_name_list+multialleic_df_concat['name'].unique().tolist()).tolist()


# In[59]:


test_marker_list=pd.Index(test_marker_list)
test_marker_list=test_marker_list.difference(multialleic_df_concat_always['marker'])


# In[ ]:





# In[ ]:


if plink is not None:
    
    plink_bialleic_list=plink_bim.index
    plink_bialleic_list=plink_bialleic_list.difference(plink_multialleic_dict.keys())
    plink_bialleic_list=plink_bialleic_list.difference(plink_multialleic_always_dict.keys())    
    
    for marker_idx,marker in enumerate(plink_bialleic_list):
        if marker_idx%10==0:
            sys.stdout.write('\r{:.2f}%'.format(100*marker_idx/len(plink_bialleic_list)))
            sys.stdout.flush()           

        
        dosage=plink.get_geno_marker(marker).astype(float)
        dosage[dosage==-1]=np.nan
                
        x_data_intercept=np.array([np.ones(plink_fam.shape[0])]).transpose()    
        x_data_dosage=plink.get_geno_marker(marker).astype(float);x_data_dosage[x_data_dosage==-1]=np.nan
        x_data_dosage=np.expand_dims(x_data_dosage,axis=1)
        x_data_covariate=covar.iloc[:,2:].values

        x_data=np.concatenate([x_data_intercept,x_data_covariate,x_data_dosage],axis=1)#[~x_y_data_nan]
        x_data_names=np.array(['const']+covar.iloc[:,2:].columns.values.tolist()+['THIS'])
        #x_data_null=np.concatenate([x_data_intercept,x_data_covariate],axis=1)#[~x_y_data_nan]
        y_data=pheno['pheno']
        
        family=(sm.families.Gaussian() if assoc=='linear' else sm.families.Binomial())
        model=sm.GLM(y_data,x_data, family=family,missing='drop')
        #model_null=sm.GLM(y_data,x_data_null, family=family,missing='drop')
        model_result=model.fit()
        
        
        for model_result_idx in range(len(model_result.params)):
            
            assoc_result_record(marker_name=marker,
                                P=model_result.pvalues.iloc[model_result_idx],
                                coef=model_result.params.iloc[model_result_idx],
                                std=model_result.bse.iloc[model_result_idx],
                                Z=model_result.tvalues.iloc[model_result_idx],
                                term=x_data_names[model_result_idx],
                                nobs=model_result.nobs,
                                note='plink bialleic'
                               
                               )
        #model_result_null=model_null.fit()            


# In[ ]:


assert len(aa_marker_data_unique_nonan)==aa_marker_data_onehot_nonan.shape[1]
assert np.isnan(aa_marker_data_onehot_nonan).sum()==0
assert (np.nan not in aa_marker_data_unique_nonan) and ('nan' not in aa_marker_data_unique_nonan)

aa_marker_data_onehot_nonan_sumrow=np.sum(aa_marker_data_onehot_nonan,axis=1)
aa_marker_data_onehot_nonan_sumcol=np.sum(aa_marker_data_onehot_nonan,axis=0)
#print(aa_marker_data_onehot_nonan_sumcol,np.argmax(aa_marker_data_onehot_nonan_sumcol))
#print(aa_marker_data_onehot_nonan_sumrow.shape,aa_marker_data_onehot_nonan_sumcol.shape)
#print(np.argmax(aa_marker_data_onehot_nonan_sumcol))
if set_nan:
    aa_marker_data_onehot_nonan[aa_marker_data_onehot_nonan_sumrow!=1,:]=np.nan
if cut_mostfrequent:
    aa_marker_data_unique_nonan=np.delete(aa_marker_data_unique_nonan, np.argmax(aa_marker_data_onehot_nonan_sumcol))
    aa_marker_data_onehot_nonan=np.delete(aa_marker_data_onehot_nonan, np.argmax(aa_marker_data_onehot_nonan_sumcol), axis=1)


# In[ ]:


for marker in plink_bialleic_list:
    name,allele=(re_exp.search(marker).group('name'),re_exp.search(marker).group('allele')) if re_exp.search(marker) is not None else (None,None)
    if name is not None and name.split('*')[0].split('_')[1]==HLA_gene:
        marker_dosage_list=marker_dosage_list_dict.get(name,[])
        marker_dosage_list.append({"marker":marker,"dosage":plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker)})
        marker_dosage_list_dict[name]=marker_dosage_list


# In[215]:


#phased_get_dosage('rs9380355','A')


# In[227]:


assoc_result_list=[]

def assoc_result_record(marker_name='',P=np.nan,nobs=np.nan,coef=np.nan,std=np.nan,Z=np.nan,chisq=np.nan,df=np.nan,term=np.nan,note=''):
    assoc_result_list.append({'marker_name':marker_name,'P':P,'nobs':nobs,'coef':coef,'std':std,'Z':Z,'chisq':chisq,'df':df,'term':term,'note':note})


# In[241]:


if plink is not None:
    
    plink_bialleic_list=plink_bim.index
    plink_bialleic_list=plink_bialleic_list.difference(plink_multialleic_dict.keys())
    plink_bialleic_list=plink_bialleic_list.difference(plink_multialleic_always_dict.keys())    
    
    for marker_idx,marker in enumerate(plink_bialleic_list):
        if marker_idx%10==0:
            sys.stdout.write('\r{:.2f}%'.format(100*marker_idx/len(plink_bialleic_list)))
            sys.stdout.flush()           

        
        dosage=plink.get_geno_marker(marker).astype(float)
        dosage[dosage==-1]=np.nan
                
        x_data_intercept=np.array([np.ones(plink_fam.shape[0])]).transpose()    
        x_data_dosage=plink.get_geno_marker(marker).astype(float);x_data_dosage[x_data_dosage==-1]=np.nan
        x_data_dosage=np.expand_dims(x_data_dosage,axis=1)
        x_data_covariate=covar.iloc[:,2:].values

        x_data=np.concatenate([x_data_intercept,x_data_covariate,x_data_dosage],axis=1)#[~x_y_data_nan]
        x_data_names=np.array(['const']+covar.iloc[:,2:].columns.values.tolist()+['THIS'])
        #x_data_null=np.concatenate([x_data_intercept,x_data_covariate],axis=1)#[~x_y_data_nan]
        y_data=pheno['pheno']
        
        family=(sm.families.Gaussian() if assoc=='linear' else sm.families.Binomial())
        model=sm.GLM(y_data,x_data, family=family,missing='drop')
        #model_null=sm.GLM(y_data,x_data_null, family=family,missing='drop')
        model_result=model.fit()
        
        
        for model_result_idx in range(len(model_result.params)):
            
            assoc_result_record(marker_name=marker,
                                P=model_result.pvalues.iloc[model_result_idx],
                                coef=model_result.params.iloc[model_result_idx],
                                std=model_result.bse.iloc[model_result_idx],
                                Z=model_result.tvalues.iloc[model_result_idx],
                                term=x_data_names[model_result_idx],
                                nobs=model_result.nobs,
                                note='plink bialleic'
                               
                               )
        #model_result_null=model_null.fit()            


# In[234]:


pd.DataFrame(assoc_result_list)


# In[215]:


model_result.nobs,model_result.params.iloc[-1],model_result.bse.iloc[-1],model_result.tvalues.iloc[-1],model_result.pvalues.iloc[-1]#.summary2()


# In[220]:


x_data.shape,len(x_data_names)


# In[151]:


for marker in plink_bialleic_list:
    plink
    #name,allele=(re_exp.search(marker).group('name'),re_exp.search(marker).group('allele')) if re_exp.search(marker) is not None else (None,None)
    #if name is not None and name.split('*')[0].split('_')[1]==HLA_gene:
    #    marker_dosage_list=marker_dosage_list_dict.get(name,[])
    #    marker_dosage_list.append({"marker":marker,"dosage":plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker)})
    #    marker_dosage_list_dict[name]=marker_dosage_list


# In[153]:


plink_bialleic_list


# In[ ]:


for idx_bim,(SNP,row) in enumerate(plink_bim.iterrows()):
    y_data=pheno['pheno'].replace(-9,np.nan)
    
    
    x_data_intercept=np.array([np.ones(plink_KCHIP_HLA_AA_SNP_1000G_fam.shape[0])]).transpose()    
    x_data_dosage=plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(SNP).astype(float);x_data_dosage[x_data_dosage==-1]=np.nan
    x_data_dosage=np.expand_dims(x_data_dosage,axis=1)
    x_data_covariate=covariate_plink_df.values
    #print(x_data_intercept.shape,covariate_plink_df.values.shape)
    x_data=np.concatenate([x_data_intercept,x_data_covariate,x_data_dosage],axis=1)#[~x_y_data_nan]
    x_data_null=np.concatenate([x_data_intercept,x_data_covariate],axis=1)#[~x_y_data_nan]
    
    family=(sm.families.Gaussian() if phenotype_type=='continuous' else sm.families.Binomial())
    model=sm.GLM(y_data,x_data, family=family,missing='drop')
    model_null=sm.GLM(y_data,x_data_null, family=family,missing='drop')
    model_result=model.fit()
    model_result_null=model_null.fit()
    
    print(model_result.summary())
    if idx_bim==2:
        break
    


# In[ ]:





# In[ ]:





# In[24]:


log.info_head("Checking missing values(observations)")


# In[ ]:


conditional_omnibus_list=[]

for HLA_gene in conditional_variant_list_HLA_AA:
    marker_dosage_list_dict={}
    for re_exp in [HLA_re_exp,AA_re_exp]:
        for marker in plink_KCHIP_HLA_AA_SNP_1000G_bim.index:
            name,allele=(re_exp.search(marker).group('name'),re_exp.search(marker).group('allele')) if re_exp.search(marker) is not None else (None,None)
            if name is not None and name.split('*')[0].split('_')[1]==HLA_gene:
                marker_dosage_list=marker_dosage_list_dict.get(name,[])
                marker_dosage_list.append({"marker":marker,"dosage":plink_KCHIP_HLA_AA_SNP_1000G.get_geno_marker(marker)})
                marker_dosage_list_dict[name]=marker_dosage_list

    for key,marker_dosage_list in marker_dosage_list_dict.items():
        marker_list=[marker_dosage['marker'] for marker_dosage in marker_dosage_list]
        dosage_array=np.array([marker_dosage['dosage'] for marker_dosage in marker_dosage_list]).transpose()
        if len(marker_list)>1:
            marker_list_cut=np.delete(marker_list,dosage_array.sum(axis=0).argmax())
            dosage_array_cut=np.delete(dosage_array,dosage_array.sum(axis=0).argmax(),axis=1)
        else:
            marker_list_cut=marker_list
            dosage_array_cut=dosage_array

        for i,marker in enumerate(marker_list_cut):
            covariate_plink_df[marker]=dosage_array_cut[:,i]

    conditional_omnibus_list+=marker_dosage_list_dict.keys()


# In[25]:


def marker_data_to_onehot(aa_marker_data):
    
    aa_marker_data_unique=np.unique(aa_marker_data)
    aa_marker_data_unique_nonan=aa_marker_data_unique[aa_marker_data_unique!='nan'].tolist()
    aa_marker_data_unique_nan=aa_marker_data_unique_nonan+['nan']
    
    aa_marker_data_int_nan=list(map(lambda x: aa_marker_data_unique_nan.index(x),aa_marker_data))
    
    aa_marker_data_onehot_nan=np.zeros((len(aa_marker_data),len(aa_marker_data_unique_nan)))
    aa_marker_data_onehot_nan[np.arange(len(aa_marker_data)),aa_marker_data_int_nan]=1
    aa_marker_data_onehot_nonan=aa_marker_data_onehot_nan[:,:-1]
    
    return aa_marker_data_unique_nonan,aa_marker_data_onehot_nonan

def prepare_onehot(aa_marker_data_unique_nonan,aa_marker_data_onehot_nonan,set_nan=True,cut_mostfrequent=True):
    assert len(aa_marker_data_unique_nonan)==aa_marker_data_onehot_nonan.shape[1]
    assert np.isnan(aa_marker_data_onehot_nonan).sum()==0
    assert (np.nan not in aa_marker_data_unique_nonan) and ('nan' not in aa_marker_data_unique_nonan)
    
    aa_marker_data_onehot_nonan_sumrow=np.sum(aa_marker_data_onehot_nonan,axis=1)
    aa_marker_data_onehot_nonan_sumcol=np.sum(aa_marker_data_onehot_nonan,axis=0)
    #print(aa_marker_data_onehot_nonan_sumcol,np.argmax(aa_marker_data_onehot_nonan_sumcol))
    #print(aa_marker_data_onehot_nonan_sumrow.shape,aa_marker_data_onehot_nonan_sumcol.shape)
    #print(np.argmax(aa_marker_data_onehot_nonan_sumcol))
    if set_nan:
        aa_marker_data_onehot_nonan[aa_marker_data_onehot_nonan_sumrow!=1,:]=np.nan
    if cut_mostfrequent:
        aa_marker_data_unique_nonan=np.delete(aa_marker_data_unique_nonan, np.argmax(aa_marker_data_onehot_nonan_sumcol))
        aa_marker_data_onehot_nonan=np.delete(aa_marker_data_onehot_nonan, np.argmax(aa_marker_data_onehot_nonan_sumcol), axis=1)
    return aa_marker_data_unique_nonan,aa_marker_data_onehot_nonan


# ## y data

# In[26]:


y_data_pheno=np.repeat(pheno['pheno'].values,2)


# In[27]:


y_data_pheno_nan=np.isnan(y_data_pheno)


# In[28]:


if np.var(y_data_pheno[~y_data_pheno_nan])==0:
    log.error("No variance in y_data")
    raise


# In[29]:


log.info("missing values in y_data_pheno: {}".format(y_data_pheno_nan.sum()))


# ## x data

# In[30]:


x_data_intercept=np.array([np.ones(2*fam.shape[0])]).transpose()


# In[31]:


x_data_covar=np.repeat(covar.iloc[:,2:].values,2,axis=0)


# In[32]:


x_data_covar_nan=np.any(np.isnan(x_data_covar),axis=1)
log.info("missing values in x_data_covar: {}".format(x_data_covar_nan.sum()))


# In[33]:


x_data_covar_varcheck=np.array([np.var(covar) for covar in x_data_covar[~x_data_covar_nan].transpose()])==0

if np.any(x_data_covar_varcheck):
    log.info("Removed covariates of no variance"+", ".join(covar.columns[2:][x_data_covar_varcheck].tolist()))
    x_data_covar=np.delete(x_data_covar,np.arange(len(x_data_covar_varcheck))[x_data_covar_varcheck],axis=1)


# In[34]:


if (x_data_covar.size!=0) and (np.linalg.matrix_rank(x_data_covar)<x_data_covar.shape[1]):
    log.info("duplicated covariates were found. (rank(covariates)< # of covarirates))")
    #raise
    #x_data_covar=np.delete(x_data_covar,np.arange(len(x_data_covar_varcheck))[x_data_covar_varcheck],axis=1)


# In[35]:


temp_list=[x_data_intercept[:,0:0]]
for aa_marker_name in condition_list:
    aa_marker_data=aa_marker_data_list[aa_marker_name_list.index(aa_marker_name)]

    aa_marker_data_unique,aa_marker_data_onehot=marker_data_to_onehot(aa_marker_data)
    aa_marker_data_unique_cut,aa_marker_data_onehot_cut=prepare_onehot(aa_marker_data_unique,aa_marker_data_onehot,set_nan=True,cut_mostfrequent=True)        
    temp_list.append(aa_marker_data_onehot_cut)
    
x_data_condition=np.concatenate(temp_list,axis=1)
x_data_condition_nan=np.any(np.isnan(x_data_condition),axis=1)  


# In[36]:


log.info("missing values in x_data_condition: {}".format(x_data_condition_nan.sum()))


# In[37]:


"""
x_data_condition=np.delete(x_data_condition,list(todel),axis=1)
todel=set()
for i in range(x_data_condition.shape[1]):
    for j in range(i+1,x_data_condition.shape[1]):
        if np.corrcoef(x_data_condition[~x_data_condition_nan][:,i],x_data_condition[~x_data_condition_nan][:,j])[1,0]>0.9:
            todel.add(i)
            todel.add(j)
x_data_condition_nan=np.any(np.isnan(x_data_condition),axis=1)  
np.linalg.matrix_rank(x_data_condition[~x_data_condition_nan]),x_data_condition.shape[1]
"""            


# In[38]:


log.info_head("Regression")


# In[ ]:


log.info('[{:3d}/{:3d}] {:10s} {:15s} {:5s} {:.5s}({}) {:.5s}'.format(
                            0,
                            len(aa_marker_name_list_aaonly),
                            'ID',
                            'residues',
                            'n_obs',
                            'chisq',
                            'df',
                            'P'
                        ))
assoc_result_list=[]

for idx,aa_marker_name in enumerate(aa_marker_name_list_aaonly):
    aa_marker_data=aa_marker_data_list[idx]
    aa_marker_data_unique,aa_marker_data_onehot=marker_data_to_onehot(aa_marker_data)
    aa_marker_data_unique_cut,aa_marker_data_onehot_cut=prepare_onehot(aa_marker_data_unique,aa_marker_data_onehot,set_nan=True,cut_mostfrequent=True)
    
    x_data_aa_marker=aa_marker_data_onehot_cut
    x_data_aa_marker_nan=np.any(np.isnan(x_data_aa_marker),axis=1)
    
    x_data_nan=np.logical_or.reduce([x_data_covar_nan,x_data_condition_nan,x_data_aa_marker_nan])
    
    y_data=y_data_pheno
    y_data_nan=y_data_pheno_nan    
    x_y_data_nan=(x_data_nan)|(y_data_nan)


    x_data_null=np.concatenate([x_data_intercept,x_data_covar,x_data_condition,x_data_aa_marker],axis=1)[~x_y_data_nan]
    x_data_alt=np.concatenate([x_data_intercept,x_data_covar,x_data_condition],axis=1)[~x_y_data_nan]
    y_data=y_data[~x_y_data_nan]
    
    
    family=(sm.families.Gaussian() if assoc=='linear' else sm.families.Binomial())
    model_null = sm.GLM(y_data,x_data_null, family=family,missing='raise')
    model_alt = sm.GLM(y_data,x_data_alt, family=family,missing='raise')
    
    try:
        result_alt = model_alt.fit()
        result_null = model_null.fit()
    except sm.tools.sm_exceptions.PerfectSeparationError as e:
        nobs=np.nan
        chisq_diff=np.nan
        df_diff=np.nan
        p_value=np.nan
    else:
        assert result_alt.nobs==result_null.nobs
        nobs=result_alt.nobs
        chisq_diff=2*(result_null.llf-result_alt.llf)
        df_diff=result_null.df_model-result_alt.df_model
        p_value=chi2.sf(chisq_diff,df_diff)        
        #p_value=1 - chi2.cdf(chisq_diff,df_diff)

    #print(result_alt.summary())
    assoc_result={'idx':idx+1,
                  'ID':aa_marker_name,
                  'residues':','.join(aa_marker_data_unique),
                  'n_obs':nobs,
                  'chisq':chisq_diff,
                  'df':df_diff,
                'P':p_value}
    
    
    assoc_result_list.append(assoc_result)
    
    log.info('[{:3d}/{:3d}] {:10s} {:15s} {:5f} {:.5f}({}) {:e}'.format(
                                assoc_result['idx'],
                                len(aa_marker_name_list_aaonly),
                                assoc_result['ID'],
                                assoc_result['residues'],
                                assoc_result['n_obs'],
                                assoc_result['chisq'],
                                assoc_result['df'],
                                assoc_result['P']
                            ))


# In[ ]:


pd.DataFrame(assoc_result_list)[['idx','ID','residues','n_obs','chisq','df','P']].to_csv(out+'.'+assoc,sep='\t',index=None)


# In[9]:


log.info("End time: "+time.strftime('%c', time.localtime(time.time())))

