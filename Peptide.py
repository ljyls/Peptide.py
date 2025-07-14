#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import os
import re

print(os.getcwd())
 
def auto_create_path(FilePath):
    if os.path.exists(FilePath): 
           print( 'dir exists'  )
    else:
           print( 'dir not exists')
           os.makedirs(FilePath) 

auto_create_path('./Peptide')
auto_create_path('./Peptide/Peptides')
auto_create_path('./Peptide/Peptides_sorting')
auto_create_path('./Peptide/Peptides_fasta')
auto_create_path('./Peptide/Mature_Peptides')

os.chdir(r'./Peptide/Peptides')

#-------------------------------------------------------------------------------Import peptide information-------------------------------------------------------------------------------


data=pd.read_excel('../../Search_file_export.xlsx', sheet_name='Proteins', header=0)
data=data.rename(columns={'Accession':'Master_accession'})
data['Master']=data['Master'].replace('Master Protein', 'IsMasterProtein')
data['Checked']=data['Checked'].replace(0.0, 'False')
data.columns=data.columns.str.replace(r'(^Abundance.*Sample$)', 'Abundance',regex=True)
data.reset_index(drop=True,inplace=True)
str_FALSE=data[(data.Checked=='False')].index.tolist()
for i in str_FALSE:
    data.loc[i+1,'Checked']=data.loc[i,'Master_accession']
data['Checked']=data['Checked'].replace('', np.nan).ffill(axis=0)

#Remove proteins with medium or low FDR confidence, if any
data.reset_index(drop=True,inplace=True)
l_end=data.shape[0]
if (data['Protein FDR Confidence: Combined']=='Medium').any():
    data_cut=data[(data['Protein FDR Confidence: Combined']=='Medium')].index.tolist()
else:
    data_cut=[l_end]
data=data.iloc[0:data_cut[0]]

data.reset_index(drop=True,inplace=True)
l_end=data.shape[0]
if (data['Protein FDR Confidence: Combined']=='Low').any():
    data_cut=data[(data['Protein FDR Confidence: Combined']=='Low')].index.tolist()
else:
    data_cut=[l_end]
data=data.iloc[0:data_cut[0]]

#+++++++++++++Split peptides by protein into files+++++++++++++
#Generate a sequence collection file (protein + peptide list) for each protein to the Peptides folder
data.reset_index(drop=True,inplace=True)
rows=data.shape[0]
target_col = 'Checked'
cols_list = []
 
for i in range(0, rows):
    temp = data[target_col][i]
    if temp not in cols_list:
        cols_list.append(temp)
 
for col in cols_list:
    new_df=pd.DataFrame()
 
    for i in range(0, rows):
        if data[target_col][i] == col:
            new_df = pd.concat([new_df, data.iloc[[i], :]], axis=0, ignore_index=True)
 
    new_df.to_csv( f'{str(col)}.txt', header=True, index=False, sep='\t', na_rep='')

os.replace("False.txt","../Master_proteins.txt")

dfP_all=pd.read_csv('../Master_proteins.txt',header=0,usecols=['Master','Master_accession','Sequence'],sep='\t')

dfP=dfP_all[['Master_accession','Sequence']].copy()
dfP=dfP.rename(columns={'Master_accession':'Labels'})

dfP_all=dfP_all.rename(columns={'Sequence':'Protein'})
dfP_all['Protein_length']=dfP_all['Protein'].apply(len)

#First filtering
def odd_count(content):
    pattern=re.compile(r'(\d+)xDehydro\s\[C')
    result=pattern.findall(content)
    result=list(map(eval, result))
    odd = []
    for i in result:
        if i % 2 != 0:
            odd.append(i)
    return len(odd)

filenames=os.listdir()
for i in filenames:
    pi=i.rstrip('.txt')
    df=pd.read_csv(f'{pi}.txt',header=1,sep='\t')
    df.rename(columns={pi:"Master_accession"}, inplace=True)
    df.columns=df.columns.str.replace(r'(^Abundance.*Sample$)', 'Abundance',regex=True)
    df['Modifications'].fillna('N',inplace=True)
    df['odd_count']=df['Modifications'].apply(lambda x : odd_count(x))
    df=df.loc[~(df['odd_count'] > 0)]#Delete peptides with ‘Odd xDehydro’ in each protein file
    df=df.drop('odd_count',axis=1)
    df=df.loc[(df['PSM Ambiguity'] == 'Unambiguous')]# “Unambiguous” for “PSM Ambiguity” and “High” for “Confidence (by Search Engine): Sequest HT” are retained
    df=df.loc[(df['Confidence (by Search Engine): Sequest HT'] == 'High')]
    df.sort_values(['Abundance','# PSMs'],ascending=[False,False],inplace=True, na_position='last')
    df.to_csv( f'{pi}.txt', header=True, index=False, sep='\t', na_rep='')

#Delete protein files without peptide
filenames=os.listdir()
for i in filenames:
    df_drop=pd.read_csv(i,header=0,sep='\t')
    row_count=df_drop.shape[0]
    if row_count==0:
        os.remove(i)

#combine the peptide isoforms into peptides without considering modification of peptide isoforms and accumulated their abundances
#generate the most abundant peptide among all exclusive peptides of each protein, which is deemed the candidate mature peptide for that protein. The exclusive peptide is required to be not classified as “No Quan Values”, “Not Reliable”, or “Shared” in the “Quan Info” column
filenames=os.listdir()
for i in filenames:
    pi=i.rstrip('.txt')
    df=pd.read_csv(f'{pi}.txt',header=0,sep='\t')   
    df=df[['Master_accession','Positions in Proteins','Sequence','Abundance','Modifications','# PSMs','Quan Info']]
    df.reset_index(drop=True,inplace=True)
    df.index=df.index+1
    df['index_Quan Info_split']=df.index
    df['Quan Info']=df['Quan Info'].replace(['No Quan Values','Not Reliable','Shared'],['N','N','S'])
    df['Quan Info'].fillna('E', inplace=True)
    df['Quan Info']=df['index_Quan Info_split'].map(str)+"-"+df['Quan Info'].map(str)
    df['Modifications']=df['index_Quan Info_split'].map(str)+"-"+df['Modifications'].map(str)
    df=df.drop('index_Quan Info_split',axis=1)
    df['Quan Info']=df.groupby(by=['Sequence'])['Quan Info'].transform(lambda x : '^'.join(x))
    df['Modifications']=df.groupby(by=['Sequence'])['Modifications'].transform(lambda x : '^'.join(x)) 
    dfQ=df[['Sequence','Quan Info','Modifications']].copy()
    dfQ.drop_duplicates(keep='first',inplace=True)
    dfPos=df[['Master_accession','Positions in Proteins','Sequence']].groupby(by=['Sequence']).agg(lambda x: '; '.join(set(x))).reset_index()
    dfQ=pd.merge(dfPos,dfQ,how='left',on='Sequence')
    dfQ.drop_duplicates(keep='first',inplace=True)

    dfsort=df[['Sequence','Abundance','# PSMs']].copy()
    dfsort['Abundance']=dfsort['Abundance'].astype(float)
    dfsort['# PSMs']=dfsort['# PSMs'].astype(int)
    dfsort=dfsort.groupby(by=['Sequence']).agg({'# PSMs':'sum','Abundance':'sum'}).reset_index()
    dfsort['Sequence_length']=dfsort['Sequence'].apply(len)
    dfsort.sort_values(['# PSMs','Abundance','Sequence_length'],ascending=[False,False,False],inplace=True, na_position='last')
    dfsort.reset_index(drop=True,inplace=True)
    dfsort.index=dfsort.index+1
    dfsort['# PSMs_sort']=dfsort.index
    dfsort['# PSMs_sort']=dfsort['# PSMs_sort'].map(str)+"-"+dfsort['# PSMs'].map(str)
    
    dfsort.sort_values(['Sequence_length','Abundance','# PSMs'],ascending=[False,False,False],inplace=True, na_position='last')
    dfsort.reset_index(drop=True,inplace=True)
    dfsort.index=dfsort.index+1
    dfsort['Sequence_length_sort']=dfsort.index
    dfsort['Sequence_length_sort']=dfsort['Sequence_length_sort'].map(str)+"-"+dfsort['Sequence_length'].map(str)  
    
    dfsort.sort_values(['Abundance','# PSMs','Sequence_length'],ascending=[False,False,False],inplace=True, na_position='last')
    dfsort['AbundanceE']=dfsort['Abundance'].apply(lambda x: '{:.1e}'.format(x))
    dfsort.reset_index(drop=True,inplace=True)
    dfsort.index=dfsort.index+1
    dfsort['Abundance_sort']=dfsort.index
    dfsort['Abundance_sort']=dfsort['Abundance_sort'].map(str)+"-"+dfsort['AbundanceE'].map(str)
    dfsort=dfsort.drop('AbundanceE',axis=1)
    
    dfsort=pd.merge(dfsort,dfQ,how='left',on='Sequence')

    dfsort['Labels']=dfsort[['Abundance_sort','# PSMs_sort','Sequence_length_sort','Quan Info']].agg('|'.join, axis=1)
    dfsort.to_csv(f'../Peptides_sorting/Sort_{pi}.txt', index=False, header=True, sep='\t', na_rep='')
    
    dfsort=dfsort[['Labels','Sequence']]
    dfPi=dfP[dfP['Labels']==pi]
    df_fasta=pd.concat([dfPi, dfsort], axis=0)
    df_fasta.to_csv(f'../Peptides_fasta/{pi}.txt', index=False, header=False, sep='\t', na_rep='')

#Second filtering
os.chdir(r'../Peptides_sorting')

filenames=os.listdir()

dfPlist=None
for i in filenames:
    pi=i.rstrip('.txt')
    dfP1=pd.read_csv(f'{pi}.txt',header=0,sep='\t')
    dfP1=dfP1.loc[dfP1['Quan Info'].str.contains('E')]
    dfP1=dfP1.iloc[:1]
    if dfPlist is None:
       dfPlist=dfP1
    else:
      dfPlist=pd.concat([dfPlist,dfP1], axis=0)

dfPlist.rename(columns={'Sequence':'Peptide','Abundance':'Peptide_abundance','Positions in Proteins':'Peptide_positions_in_Proteins','Labels':'Peptide_labels_in_Master_Protein','Modifications':'Peptide_Modifications'}, inplace=True)      
#dfPlist.to_csv(f'../Max_abundance_exclusive_peptide.txt', index=False, header=True, sep='\t', na_rep='')
dfPlist['Peptide_length']=dfPlist['Peptide'].apply(len)
dfPlist=dfPlist[['Peptide','Peptide_length','Peptide_abundance','Master_accession','Peptide_positions_in_Proteins','Peptide_labels_in_Master_Protein','Peptide_Modifications']]

dfP_all=pd.merge(dfP_all,dfPlist,how='left',on='Master_accession')
dfP_all.dropna(subset=['Peptide'],inplace=True)
dfP_all=dfP_all[['Master_accession','Protein','Protein_length','Peptide','Peptide_length','Peptide_labels_in_Master_Protein','Peptide_Modifications','Peptide_positions_in_Proteins','Peptide_abundance']]
dfP_all.sort_values(['Protein_length'],ascending=[False],inplace=True, na_position='last')
dfP_all.drop_duplicates(subset=['Peptide'],keep='first',inplace=True)
dfP_all.sort_values(['Master_accession'],ascending=[True],inplace=True, na_position='last')
dfP_all.fillna(value=0,inplace=True)
dfP_all[['Protein_length','Peptide_length','Peptide_abundance']]=dfP_all[['Protein_length','Peptide_length','Peptide_abundance']].astype(int)

dfT=dfP_all[['Master_accession','Peptide_positions_in_Proteins']].copy()
dfT=dfT.join(dfT['Peptide_positions_in_Proteins'].str.split('; ', expand=True).stack().reset_index(level=1, drop=True).rename('Peptide_positions_in_Proteins_split'))
dfT=dfT.drop('Peptide_positions_in_Proteins',axis=1)
dfT.rename(columns={'Peptide_positions_in_Proteins_split':'Peptide_positions_in_Proteins'}, inplace=True)
dfT['Peptide_positions_in_Proteins']=dfT['Peptide_positions_in_Proteins'].str.strip()
dfT.dropna(subset=['Peptide_positions_in_Proteins'],inplace=True)
dfT.drop_duplicates(keep='first',inplace=True)
dfT.sort_values(['Master_accession','Peptide_positions_in_Proteins'],ascending=[True,True],inplace=True, na_position='last')
dfT.reset_index(drop=True,inplace=True)
dfT=dfT[['Peptide_positions_in_Proteins','Master_accession']].groupby(by=['Master_accession']).agg(lambda x: ';'.join(x)).reset_index()

dfP_all=dfP_all.drop('Peptide_positions_in_Proteins',axis=1)
dfP_all=pd.merge(dfP_all,dfT,how='left',on='Master_accession')
dfP_all['Peptide_abundance']=dfP_all['Peptide_abundance'].astype('Int64')
#dfP_all.to_csv('../Peptide_filtered.txt', index=False, header=True, sep='\t', na_rep='')

print('done ++++++++++++++++++Peptide_filtered++++++++++++++++++')

#-----------------------------Remove mature peptide sequences with inclusion relationships, if any--------------------------------------

os.chdir(r'../Mature_Peptides')
dfP_all.sort_values(['Peptide_length','Peptide_abundance','Master_accession'],ascending=[False,False,True],inplace=True, na_position='last')

dfV=dfP_all[['Master_accession','Peptide']].copy()
dflist=pd.DataFrame(dfV['Master_accession'])
dfV.drop_duplicates(keep='first',inplace=True)
dfV.reset_index(drop=True,inplace=True)
dfV_clone=dfV.copy()
dfV_clone=dfV_clone.rename(columns={'Master_accession':'Original_accession'})
dfV_treat=dfV.copy()

for i, row in dfV.iterrows():
    tmp=dfV_treat[dfV_treat['Peptide'].str.contains(row['Peptide'])]
    dfV_clone.loc[i, 'Master_accession']=tmp['Master_accession'].values[0]
    
dfV_clone=pd.DataFrame(dfV_clone['Master_accession'])
dfV_clone.drop_duplicates(keep='first',inplace=True)

dflist=dflist.append(dfV_clone)
dflist=dflist.append(dfV_clone)
dflist.drop_duplicates(keep=False,inplace=True)
#dflist.to_csv('./discarded.list', index=False, header=None, sep='\t', na_rep='')

df_discard=pd.merge(dfP_all,dflist,how='right',on='Master_accession')
#df_discard.to_csv('../Peptide_removed_by_refining.txt', index=False, header=True, sep='\t', na_rep='')

dfP_all=pd.merge(dfP_all,dfV_clone,how='right',on='Master_accession')
dfP_all.to_csv('./Candidate_Mature_Peptides.table', index=False, header=True, sep='\t', na_rep='')
#dfP_all[['Master_accession','Peptide']].to_csv('./Mature_Peptides_sequenceONLY.table', index=False, header=None, sep='\t')

#--------------------------------Delete the files where the filtered proteins are located-----------------------------------------------

os.chdir(r'../')

data=[]
for file in sorted(os.listdir('Peptides')):
    data.append(file)
list1=pd.DataFrame(data, columns=['Master_accession'])
list1.drop_duplicates(keep='first',inplace=True)

list2=pd.DataFrame(dfP_all['Master_accession']+'.txt')
list2.drop_duplicates(keep='first',inplace=True)

list1=list1.append(list2)
list1=list1.append(list2)
list1.drop_duplicates(keep=False,inplace=True)
dflist=list1['Master_accession'].values.tolist()

from pathlib import Path

p=Path(r'./')
for i in dflist:
    for file in p.rglob('*'+i+'*'):    
        if os.path.isfile(file):
            os.remove(file)


print('done ++++++++++++++++++Peptide_refined++++++++++++++++++')

