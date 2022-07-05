#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 14:12:50 2021

@author: alexandra
"""

import pickle
import pandas as pd
from matplotlib import pyplot as plt 
from scipy.stats import mannwhitneyu, normaltest
import seaborn as sbn
from statannot import add_stat_annotation
import numpy as np

import os
from os import listdir
from os.path import isfile, join

from scipy import stats

import itertools

sbn.set(style="whitegrid")  

#amplitudedb[amplitudedb['amp pc']>50][['mouse number','trial','bandname','tmin','amp pc']]


All=False
box=False
hist=False
violin=True

colors=['darkorchid','mediumseagreen','darkorange','blue']
my_pal = {"REM": "darkorchid", 
          "NREM": "seagreen", 
          "IS":"darkorange", 
          "Baseline":"royalblue",
          "REM WT": "darkorchid", 
          "NREM WT": "seagreen", 
          "IS WT":"darkorange", 
          "Baseline WT":"royalblue",
          "REM AQP4KO": "mediumorchid", 
          "NREM AQP4KO": "mediumseagreen", 
          "IS AQP4KO":"moccasin", 
          "Baseline AQP4KO":"cornflowerblue",
          "Quiet":"dodgerblue",
          "Locomotion":"tomato",
          "Whisking":"lightsalmon"}

legends={'period':'Oscillation period (s)','amp':'Oscillation amplitude (um)','mean':'Radius mean value','amp pc':'Oscillation amplitude (pc)','amp norm':'Normalized amplitude','period norm':'Normalized period','mean norm':'Normalized mean radius' }
        



stage_list={'stage':['Locomotion','Whisking','Quiet','Baseline', 'NREM','IS','REM'],'sleep':['Baseline', 'NREM','IS','REM'], 'stage mouse':['Baseline WT','Baseline AQP4KO', 'NREM WT','NREM AQP4KO','IS WT','IS AQP4KO','REM WT','REM AQP4KO']}


allstages=[('Locomotion','Quiet'),
                       ('Whisking','Quiet'),
                       ('Quiet','Baseline'),
                     ('Baseline', 'NREM'),
                     ('Baseline', 'IS'),
                     ('Baseline', 'REM')]

awakestages=[('Locomotion','Quiet'),
             ('Whisking','Quiet'),
             ('Locomotion','Whisking')]


sleepstages=[('Baseline', 'NREM'),
              ('Baseline', 'IS'),
              ('Baseline', 'REM')]




# set the data directory 
dir='/home/alexandra/Documents/Python/linescan-analysis/output/databases/PenetratingArteriolesWT6/'
#dir='/home/alexandra/Documents/Python/linescan-analysis/output/databases/oldsleepWToldmeanIF/'


outputdir='../output/statistics/PenetratingArteriolesWT6/'

# chose the analysis 
comparison_variable='stage' #'stage mouse'

# chose the variables
bandnamelist=['cardiac']
variablelist=['amp pc','mean','mean norm']


# Directory to store output data
if not os.path.exists(outputdir):
    os.makedirs(outputdir)


for scanline in ['PVS']: #,'lumen','endfoot'   

    # Directory to store output data
    if not os.path.exists(outputdir+scanline):
        os.makedirs(outputdir+scanline)

    ###get data
    
    # # just one scanline
    
    # mousekind='AQP4KO'
    # mousenumber='04'
    # trial='16'
    # nofinal=''
    # study=mousekind+'-'+mousenumber+'-'+trial+'-'+nofinal

    # amplitudedb=pd.read_pickle('../output/databases/amplitude_'+mousekind+'-'+mousenumber+'-'+trial+'-'+nofinal+'-'+scanline+'.pkl')
        
    
    # #several linescans

    
    # get the list of files
    
    files = [f for f in listdir(dir) if isfile(join(dir, f))]
    
    studylist=[]
    
    for f in files : 
        data=pd.read_pickle(dir+f)
        studylist.append(data)
        
    #concanate
    study='stage'

    amplitudedb=pd.concat([data for data in studylist], ignore_index=True)
    
    
    
        # add new variables
    amplitudedb['amp pc']=amplitudedb['amp']/amplitudedb['mean']*100
    
    amplitudedb['stage mouse']=amplitudedb['stage']+' '+amplitudedb['mouse kind']
    

    # filter

    filter=(amplitudedb['linescan']==scanline)
    
    filter=filter&(amplitudedb['amp pc']>0)
    filter=filter&(amplitudedb['hasnan']==False)
    
    filter=filter&(amplitudedb['corrlumen']>0.8)
    filter=filter&(amplitudedb['correndfoot']>0.8)
    
    #filter rigid motion 
    #rigidmotion= (amplitudedb['bandname']=='VLF')*1.5 + (amplitudedb['bandname']=='LF')*0.8 + (amplitudedb['bandname']=='resp')*0.4 + (amplitudedb['bandname']=='cardiac')*0.4
    #theta=np.arcsin(rigidmotion/amplitudedb['meanlumen'])
    #rigidartefact=amplitudedb['meanlumen']*(1-np.cos(theta))
    #filter=filter&(amplitudedb['amp']>rigidartefact)
    
    #filter=filter&(amplitudedb[filter]['mouse number']!='06')&(amplitudedb[filter]['mouse number']!='09')
    
    # remove outliers
    #filter=filter&(np.abs(stats.zscore(np.log(amplitudedb[['period norm','mean norm']]))) < 2).any(axis=1)
    
    #filter=filter&((amplitudedb['tend']-amplitudedb['tini'])>10)
    
    # remove the cardiac pulsation when too much noise
    #filter=filter& np.logical_not(((amplitudedb['bandname']=='cardiac')|(amplitudedb['bandname']=='resp'))&(amplitudedb['HFerror']>0.1))
    #filter=filter& np.logical_not(((amplitudedb['bandname']=='cardiac')|(amplitudedb['bandname']=='resp'))&((amplitudedb['HFerrormax'])>0.8))
    #filter=filter& np.logical_not((amplitudedb['bandname']=='cardiac')&(amplitudedb['tortuosity']>1.2))
    #filter=filter& np.logical_not((amplitudedb['bandname']=='cardiac')&(amplitudedb['rawvar']>0.15))
    #filter=filter& np.logical_not((amplitudedb['bandname']=='cardiac')&(amplitudedb['slopevar']>1.0))
    ##filter=filter& np.logical_not((amplitudedb['bandname']=='cardiac')&(amplitudedb['ampvar']>1))
    ##filter=filter& np.logical_not((amplitudedb['bandname']=='cardiac')&(amplitudedb['periodvar']>0.5))

    
    # remove outliers
    filter=filter& np.logical_not((amplitudedb['bandname']=='cardiac')&(amplitudedb['outlier']))
    # stronger requierement of correlation for cardiac
    filter=filter& np.logical_not((amplitudedb['bandname']=='cardiac')&(amplitudedb['corrlumen']>0.95))
    filter=filter& np.logical_not((amplitudedb['bandname']=='cardiac')&(amplitudedb['correndfoot']>0.95))
    
    
    #filtermice=((amplitudedb['mouse kind']=='WT')&(amplitudedb['mouse number']=='10'))|((amplitudedb['mouse kind']=='AQP4KO')&(amplitudedb['mouse number']=='04'))
    #filter=filter &filtermice
    
    #filter=filter &(np.abs(amplitudedb['amp pc'])<=100)
    
    #filter=filter&(amplitudedb['mouse kind']=='AQP4KO')
    
    
    if sum(filter) ==0:
        break
    
    amplitudedb=amplitudedb[filter]


            
    for bandname in bandnamelist:
        
        #stages=['Baseline', 'NREM', 'IS', 'REM', 'Quiet','Locomotion', 'Whisking']
        
        stages=np.unique([item for t in combinations[comparison_variable] for item in t])
        combis=[combi for combi in combinations[comparison_variable] if (combi[0] in stages)&(combi[1] in stages)]
        
        
    
        for variable in variablelist:
            
            if hist :
        
                plt.figure() 
                for stage,c in zip(stage_list[comparison_variable],colors):
                    filter=(amplitudedb[comparison_variable]==stage)&(amplitudedb['bandname']==bandname)
                    
                    dataplot=amplitudedb[filter]
                    
                    
                    width=(max(amplitudedb[variable][amplitudedb['bandname']==bandname])-min(amplitudedb[variable][amplitudedb['bandname']==bandname]))/15
                    sbn.histplot(data=dataplot,x=variable,palette="crest",color=c,label=stage,stat="percent",binwidth=width)
                
                
                
                plt.xlabel(legends[variable])
                plt.legend()
                plt.title(bandname+' '+study+' '+scanline)
                
                plt.savefig(outputdir+linescan+'/'+study+'-'+bandname+'_'+variable+'_hist.png')
                #plt.savefig(outputdir+linescan+'/'+study+'-'+bandname+'_'+variable+'_hist.svg')
                
                
            
            
            filter=(amplitudedb['bandname']==bandname)
            filter=filter&(amplitudedb['stage'].isin(stages))
            dataplot=amplitudedb[filter]
            #dataplot.replace([np.inf, -np.inf], np.nan, inplace=True)
            dataplot.dropna(axis=0,how='any')
            
            
            if box :
            
                plt.figure()
                
                yposlist = (dataplot.groupby([comparison_variable])[variable].median()).tolist()
                Nlist = (dataplot.groupby([comparison_variable])[variable].size()).tolist()                
                xposlist = np.array(range(len(yposlist)))-0.01
                stringlist = ['%.2f, N=%i'%(m,s) for m,s in zip(yposlist,Nlist)]
                
                
                ax = sbn.boxplot(x=comparison_variable, y=variable, data=dataplot,palette=my_pal,order=list(np.unique(dataplot[comparison_variable])))  
                plt.ylabel(legends[variable])
                plt.title(bandname+' '+study+' '+scanline)
                
                for i in range(len(stringlist)):
                    ax.text(xposlist[i], yposlist[i], stringlist[i])
    
                #combinations=list(itertools.combinations(list(np.unique(dataplot[comparison_variable])), 2))
                
                
                add_stat_annotation(ax, data=dataplot, x=comparison_variable, y=variable, order=list(np.unique(dataplot[comparison_variable])),box_pairs=combinations[comparison_variable],
                                    test='Mann-Whitney', comparisons_correction=None,  text_format='star',verbose=2)
                
                plt.tight_layout()
                
                plt.savefig(outputdir+scanline+'/'+study+'-'+bandname+'_'+variable+'_box.png')
                #plt.savefig(outputdir+linescan+'/'+study+'-'+bandname+'_'+variable+'_box.svg')
                
            if violin :
            
                plt.figure(figsize=(10, 5)) 
                
                
                yposlist = (dataplot.groupby([comparison_variable])[variable].median()[stages]).tolist()
                Nlist = (dataplot.groupby([comparison_variable])[variable].size()[stages]).tolist()                
                maxlist = (dataplot.groupby([comparison_variable])[variable].max()[stages]).tolist()                

                xposlist = np.array(range(len(yposlist)))
                stringlist = ['%.2f, N=%i'%(m,s) for m,s in zip(yposlist,Nlist)]
                
                
                ax=sbn.violinplot(x=comparison_variable, y=variable, data=dataplot, order=stages,alpha=0.3, palette=my_pal, cut=0)  
               # sbn.swarmplot(data=dataplot,x="stage", y=variable, color="white", edgecolor="gray")
                ax.set_xticklabels(ax.get_xticklabels(),rotation = 30)
                plt.ylabel(legends[variable])
                plt.title(bandname+' '+study+' '+scanline)
                
                for i in range(len(stringlist)):
                    ax.text(xposlist[i]-0.3, maxlist[i]*1.02, stringlist[i])
    
               # combinations=list(itertools.combinations(list(np.unique(dataplot[comparison_variable])), 2))
                
                add_stat_annotation(ax, data=dataplot, x=comparison_variable, y=variable, order=stages,box_pairs=combis,
                                    test='Mann-Whitney',  comparisons_correction=None, text_format='star',verbose=2)
                
                plt.tight_layout()
                
                #plt.show()
                
                plt.savefig(outputdir+scanline+'/'+study+'-'+bandname+'_'+variable+'_violin.png')
                #plt.savefig(outputdir+linescan+'/'+study+'-'+bandname+'_'+variable+'_violin.svg')
                
plt.show()