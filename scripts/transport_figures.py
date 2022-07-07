#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 10:34:36 2022

@author: alexandra
"""
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import os
import shutil
from scipy.signal import savgol_filter
from os import path

from scipy import special

def find_id (row):
    return row['job'].split('-id')[1]

##########################################################
### Parameters of the script : can be modified
##########################################################

# define the folder with the simulations results for the transport analysis
rep='../data/simulations/'
analysis='-d7e-08-l6e-02'
rep=rep+'intakeRandomWT10t200area'+analysis+'/'

#associated diffusion coefficient
#D=1.680000e-07
D=6.800000e-08 #cm2/s

#chose the state to plot
stages=['baseline','stageNREM']#'stageIS',,'stageREM'

#chose the frequency bands to plot
bandnames=[ 'Card','LFVLFCard','LFVLF']

#name of the figure
namefigure='transport_'+'_d%.1e'%D

# output folder
outputfolder='../output/figures/'


#style

plt.rcParams.update({'font.size': 18})

my_pal = {"baseline":"gray","stageREM": "darkorchid", "stageNREM": "mediumseagreen", "stageIS":"darkorange","stageAwakening":"blue"}

##########################################################
### script
##########################################################

# Directory to store output data             
if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)

# get all the log files
from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir(rep) if isfile(join(rep, f))]
islog=[string.endswith('.log') for string in onlyfiles]
whereislog=np.where(islog)[0]
logfiles=[onlyfiles[index] for index in whereislog]

import re
end=[re.findall(r'id\S+PVSinfo.log',string) for string in logfiles]

vesselIDS=[re.findall(r'id\S+PVSinfo.log',string)[0].replace('_PVSinfo.log','') for string in logfiles if string[0:4]=='disp']

vesselIDS=np.unique(vesselIDS)

errorfiles=[]

factortime=1#1/60  

plt.figure()

# check computed diffusion

file=[string for string in logfiles if string[0:4]=='diff']

if len(file) :
    job=file[0].replace('_PVSinfo.log','')
    file=job+'_concentration.txt'
    Data=np.loadtxt(rep+file,delimiter=',')

    time=Data[1:,0]
    concentration=Data[1:,1:]
    x=Data[0,1:]
    
    xth=[0]
                
    for j,t in enumerate(time[1:]) :
                    
        th=0.1 #concentration[j+1,:].max()/10
                    
        ith=np.where(concentration[j+1,:]<th)[0][0]
                    
        xth.append(x[ith])
                    
                
                
    # interpolate on a time scale
    spantime=np.linspace(0,200,100)
            
    smoothedxth=savgol_filter(xth,81,3) 
    spanxth=np.interp(spantime,time,smoothedxth)
    #plt.plot(spantime*factortime, spanxth*1e4,c='c',alpha=1,linewidth=3,label='simulated diffusion')
            


for bandname in bandnames :

    plt.figure()
    
    for stage in stages :
        
        list_xth=[]
        list_time=[]
        for vessel in vesselIDS :
            
            job='disp'+analysis+'-'+stage+'-'+bandname+'-'+vessel
            print(job)
            print(rep)
            
            file=job+'_concentration.txt'

            try :
                Data=np.loadtxt(rep+file,delimiter=',')

                time=Data[1:,0]
                concentration=Data[1:,1:]
                x=Data[0,1:]
                
                
            except :
                print('Error during the loading of '+file)
                errorfiles.append([file])
                continue
            
            xth=[0]
            
            for j,t in enumerate(time[1:]) :
                
                th=0.1 #concentration[j+1,:].max()/10
                
                ith=np.where(concentration[j+1,:]<=th)[0][0]
                
                xth.append(x[ith])
                
                
            #plt.plot(time,xth,'k',alpha=0.3)
            list_xth.append(xth)
            list_time.append(time)
            
            
        # interpolate on a time scale
        spantime=np.linspace(0,200,50)
        #convert to array
        xtharray=np.zeros((len(list_xth),len(spantime)))
        
        for i,xth in enumerate(list_xth) :
            smoothedxth=savgol_filter(list_xth[i],5,3) 
            xtharray[i,:]=np.interp(spantime,list_time[i],smoothedxth)
            plt.plot(spantime*factortime,xtharray[i,:]*1e4,c=my_pal[stage],alpha=0.1)
            smoothedxth=savgol_filter(list_xth[i],101,3) 
            xtharray[i,:]=np.interp(spantime,list_time[i],smoothedxth)
            
            
        xthmean=np.median(xtharray,axis=0)
        xthq1=np.percentile(xtharray,10,axis=0)
        xthq3=np.percentile(xtharray,90,axis=0)
        
        
        
        plt.fill_between(spantime*factortime,xthq1*1e4,xthq3*1e4,alpha=0.3,color=my_pal[stage])
        plt.plot(spantime*factortime,xthmean*1e4,label='dispersion '+bandname+'_'+stage,color=my_pal[stage],linewidth=3)
        
    
    spantimeth=np.linspace(0,200,1000)    
    diffusion=special.erfcinv(0.1)*2*np.sqrt(D*spantimeth)
    # diffusion06=special.erfcinv(0.1)*2*np.sqrt(D*1.3*spantimeth)
    plt.plot(spantimeth*factortime,diffusion*1e4,'k',label='diffusion')
    #plt.plot(spantimeth*factortime,diffusion06*1e4,'k',label='diffusion enhanced 30 pc')
    
    
    
    # #plt.plot(spantime*factortime,spantime*0.5,'k-.',label='net flow 0.1 um /s')
    plt.plot(spantime*factortime,spantime*1,'k-.',label='net flow 1 um /s')
    plt.plot(spantime*factortime,spantime*5,'k:',label='net flow 5 um /s')
    plt.plot(spantime*factortime,spantime*10,'k--',label='net flow 10 um /s')
    
    plt.ylim([0,130])
    plt.xlim([0,120])
    
    plt.xlabel('time (s)')
    plt.ylabel('depth (um)')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
          fancybox=True, shadow=True, ncol=3)
        
            
    plt.savefig(outputfolder+namefigure+'-'+bandname+'.png' , bbox_inches='tight') 
    plt.savefig(outputfolder+namefigure+'-'+bandname+'.pdf' , bbox_inches='tight')
    
    #plt.title(stage + ' Free diffusion 7e-8cm2/s')
    plt.show()

            
            

