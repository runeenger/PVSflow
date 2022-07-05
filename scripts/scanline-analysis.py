#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 14:21:29 2021

@author: Alexandra Vallet



This scritps is used to post-process the traces of the linescan data. 
The raw data consists in the position of the two borders of the vessel and the two borders of the atrocytes with time.

Here we compute the radius and area of vessels, astrocyte endfeet tube and PVS. 
We reduce the noise and separate the signal into several frequency bands.
We then extract  amplitude and period for every individual pulsation.
The results are saved in a database.
"""

# Import libraries

import numpy as np
from matplotlib import pyplot as plt 
plt.ioff()

import pandas as pd
import scipy.signal as sg
import os

from scipy.signal import savgol_filter

import src.linescan.datanalysis as da

from matplotlib.patches import Rectangle

from os import listdir
from os.path import isfile, join

from scipy.stats import pearsonr

import regex as re

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


##########################################################
### Parameters of the script : can be modified
##########################################################


# set the data directory 
dir='../data/raw/160322_4traces/Veins/'
#dir='/../data/raw/160322_4traces/Pials/'
#dir='../data/raw/160322_4traces/PenetratingArterioles/'

# name the output database
name='Veins'

# format for the figures
outputformat='.pdf'

# Definition of the frequency bands. 
# Cutoff in Hertz

cardiac={'bandname':'cardiac','cutoff1':4,'cutoff2':15}

resp={'bandname':'resp','cutoff1':1,'cutoff2':4}

lowfreq={'bandname':'LF','cutoff1':0.3,'cutoff2':1}

verylowfreq={'bandname':'VLF','cutoff1':0.1,'cutoff2':0.3}

continuous={'bandname':'continuous','cutoff1':0.0,'cutoff2':0.1}


# Selection of the signal to analyse
linescanlist=['rigid lumen', 'lumen','rigid endfoot','PVS', 'endfoot','Area']

# Selection of the frequency bands to analyse
bandfreqencies=[verylowfreq,lowfreq,resp,cardiac]




#################################################
###  Script
#################################################


# import data
outputdirdatabase='../output/databases/'+name+'/'


# get the list of files
files = [f for f in listdir(dir) if isfile(join(dir, f))]

studylist=[]

for f in files :

    split=re.split(r'\s|\.',f)
    mouse=split.pop(0)
    number=split.pop(0)
    date=split.pop(0)
    trial=split.pop(0)
    split.pop(-1)
    nofinal=split.pop(-1)
    vessel=' '.join(split)
    path=dir+f
    
    studylist.append({'trial':trial,'mouse':mouse, 'number':number, 'vessel':vessel, 'vesselID':nofinal, 'file':path}
)
    


for study in studylist :
    plt.close('all')
    
    mousenumber=study['number']
    mousekind=study['mouse']
    trial=study['trial']
    vesseltype=study['vessel']
    file=study['file']
    nofinal=study['vesselID']
    
    print('Analysis of :',mousekind+' '+mousenumber+' '+trial+' '+vesseltype+' '+nofinal)
       
    data=pd.read_csv(file, decimal='.', delimiter=",")
    
    data['file']=file
    
    
    data['lumen']=-(data['lumen_upper']-data['lumen_lower'])/2# equiv to radius
    data['rigid lumen']=(data['lumen_upper']+data['lumen_lower'])/2
 
    try :
        data['endfoot']=-(data['endfoot_upper']-data['endfoot_lower'])/2# equiv to radius
        data['rigid endfoot']=(data['endfoot_upper']+data['endfoot_lower'])/2
    except :
        data['endfoot']=np.nan
        data['rigid endfoot']=np.nan
    
    try :
        data['PVS']=-(data['lumen']-data['endfoot']) 
        data['Area']=(np.pi*data['endfoot']**2-np.pi*data['lumen']**2)
    except :
        data['PVS']=np.nan
        data['Area']=np.nan
  
    
    #Smooth the signal (remove noise artifact but keep the cardiac details)
    for linescan in linescanlist :
        data[linescan]=savgol_filter(data[linescan],11,3)  
  
    coupledtraces={}
    
    coupledtraces['lumen']=[data['lumen_upper'],data['lumen_lower']]
    try :
        coupledtraces['endfoot']=[data['endfoot_upper'],data['endfoot_lower']]
    except :
        coupledtraces['endfoot']=[data['endfoot'],data['endfoot']]
    


    # get sequences
    REM_list, NREM_list, IS_list, BASE_list, LOCO_list,WHISK_list,QUIET_list,AWAKENING_list=da.extract_sequences(data)

    outputdir='../output/amp_analysis/'+name+'/'+mousekind+'-'+mousenumber+'-'+trial+'-'+nofinal+'/'
    
    # Directory to store output data             
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
        
    

    for linescan in linescanlist :

                
        # Global spectrum
        signal=data[linescan][:].values
        
        signal = signal[~np.isnan(signal)]
        
        time=data['t'][:].values
        time_step=time[1]-time[0]
        fs=1/time_step
    
        frequency, Pxx_spec = sg.periodogram(signal,fs, 'hanning', scaling='spectrum')
    
        fig, ax = plt.subplots()
        ax.plot(frequency,Pxx_spec)
        
        maxlim=Pxx_spec.max()
        
        # Create a Rectangle patch
        rect = Rectangle((cardiac['cutoff1'],0 ),cardiac['cutoff2'] - cardiac['cutoff1'], maxlim, facecolor='r', alpha=0.4)
        # Add the patch to the Axes
        ax.add_patch(rect)
        
        
        # Create a Rectangle patch
        rect = Rectangle((resp['cutoff1'],0 ), resp['cutoff2'] - resp['cutoff1'], maxlim, facecolor='g', alpha=0.4)
        # Add the patch to the Axes
        ax.add_patch(rect)
        
        
        # Create a Rectangle patch
        rect = Rectangle((lowfreq['cutoff1'],0 ), lowfreq['cutoff2'] - lowfreq['cutoff1'], maxlim, facecolor='b', alpha=0.4)
        # Add the patch to the Axes
        ax.add_patch(rect)
        
        # Create a Rectangle patch
        rect = Rectangle((verylowfreq['cutoff1'],0 ), verylowfreq['cutoff2'] - verylowfreq['cutoff1'], maxlim, facecolor='m', alpha=0.4)
        # Add the patch to the Axes
        ax.add_patch(rect)
        
        # Create a Rectangle patch
        rect = Rectangle((continuous['cutoff1'],0 ), continuous['cutoff2'] - continuous['cutoff1'], maxlim, facecolor='k', alpha=0.4)
        # Add the patch to the Axes
        ax.add_patch(rect)
        
        plt.xscale('log')
        plt.yscale('log')
        
        plt.savefig(outputdir+'/full-spectrum'+outputformat)
        plt.close()
              
    
    
    for stage, seqlist in zip(['IS','NREM','REM','Baseline','Locomotion','Whisking','Quiet','Awakening'],[IS_list,NREM_list,REM_list,BASE_list,LOCO_list,WHISK_list,QUIET_list,AWAKENING_list]):
        
        namedatabase=outputdirdatabase+'detailled/'+name+mousekind+'-'+mousenumber+'-'+trial+'-'+nofinal+'-'+stage+'.pkl'
        
        if os.path.isfile(namedatabase) :
            continue  

        if len(seqlist)==0:
            continue
        
        # Analysis
        ampdata_list=[]    
        averagedata_list=[]
        
        averageamp={}
        averagemean={}
        averageperiod={}
        averagecorrlumen={}
        averagecorrendfoot={}
        
        for FB in bandfreqencies :
            averageamp[FB['bandname']]={}
            averageperiod[FB['bandname']]={}
            
        
        
        seqno=0
        for sequence in seqlist:
            
            for linescan in linescanlist :

                print('Analysis of '+linescan)
                
                if np.isnan(data[linescan][0]) :
                    continue
                

                #signal
                signal=data[linescan][sequence.ibegin:sequence.iend].values

                
                #check for nan values
                nans, indexes= da.nan_helper(signal)
                
                #if there is more than 50% of missing data we skip
                if sum(nans)/len(signal)>=0.5 :
                    continue
                
                # remove the nan values at the begining or the end of the sequence
                ibegin=sequence.ibegin+np.where(~nans)[0][0]
                iend=sequence.ibegin+np.where(~nans)[0][-1]
                
                traces={}
                
                traces['lumen1']=coupledtraces['lumen'][0][ibegin:iend].values
                traces['lumen2']=coupledtraces['lumen'][1][ibegin:iend].values
                traces['endfoot1']=coupledtraces['endfoot'][0][ibegin:iend].values
                traces['endfoot2']=coupledtraces['endfoot'][1][ibegin:iend].values
                
                vessel=data['lumen'][ibegin:iend].values
                
                signal=data[linescan][ibegin:iend].values
                
                
                averagemean[linescan]=np.median(signal)
                
                #check for nan values
                for tracelabel in traces :
                    if sum(np.isnan(traces[tracelabel]))<len(traces[tracelabel]) :
                        traces[tracelabel],nans=da.remove_nan(traces[tracelabel])
                              
                vessel,nans=da.remove_nan(vessel)
                signal,nans=da.remove_nan(signal)
                

                time=data['t'][ibegin:iend].values-data['t'][ibegin]
                
                time_step=time[1]-time[0]
                fs=1/time_step
                
                frequency, Pxx_spec = sg.periodogram(signal,fs, 'hanning', scaling='density')
                
                
                fig, ax = plt.subplots()
                ax.plot(frequency,Pxx_spec)
                #plt.yscale('log')
                #plt.xlim([0,12])
                #plt.ylim([0,0.1])
        
                # Create a Rectangle patch
                rect = Rectangle((cardiac['cutoff1'],0 ),cardiac['cutoff2'] - cardiac['cutoff1'], 100, facecolor='r', alpha=0.4)
                # Add the patch to the Axes
                ax.add_patch(rect)


                # Create a Rectangle patch
                rect = Rectangle((resp['cutoff1'],0 ), resp['cutoff2'] - resp['cutoff1'], 100, facecolor='g', alpha=0.4)
                # Add the patch to the Axes
                ax.add_patch(rect)

                
                # Create a Rectangle patch
                rect = Rectangle((lowfreq['cutoff1'],0 ), lowfreq['cutoff2'] - lowfreq['cutoff1'], 100, facecolor='b', alpha=0.4)
                # Add the patch to the Axes
                ax.add_patch(rect)
                
                # Create a Rectangle patch
                rect = Rectangle((verylowfreq['cutoff1'],0 ), verylowfreq['cutoff2'] - verylowfreq['cutoff1'], 100, facecolor='m', alpha=0.4)
                # Add the patch to the Axes
                ax.add_patch(rect)
                
                # Create a Rectangle patch
                rect = Rectangle((continuous['cutoff1'],0 ), continuous['cutoff2'] - continuous['cutoff1'], 100, facecolor='k', alpha=0.4)
                # Add the patch to the Axes
                ax.add_patch(rect)
                
                plt.xscale('log')
                plt.yscale('log')
                
                #plt.xlim([0,12])
                #plt.ylim([0,100])
                
                plt.savefig(outputdir+'/spectrum_'+linescan+'-'+stage+str(seqno)+''+outputformat)
                plt.close()
                
                #peak to peak analysis
                # in each frequency band

                fig, axs = plt.subplots(6,1,figsize=(16,9), gridspec_kw={'height_ratios': [1,1,1, 1,1,1]})
                axs[0].plot(time,signal)


                signal_continuous=da.bandpassfilter(signal-np.mean(signal), fs, continuous['cutoff1'], continuous['cutoff2'])
                
                axs[1].plot(time,signal_continuous)
                
                ifreq=0
                
                reconstruction=np.mean(signal)+signal_continuous
                
                plt.rcParams["axes.linewidth"] = 2
                plt.rcParams["patch.linewidth"] = 2
                
                for FB in bandfreqencies :

                #if sequence.end-sequence.begin > 2/cutoff2 :
                    
                    print(stage+'-'+str(seqno)+'-'+FB['bandname'])
                    
                    signal_filtered=da.bandpassfilter(signal-np.mean(signal), fs, FB['cutoff1'], FB['cutoff2'])
                    
                    reconstruction+=signal_filtered
                    
                    axs[ifreq+2].plot(time,signal_filtered)
                    ifreq+=1
                axs[0].plot(time,reconstruction)   
                
                axs[1].set_ylim([-0.8,0.8])
                axs[2].set_ylim([-0.8,0.8])
                axs[3].set_ylim([-0.8,0.8])
                axs[4].set_ylim([-0.15,0.15])
                axs[5].set_ylim([-0.15,0.15])
                
                for ax in axs:
                     ax.set_xlim([0,15])
                
                for ax in axs[0:-1]:
                    ax.set_xticks([])
                   
                plt.savefig(outputdir+'/decomposition_'+linescan+'-'+stage+str(seqno)+''+outputformat)
                plt.close()
                
                for FB in bandfreqencies :
                    signal_filtered=da.bandpassfilter(signal-np.mean(signal), fs, FB['cutoff1'], FB['cutoff2'])
                    
                    # keep only cardiac FB for the traces
                    for tracelabel in traces :
                        traces[tracelabel]=da.highpassfilter(traces[tracelabel]-np.mean(traces[tracelabel]), fs, cardiac['cutoff1'])

                    
                    
                    #ampl,period,tvalley,tpeak,mean=da.amp_analysis(np.mean(signal)+signal_continuous,signal_filtered,time,FB['cutoff1'], FB['cutoff2'],export=linescan+'-'+stage+str(seqno)+'-'+FB['bandname'],outputdir=outputdir)
                    ampl,period,tvalley,tpeak,mean=da.amp_analysis(signal,signal_filtered,time,FB['cutoff1'], FB['cutoff2'],export=linescan+'-'+stage+str(seqno)+'-'+FB['bandname'],outputdir=outputdir, outputformat=outputformat)

                    
                    if FB['bandname']=='cardiac':
                        amplt=np.sqrt(ampl)
                        periodt=np.log(period)
                        outliers=(amplt<=np.mean(amplt)-2*np.std(amplt))|(amplt>=np.mean(amplt)+2*np.std(amplt))
                        outliers=outliers|(periodt<=np.mean(periodt)-2*np.std(periodt))|(periodt>=np.mean(periodt)+2*np.std(periodt))

                    else : 
                        outliers=np.ones_like(ampl)*np.nan


                    
                    if len(mean)==0 :
                        continue
                    else : 
                        
                        # The data might be too noisy so we measure metrics of reliability
                        corr_lumen=[]
                        corr_endfoot=[]
                        has_nan=[]
                        median_lumen=[]
                        for i in range(0,len(ampl)) :
                            filter=(time<=tvalley[i]+period[i])&(time>=tvalley[i])
                            
                            has_nan.append(nans[filter].max())
                            median_lumen.append(np.median(vessel[filter]))
                            r_lumen,p=pearsonr(traces['lumen1'][filter],traces['lumen2'][filter])
                            try :
                                r_endfoot,p=pearsonr(traces['endfoot1'][filter],traces['endfoot2'][filter])
                            except :
                                r_endfoot=np.nan
                            corr_lumen.append(r_lumen)
                            corr_endfoot.append(r_endfoot)
                            
                            
                        ### tortuosity and variability of the highpass filtered signal
                        cuttedsignal=da.highpassfilter(signal-np.mean(signal), fs, 4)

                        tortuosity=[]
                        rawvar=[]
                        for i in range(0,len(ampl)) :
                            filter=(time<=tvalley[i]+period[i])&(time>=tvalley[i])
                            tortuosity.append(da.arc_length(time[filter], cuttedsignal[filter]))
                            rawvar.append(cuttedsignal[filter].std()) 
                            
                            
                        if FB['bandname']=='cardiac':
                            varamp=[]
                            varperiod=[]
                            slopevar=[]
                            for i in range(0,len(ampl)) :
                                #measure the variability of the cardiac pulsation over 1 second
                                filter=(tvalley<=tvalley[i]+0.5)&(tvalley>=tvalley[i]-0.5)
                                slope=ampl[filter]/(tpeak[filter]-tvalley[filter])
                                slopevar.append(slope.std())
                                varamp.append(ampl[filter].std())
                                varperiod.append(period[filter].std())
                        else :
                            varamp=np.ones_like(ampl)*np.nan
                            varperiod=np.ones_like(ampl)*np.nan
                            slopevar=np.ones_like(ampl)*np.nan
                            
                           
                        if (np.array(mean)<0).any():
                            print('problem : there is at least a negative value in the mean radius !!!')
                            
                        #if ((np.array(ampl)/np.array(mean)).max()>0.2)&(FB['bandname']=='cardiac') :
                        #    print('oups')
                            
    
                        #store data    
                        for a,p,m,t1,t2,cl,ce,ampvar,periodvar,tor,var,svar,out,nan,ml in zip(ampl,period,mean,tvalley,tpeak,corr_lumen,corr_endfoot,varamp,varperiod,tortuosity,rawvar,slopevar,outliers,has_nan,median_lumen) :
                            dictamp={'file':file,'mouse kind':mousekind,'mouse number':mousenumber,'trial':trial, 'vessel':vesseltype, 'vesselID':nofinal,'linescan':linescan,'cutoff1':FB['cutoff1'],'cutoff2':FB['cutoff2'],'sequence':seqno,'stage':stage,'tini':sequence.begin,'tend':sequence.end, 'bandname':FB['bandname']}
                            dictamp['amp']=a
                            dictamp['period']=p
                            dictamp['corrlumen']=cl
                            dictamp['correndfoot']=ce
                            dictamp['periodvar']=ampvar
                            dictamp['ampvar']=periodvar
                            dictamp['slopevar']=svar                            
                            dictamp['tortuosity']=tor
                            dictamp['rawvar']=var
                            dictamp['mean']=m
                            dictamp['outlier']=out
                            dictamp['tmin']=t1+data['t'][ibegin]
                            dictamp['tmax']=t2+data['t'][ibegin]
                            dictamp['medianepisode']=np.median(signal)
                            dictamp['hasnan']=nan
                            dictamp['medianlumen']=ml
                            ampdata_list.append(dictamp)
                            
                        #compute the median values over the episode for a more concise database
                        averageamp[FB['bandname']][linescan]=np.median(ampl)
                        averageperiod[FB['bandname']][linescan]=np.median(period)
                        averagecorrendfoot[FB['bandname']]=np.median(corr_lumen)
                        averagecorrlumen[FB['bandname']]=np.median(corr_endfoot)
                            
                        

            # store the averages values over the sequence
            dictaverage={'file':file,'mouse kind':mousekind,'mouse number':mousenumber,'trial':trial, 'vessel':vesseltype, 'vesselID':nofinal,'sequence':seqno,'stage':stage,'tini':sequence.begin,'tend':sequence.end}
            for linescan in linescanlist:
                for FB in bandfreqencies :
                    try :
                        dictaverage['amp '+linescan+' '+FB['bandname']]=averageamp[FB['bandname']][linescan]
                    except :
                        dictaverage['amp '+linescan+' '+FB['bandname']]=np.nan
                        
                    try :    
                        dictaverage['period '+linescan+' '+FB['bandname']]=averageperiod[FB['bandname']][linescan]
                    except :
                        dictaverage['period '+linescan+' '+FB['bandname']]=np.nan
                        
                    try :
                        dictaverage['corrlumen '+FB['bandname']]=averagecorrlumen[FB['bandname']]
                    except :
                        dictaverage['corrlumen '+FB['bandname']]=np.nan
                        
                    try :
                        dictaverage['correndfoot '+FB['bandname']]=averagecorrendfoot[FB['bandname']]
                    except :
                        dictaverage['correndfoot '+FB['bandname']]=np.nan
                try :        
                    dictaverage['mean '+linescan]=averagemean[linescan]
                except :
                    dictaverage['mean '+linescan]=np.nan
                    
            averagedata_list.append(dictaverage)
            
            seqno+=1
        
        amplitudedb = pd.DataFrame(ampdata_list)
        averagedb = pd.DataFrame(averagedata_list)
        
        # ### For normalisation we compute meanvalues during the baseline stage
        
        # amplitudedb['median base radius']=np.nan
        # amplitudedb['median base period']=np.nan
        

        # # The mean radius is taken over the time scales of the LF oscillations
        # filter=((amplitudedb['stage']=='Baseline')|(amplitudedb['stage']=='Quiet'))
        # meanradius=np.median(amplitudedb['medianepisode'][filter]) 
        # amplitudedb['median base radius']=meanradius
           
        # # for each frequency band we compute the mean period during the baseline stage
        # for FB in bandfreqencies :
        #     bandname=FB['bandname']
        #     filter=((amplitudedb['stage']=='Baseline')|(amplitudedb['stage']=='Quiet'))&(amplitudedb['bandname']==bandname)
        #     meanperiod=np.median(amplitudedb['period'][filter])
        #     amplitudedb['median base period'][(amplitudedb['bandname']==bandname)]=meanperiod

        # # Normalised values compared to the baseline stage    
        # amplitudedb['amp norm']=amplitudedb['amp']/amplitudedb['median base radius']
        # amplitudedb['mean norm']=amplitudedb['mean']/amplitudedb['median base radius']
        # amplitudedb['period norm']=amplitudedb['period']/amplitudedb['median base period']
            
            
        # Directory to store output data
        if not os.path.exists(outputdirdatabase):
            os.makedirs(outputdirdatabase)  

        if not os.path.exists(outputdirdatabase+'detailled/'):
            os.makedirs(outputdirdatabase+'detailled/')  
        
        if not os.path.exists(outputdirdatabase+'averaged/'):
            os.makedirs(outputdirdatabase+'averaged/')              
            
        
        print('save:'+name+mousekind+'-'+mousenumber+'-'+trial+'-'+nofinal+'-'+stage+'.pkl')
        amplitudedb.to_pickle(outputdirdatabase+'detailled/'+name+mousekind+'-'+mousenumber+'-'+trial+'-'+nofinal+'-'+stage+'.pkl')
        averagedb.to_pickle(outputdirdatabase+'averaged/'+name+mousekind+'-'+mousenumber+'-'+trial+'-'+nofinal+'-'+stage+'.pkl')

