#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 14:12:40 2021

@author: alexandra


This script convert the pickle database (python format) to a CSV format.
It will clean the data before export.
"""

import pickle
import pandas as pd
import os
from os import listdir
from os.path import isfile, join

import numpy as np

##########################################################
### Parameters of the script : can be modified
##########################################################

# read all the databases in the folder :
dir='../output/databases/Veins/detailled/'
# name of the output file
name='VeinsWT10.csv'



#################################################
###  Script
#################################################

files = [f for f in listdir(dir) if isfile(join(dir, f))]
    
studylist=[]
    
for f in files : 
        data=pd.read_pickle(dir+f)
        studylist.append(data)
        

amplitudedb=pd.concat([data for data in studylist], ignore_index=True)


# filter the data

filter=(amplitudedb['amp']>0) # keep only positive amplitudes
filter=filter&(amplitudedb['mean']>0) # keep only positive mean values
filter=filter&(amplitudedb['hasnan']==False) # remove nan data
# remove outliers
filter=filter& np.logical_not((amplitudedb['bandname']=='cardiac')&(amplitudedb['outlier']))


amplitudedb=amplitudedb[filter]

print(np.unique(amplitudedb[['linescan']]))
print(np.unique(amplitudedb[['stage']]))


# export 

# Directory to store output data
if not os.path.exists('../output/'):
        os.makedirs('../output/')  

if not os.path.exists('../output/'+'csv/'):
        os.makedirs('../output/'+'csv/') 


amplitudedb.to_csv('../output/csv/'+name, index=False)