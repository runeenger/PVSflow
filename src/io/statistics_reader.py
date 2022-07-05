#! /usr/bin/env python3
################################
### Module to read the statistics analysis output files
##################################


import numpy as np

def ReadFixedEffect(file):
    """ Read the mean estimates from a file containing the statistical analysis of the peak to peak data.
        input : 
            - file : a string providing the file where the information must be read 
                    (one file corresponds to a type of mouse, type of vessel, type of measure, frequency band selected)

        output :  a dict. The keys corresponds to the stages. The values are float of the estimates.
                For the stage 'baseline', the value is the Intercept value.
                After discussion with Rune, this is not used anymore : If the stage is not baseline, then the value is a different value from baseline only if the p value for this stage is < 0.05
                meaning that this stage is significantly different from the baseline.
                
    """
    

    d={}
    with open(file) as f:
        line=f.readline()# skip first line
        # read data
        line=f.readline()
        lines=[]
        while line:
            line=line.split('\n')[0]
            lines.append(line.split(' '))
            line = f.readline()
    
    #get the intecept value (correspond to the log of the variable)
    intercept=float(lines[0][1])

    #the baseline value correspond to the exp of the intercept
    d['baseline']=np.exp(intercept)

    # the other stages we add the estimate to the intecept.   NOT USED ANYMORE : only if significative difference with baseline
    for line in lines[1::] :
        d[line[0].replace('"', '')]=np.exp(intercept+float(line[1]))#*(float(line[3])<0.05))

    return d


def ReadRandomEffect(file):
    """ Read the mean estimates from a file containing the statistical analysis of the peak to peak data.
        input : 
            - file : a string providing the file where the information must be read 
                    (one file corresponds to a type of mouse, type of vessel, type of measure, frequency band selected)

        output :  a dict. The keys corresponds to the stages. 
        The values are np arrays of the estimates (floats) for each vessel.

    """

    d={}
    with open(file) as f:
        # first line contains the name of the stages
        stages=f.readline()
        stages=stages.split('\n')[0]
        stages=stages.split(' ')
        # read data
        line=f.readline()
        lines=[]
        while line:
            line=line.split('\n')[0] 
            lines.append(line.split(' '))
            line = f.readline()
    

    # intercepts for each vessel
    d['baseline']={}
    for line in lines:
        vesselid=line[0]
        vesselid=vesselid.replace('"', '')
        vesselid=vesselid.replace('.', '-')
        #get the intecept value (correspond to the log of the variable)
        intercept=float(line[1])
        d['baseline'][vesselid]=np.exp(intercept)
        #the baseline value correspond to the exp of the intercepts
    
    for columnno,stage in enumerate(stages[1::]) :
        stage=stage.replace('"', '')
        d[stage]={}
        #pvalues=[]
        for line in lines:
            #get the estimates value (correspond to the log of the variable)
            intercept=float(line[1])
            estimate=line[columnno+2]# first two columns are for vessel id and intercept
            vesselid=line[0]
            vesselid=vesselid.replace('"', '')
            vesselid=vesselid.replace('.', '-')
            try :
                estimate=float(estimate)
            except :
                estimate=0.
                print('aie aie aie')
                
            ### WARNING : The random effect files do not need estimate to be added to intercept !!!
            d[stage][vesselid]=np.exp(estimate)

    

    return d