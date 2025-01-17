#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 11:15:51 2021

@author: alexandra

compute the apparent diffusion coefficient from the simulated profile of concentration in the PVS.

"""

import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, exp, log, pi

import os
import shutil

from os import path

from src.toolbox import U

from scipy.interpolate import UnivariateSpline
from sklearn import linear_model
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error, r2_score

import warnings
warnings.simplefilter(action='ignore', category=RuntimeWarning)
#todo : do not plot in the interface. Just save the figure.

def FWHM(X,Y):
    spline = UnivariateSpline(X, Y-np.max(Y)/2, s=0)
    roots=spline.roots() # find the roots
    if len(roots)==1 :
        #print('Warning : there is only one root to the spline, we return 2*dist(max)')
        fwhm=2*abs(X[np.argmax(Y)]-roots[0])
    elif len(roots)==2 :
        fwhm=abs(roots[1]-roots[0])
    else :
        #print('Error : unexpected number of roots of the fitted spline')
        fwhm=X.max()-X.min()
    
    return fwhm


# With only one linear fit 
def estimate_diff_FWHM_fit(spanTime,spanFWHM) :
    spanSigma=spanFWHM/(2*sqrt(2*log(2)))
    
    X = np.array(spanTime).reshape((-1, 1))
    y = np.array(spanSigma**2)/2

    regressor = linear_model.LinearRegression()
    regressor.fit(X, y)
    a = regressor.coef_[0]
    b = regressor.intercept_
    plt.figure()
    plt.plot(spanTime,spanSigma**2/2,'*')
    plt.plot(spanTime,a*(spanTime-(-b/a)))
    plt.xlabel('time (s)')
    plt.ylabel(u'$\sigma^2/2$')
    plt.title('linear fit for the FWHM method')

    return regressor.coef_[0]

# with the slope, can vary with time
def estimate_diff_FWHM_slope(spanTime,spanFWHM) :
    spanSigma=spanFWHM/(2*sqrt(2*log(2)))
    y = np.array(spanSigma**2)/2

    slope=np.diff(y)/np.diff(spanTime)
    
    plt.plot(spanTime[0:-1],slope,'*')
    plt.xlabel('time (s)')
    plt.ylabel('$Deff$')
    plt.title('Slope of \sigma^2/2 vs t')
    
    return slope

def gaussian (x,t,s,D,L,xi) :
    return s/sqrt(2*D*t)*np.exp(-(x-xi)**2/(4*D*t))


def estimate_diff_fit(spanC,spanX,t,sigma,L) :
    
    
    def fitfunction(x,D,xi) :
        tini=sigma**2/(2*D)
        return  gaussian (spanX,t+tini,sigma,D,L,xi)
    
    popt, pcov = curve_fit(fitfunction, spanX,spanC,bounds=([1e-9,0],[1e-4,L]))
    
    fit = fitfunction(spanX,*popt)
    
    err=mean_squared_error(spanC, fit) 
    
    D=popt[0]
    xi=popt[1]
        
    return D,xi, err

##########################################################
### Parameters of the script : can be modified
##########################################################

# localisation of the simulation folder
rep='../data/simulations/'
analysis='-d7e-08-l6e-02'  # chose different length or diffusion coefficient if needed
#rep=rep+'dispersionRandomWT10t40area'+analysis+'/'
#rep=rep+'dispersionRandomveinsWT10'+analysis+'/'
rep=rep+'dispersionSMC50WT10'+analysis+'/'
# name of the output file
#outputname='RandomVeinsWT10.csv'
#outputname='RandomWT10.csv'
outputname='SMC50WT10new.csv'

# set a condition on time analysis (to avoid analysis when tracer reaches the boundary conditions)
conditiontime=False
conditionalpha=False

## stages and frequency bands to analyse
stages=['baseline','stageIS','stageNREM','stageREM','stageAwakening']
#bandnames=['card-v1e-03','card-v5e-03', 'card-v1e-02','card','resp','LF','VLF']
bandnames=['LF','VLF']
bandnames=['card-v1e-03','card-v5e-03', 'card-v1e-02','LF','VLF']

##########################################################
### script
##########################################################

#Creation of the database
Database=[]
datalabel=['job','stage','bandname','vesselID', 'Rv0', 'Rpvs', 'L', 'DX', 'dt', 'rho', 'mu', 'D', 'sigma', 'xi', 'f', 'umax','pmax', 'Pe', 'Re', 'Wo', 'Fo' , 'A', 'beta', 'dPdx', 'T','nPeriod' , 'tend', 'FWHMend','DestFWHM', 'Destfit', 'RFWHM', 'Rfit', 'amp','thetaa', 'tau']

#### Analysis

file='concentration.txt'


# get all the log files
from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir(rep) if isfile(join(rep, f))]
islog=[string.endswith('.log') for string in onlyfiles]
whereislog=np.where(islog)[0]
logfiles=[onlyfiles[index] for index in whereislog]

import re
end=[re.findall(r'id\S+PVSinfo.log',string) for string in logfiles]

vesselIDS=[re.findall(r'id\S+PVSinfo.log',string)[0].replace('_PVSinfo.log','') for string in logfiles]

vesselIDS=np.unique(vesselIDS)




errorfiles=[]#

for stage in stages :
    for bandname in bandnames :
        for vessel in vesselIDS :
            
            job='disp'+analysis+'-'+stage+'-'+bandname+'-'+vessel
            print(job)
            print(rep)
            
            file=job+'_concentration.txt'

            try :
                Data=np.loadtxt(rep+file,delimiter=',')

                t=Data[1:,0]
                
                t[1]
                
            except :
                print('Error during the loading of '+file)
                errorfiles.append([file])
                continue
            
            
            print('\n')
            print('*** Analysis of :'+job)
            
            
            outputfolder='../output/disp_analysis/'+job+'/'
            
            if not os.path.exists(outputfolder):
                os.makedirs(outputfolder)
                
            import glob
            
            files = glob.glob('outputfolder'+'*')
            for f in files:
                os.remove(f)
            
            
            
            outputfile = open(outputfolder+'post-process.txt', "w")
            
            outputfile.write('#'*20)
            outputfile.write('\n# Dispersion analysis of '+job)
            outputfile.write('\n'+'#'*20)
            

            
            
            

            concentration=Data[1:,1:]

            pressure=np.loadtxt(rep+job+'_pressure.txt',delimiter=',')[1:,1:]

            velocity=np.loadtxt(rep+job+'_velocity.txt',delimiter=',')[1:,1:]
            
            #### Simulation parameters
            
            import re
            scinot = re.compile('[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)')
            
            # Should get D, L, rho, nu, f, amp, Rpvs, Rv0,sigma and x0 from the log file
            file=job+'_PVSinfo.log'
            
            
            with open(rep+file) as fl:
                line = fl.readline()
                while line:
                    line = fl.readline()
                    if line[0:6] == 'Vessel' :
                        Rv0=float(re.findall(scinot, line)[0])
                    elif line[0:6] == 'PVS ra' :
                        Rpvs=float(re.findall(scinot, line)[0])  
                    elif line[0:6] == 'PVS le' :
                        L=float(re.findall(scinot, line)[0])   
                    elif line[0:6] == 'densit' :
                        rho=float(re.findall(scinot, line)[0])
                    elif line[0:6] == 'dynami' :
                        mu=float(re.findall(scinot, line)[0])   
                    elif line[0:6] == 'Free d' :
                        D=float(re.findall(scinot, line)[0]) 
                    elif line[0:6] == 'STD of' :
                        sigma=float(re.findall(scinot, line)[0])  
                    elif line[0:6] == 'Initia' :
                        xi=float(re.findall(scinot, line)[0])
                    elif line[0:6] == 'fi (Hz':
                        f=float(re.findall(scinot, line)[0]) # Check if it works for several frequencies
                    elif line[0:6] == 'ai (di':
                        amp=float(re.findall(scinot, line)[0])*100 # Check if it works for several frequencies
                    elif line[0:6] == 'time s':
                        time_step=float(re.findall(scinot, line)[0])
                    elif line[0:6] == 'cell s':
                        cell_size=float(re.findall(scinot, line)[0])       
                        
            outputfile.write('\n\nRv : %.2f um'%(Rv0*1e4))
            outputfile.write('\nRpvs : %.2f um'%(Rpvs*1e4))
            outputfile.write('\nLength : %.2f um'%(L*1e4))
            
            x=np.linspace(0,L,len(concentration[0,:]))
            
            outputfile.write('\n\nDensity : %.2e g/cm3'%rho)
            outputfile.write('\nViscosity : %.2e dyn s /cm2'%mu)
            nu=mu/rho
            
            outputfile.write('\n\nDiffusion coefficient: %.2e cm2/s'%(D))      
            outputfile.write('\nSTD gausian : %.2e cm'%(sigma))
            outputfile.write('\ncenter gaussian: %.2e um'%(xi*1e4))
            
            outputfile.write('\n\nfrequency : %.2e Hz'%f)
            outputfile.write('\n\namplitude : %.2e pc'%amp)
            
            outputfile.write('\nspatial resolution : %.2e um'%(cell_size*1e4))
            outputfile.write('\ntemporal resolution : %.2e s'%time_step)
            
            ### Compute dimensionless numbers
            
            umax=np.max(abs(velocity))
            pmax=np.max(abs(pressure))
            
            
            w=2*pi*f
            h=Rpvs-Rv0
            
            Pe=h*umax/D
            Re=rho*umax*h/mu
            Wo=h*sqrt(w/mu)
            Fo=D*t[-1]/xi**2
            
            outputfile.write('\nUmax : %.3f um/s'%(umax*1e4))
            outputfile.write('/npmax : %.3f pa \n'%(pmax/10))
            
            outputfile.write('\nmax Reynolds number : %.0e'%Re)
            outputfile.write('\nmax Peclet number : %.0e'%Pe)
            outputfile.write('\nWomersley number : %.0e'%Wo)
            outputfile.write('\nFourier number : %.0e'%Fo)
            
            dPdx=pmax/L
            beta=D/mu
            A=pi*(Rpvs**2-Rv0**2)
            
            
            
            # We want to look to the results at each period (net flow 0)
            dtoutput=t[1]-t[0]
            tend=t[-1]
            
            
            
            if tend>(1/f):
                T=1/f
                ishift=0 #int(T/dtoutput*1/4)
            else :
                T=t[-1]/10
                ishift=0
                
            outputfile.write('\nfinal time simulation : %.2e s'%tend)
            outputfile.write('\noutput period : %.2e s'%dtoutput)
            outputfile.write('/nperiod : %.2e s'%T)
            DX=(x[1]-x[0])*1e4
            outputfile.write('\nspatial resotion : %.2e um'%(DX))
            
            
            iperiodic=(np.arange(ishift,len(t),max(1,round(T/dtoutput)))).astype(int)
            iperiodicumax=(np.arange(round(T/dtoutput*1/2),len(t),max(1,round(T/dtoutput)))).astype(int)
            

            plt.figure()
            plt.plot(np.max(velocity,axis=1))
            plt.plot(iperiodic,np.max(velocity[iperiodic],axis=1),'*r')
            
            xiadv=x[np.argmax(concentration[iperiodicumax[0]])]
           
            
            span_FWHM=[]
            for c in concentration :
                span_FWHM.append(FWHM(x,c))
                #plt.plot(x,c)
                
            span_FWHM=np.array(span_FWHM) 
            #plt.xlim([35e-4,65e-4])
            
             
            ### Estimation of D from FWHM
            plt.figure()
            Dest=estimate_diff_FWHM_fit(t[iperiodic[::]],span_FWHM[iperiodic[::]])
            DestFWHM=Dest
            
            DestFWHM=max(Dest,D)
            
            tmax=(L/4)**2/DestFWHM/2
            
            if conditiontime :
                #condition diffusion
                while t[iperiodic[-1]]>tmax : 
                    outputfile.write('\n* warning : limitation of tend due to diffusion')
                    it=np.where(t>=tmax)[0][0]
                    ii=np.where(iperiodic>=it)[0][0]
                    iperiodic=iperiodic[0:ii]
  
                    plt.figure()
                    Dest=estimate_diff_FWHM_fit(t[iperiodic[::]],span_FWHM[iperiodic[::]])
                    DestFWHM=max(Dest,D)
                    tmax=(L/2)**2/DestFWHM/2
                    

            thetaa=(xiadv+np.sqrt(2*DestFWHM*t[iperiodic]))/L
            
            #amp*L/(L/2-2*np.sqrt(2*DestFWHM*t[iperiodic]))
            
            if conditionalpha :
                #condition advection diffusion
                while thetaa[-1]>=0.8 : 
                    outputfile.write('\n* warning : limitation of tend due to theta a')
                    ii=np.where(thetaa>=0.8)[0][0]

                    if ii==0 :
                        iperiodic=[0]
                        break
                    
                    iperiodic=iperiodic[0:ii]
                    Dest=estimate_diff_FWHM_fit(t[iperiodic[::]],span_FWHM[iperiodic[::]])
                    DestFWHM=max(Dest,D)
                    #thetaa=amp*L/(L/2-2*np.sqrt(2*DestFWHM*t[iperiodic]))
                    thetaa=(xiadv+np.sqrt(2*DestFWHM*t[iperiodic]))/L

            nPeriod=len(iperiodic)-1
            outputfile.write('\nnumber of period analysed: %i'%nPeriod)
            
            if nPeriod==0 :
                outputfile.write('\nOscillatory dispersion analysis aborted')
                continue
            

            plt.savefig(outputfolder+'disp_FWHM.png')                 
            plt.close()
            
            outputfile.write('\n\nfinal time analysis : %.2e s'%t[iperiodic[-1]])
            
            
            

            plt.figure()
            plt.plot(t[iperiodic],span_FWHM[iperiodic]*1e4,'*')
            plt.xlabel('time (s)')
            plt.ylabel('FWHM (um)')
            plt.title('computed FWHM')
            
            plt.savefig(outputfolder+'FWHM.png')
            
            tau=t[iperiodic[-1]]/(h**2/D)
            outputfile.write('\ntau: %.2e'%tau)
            outputfile.write('\ntheta a : %.1e'%thetaa[-1])
            
            outputfile.write('\n\nEstimation of D from FWHM : %.2e'%DestFWHM)
            
            
            
            #When the apparent Diffusion coefficient is smaller with oscillation than diffusion alone at small time. What is going on? 
            
            ### Estimation of D from fit
            
            # we take the time at last period
            plt.figure()
            spanDest=[]
            for itime in iperiodic[1::] : 
                #print('\nEstimate at time %f s'%t[itime])
            
                Dest, Xi, err = estimate_diff_fit(concentration[itime],x,t[itime],sigma,L)
            
                #print('Estimation of D from fit : %.2e'%Dest)
                #print('Fit error : %.2e'%err)
                spanDest.append(Dest)
                plt.plot(t[itime],Dest,'*')
                
            plt.xlabel('time (s)')
            plt.ylabel('D (cm2/s)')
            plt.title('Estimate of D with fit method')
            
            plt.savefig(outputfolder+'disp_gauss.png') 
            plt.close()
            
            Destfit=spanDest[-1]
            outputfile.write('\nEstimation of D from fit : %.2e'%Destfit)
            
            
            outputfile.write('\n\nEstimation of R from FWHM : %.2e'%(DestFWHM/D-1))
            outputfile.write('\nEstimation of R from fit : %.2e'%(Destfit/D-1))
            
            
            plt.figure()
            color='rgbkcmy'
            icol=0
            for i in range(0,len(spanDest)) :
                plt.plot(x,concentration[iperiodic[i+1]],'.',color=color[icol%len(color)])
                plt.plot(x,gaussian(np.array(x),t[iperiodic[i+1]]+sigma**2/(2*spanDest[i]),sigma,spanDest[i],L,xi),color=color[icol%len(color)],alpha=0.5)
                icol+=1
            plt.xlabel('x (cm)')
            plt.ylabel('concentration')
            plt.title('Fit of simulation results')
            #plt.xlim([L-50e-4,L])
            
            plt.savefig(outputfolder+'concentration.png') 
            plt.close()
            
            outputfile.close()
            #Update database
            Database.append([job,stage,bandname,vessel, Rv0, Rpvs, L, cell_size, time_step, rho, mu, D, sigma, xi, f, umax,pmax, Pe, Re, Wo, Fo , A, beta, dPdx, T,nPeriod , t[iperiodic[-1]], span_FWHM[iperiodic[-1]],DestFWHM, Destfit,DestFWHM/D-1, Destfit/D-1, amp, thetaa[-1],tau])
            

# save the database

file_name='../output/disp_analysis/disp'+analysis+outputname
labelstring=''
for d in datalabel:
    labelstring+=d+','
    
labelstring=labelstring[0:-1]

f = open(file_name, "w")
f.write(labelstring+'\n')
f.close()

formatstring='%s'+',%s'*3+' ,%.6e'*31

f = open(file_name, "a")
for d in Database :
    #print(formatstring%tuple(d))
    f.write(formatstring%tuple(d)+'\n')
f.close()

print('Error during the loading of :')
print(errorfiles)




