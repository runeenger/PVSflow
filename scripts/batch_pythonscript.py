#! /usr/bin/env python3


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 14:21:29 2021

@author: Alexandra Vallet



This scritps allow to write slurm files and a batch file in order to launch the PVS simulations on the supercomputer.

The simulations can also be run locally in the docker container, by using the command lines written at the end of the slurm files.

"""

    
import numpy as np
import os
import time


from src.toolbox import U
from src.io.slurm_io import get_slurmtemplate,base_PVSBraincommandline, write_state_slurm
from src.io.statistics_reader import ReadRandomEffect

# Physical properties
# The molecular dimensions of some Dextran fractions
dextranradius={10:2.36,40:4.45,50:4.95,70:5.8,100:6.9,200:9.5,500:14.7,1000:19.9,2000:27, 0:5.8*10} # in nm
# as one correspond to the 70 kDA * 10 to take into account obstacles.    



def print_checkdt(tend,dt,Noutput,toutput):
    print('# n time steps : %i'%(tend/dt))
    print('# n output : %i'%(Noutput))
    print('# n step/output : %f'%(toutput/dt))
    print('# dt : %5e'%dt)

    if Noutput>500 :
                            print('# n time steps : %i'%(tend/dt))
                            print('# n output : %i'%(Noutput))
                            print('# n step/output : %f'%(toutput/dt))
                            print('# dt : %5e'%dt)

                            if Noutput>500 :
                                print('f',f)
                                print('tend',tend)
                                print('dt',dt)
                                print('N',Noutput)
                                print(Noutput)
                                print('Warning : too many outputs!\n')
                            if int(toutput/dt)== 0 :
                                print('dt',tend)
                                print('toutput',toutput)
                                print('Warning : bad choice of dt!\n')
                            if dt>toutput:
                                print('dt',tend)
                                print('toutput',toutput)
                                print('Warning : bad choice of dt!\n')
                            if dt>=tend:
                                print('dt',tend)
                                print('toutput',toutput)
                                print('Warning : bad choice of dt!\n')
                            if toutput>=tend:
                                print('tend',tend)
                                print('Warning : bad choice of t output!',toutput)
 
                        
                            print('\n') 
          
                            
    
if __name__ == '__main__':
    

    ##########################################################
    ### Parameters of the script : can be modified
    ##########################################################

    # We would like to cover those parameters range
    # Length of the PVS
    spanlpvs=[400e-4,600e-4]  
    # Molecular diffusion coefficient
    spandiffusion=[0.84e-7*2,0.068e-6]
    # Assumption for the cardiac peak velocity
    spanVcard=[10e-4,50e-4,100e-4]
    
    # factor to be used on the PVS thickness h0, or the area (normally should be 1)
    corrfactorh0=1
    corrfactorarea=1
    
    # Boolean to print some verification about the time parameters
    Checkdt= True


    ##########################################################
    ### script
    ##########################################################


    # Dispersion analysis
    # --------------------
    # for each state we need mean lumen radius, mean PVS thickness + amp period 
    # we compute each FB separatelly
    
    
    # serie name
    seriename='dispersionRandomWT10'
    
    # create a folder for the slurm files and the batch file

    if not os.path.exists('../output/'):
        os.makedirs('../output/')

    if not os.path.exists('../output/supercomputer/'):
        os.makedirs('../output/supercomputer/')

        
        
    seriedirroot='../output/supercomputer/'+seriename+'/'
    
    # awake data
    folder='../data/statistics/penetrating_arterioles_WT10/'
    vessel='PenetratingArterioles'
    mouse='WT10_'
    analysis='wake_'
    ftype='randomEffects.txt'
    
    
    awakedata={}
    awakedata['Rv']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'medianepisode_'+ftype)
    awakedata['h0']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'medianepisode_'+ftype)

    
    awakedata['ampcard']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'cardiac_amp_'+ftype)
    awakedata['ampresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'resp_amp_'+ftype)
    awakedata['ampLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'LF_amp_'+ftype)
    awakedata['ampVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'VLF_amp_'+ftype)

    awakedata['amp_endfootcard']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'cardiac_amp_'+ftype)
    awakedata['amp_endfootresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'resp_amp_'+ftype)
    awakedata['amp_endfootLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'LF_amp_'+ftype)
    awakedata['amp_endfootVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'VLF_amp_'+ftype)

    awakedata['amp_lumencard']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'cardiac_amp_'+ftype)
    awakedata['amp_lumenresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'resp_amp_'+ftype)
    awakedata['amp_lumenLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'LF_amp_'+ftype)
    awakedata['amp_lumenVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'VLF_amp_'+ftype)    


    awakedata['area']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'medianepisode_'+ftype)
    
    awakedata['amp_areacard']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'cardiac_amp_'+ftype)
    awakedata['amp_arearesp']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'resp_amp_'+ftype)
    awakedata['amp_areaLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'LF_amp_'+ftype)
    awakedata['amp_areaVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'VLF_amp_'+ftype)

    awakedata['periodcard']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'cardiac_period_'+ftype)
    awakedata['periodresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'resp_period_'+ftype)
    awakedata['periodLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'LF_period_'+ftype)
    awakedata['periodVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'VLF_period_'+ftype)
    
    
    # sleep data
    analysis='sleep_'

    
    sleepdata={}
    sleepdata['Rv']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'medianepisode_'+ftype)
    sleepdata['h0']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'medianepisode_'+ftype)
    
    sleepdata['ampcard']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'cardiac_amp_'+ftype)
    sleepdata['ampresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'resp_amp_'+ftype)
    sleepdata['ampLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'LF_amp_'+ftype)
    sleepdata['ampVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'VLF_amp_'+ftype)

    sleepdata['amp_endfootcard']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'cardiac_amp_'+ftype)
    sleepdata['amp_endfootresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'resp_amp_'+ftype)
    sleepdata['amp_endfootLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'LF_amp_'+ftype)
    sleepdata['amp_endfootVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'VLF_amp_'+ftype)

    sleepdata['amp_lumencard']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'cardiac_amp_'+ftype)
    sleepdata['amp_lumenresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'resp_amp_'+ftype)
    sleepdata['amp_lumenLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'LF_amp_'+ftype)
    sleepdata['amp_lumenVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'VLF_amp_'+ftype)   
    
    sleepdata['amp_areacard']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'cardiac_amp_'+ftype)
    sleepdata['amp_arearesp']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'resp_amp_'+ftype)
    sleepdata['amp_areaLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'LF_amp_'+ftype)
    sleepdata['amp_areaVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'VLF_amp_'+ftype)

    sleepdata['area']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'medianepisode_'+ftype)

    sleepdata['periodcard']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'cardiac_period_'+ftype)
    sleepdata['periodresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'resp_period_'+ftype)
    sleepdata['periodLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'LF_period_'+ftype)
    sleepdata['periodVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'VLF_period_'+ftype)
    
    # get the names of the stages
    # for now just sleeping
    spanstages=sleepdata['Rv'] 
    
                        
    # get the number of vessels
    numbervessels=len(sleepdata['Rv']['baseline'])
    
    def intersect(*d):
        result = set(d[0]).intersection(*d[1:])
        return result

    
    amplituderandom=[]

    for lpvs in spanlpvs:
        for d in spandiffusion :
            
            nl=int(lpvs/1e-4)
            nr=8
            sigma=2e-4
            
            # serie name
            serie=seriename+'-d%.0e'%d+'-l%.0e'%lpvs
    
            # create a folder for the slurm files and the batch file
            if not os.path.exists(seriedirroot+serie):
                os.makedirs(seriedirroot+serie)
            
            seriedir=seriedirroot+serie+'/'
            
            slurmfiles=[]
            for stage in spanstages:
                
                    
                ### We generale area deformation for the cardiac FB in order to impose the velocity
                    
                # process the cardiac time scale
                # We assume the value of the velocity due to cardiac pulsation
                for vcard in spanVcard:
                    #Umax = a w L . we fix L and take w from measurement. 
                    #It gives us the amplitude of deformation of area to impose.
                    
                    #get the vessel IDs
                    vesselid_Rv = sleepdata['Rv'][stage].keys()
                    vesselid_h0 = sleepdata['h0'][stage].keys()

                    vesselid_osc = sleepdata['amp_areacard'][stage].keys()


                    vesselIDs=intersect(vesselid_Rv,vesselid_h0,vesselid_osc) 
                   

                
                    for vesselID in vesselIDs :
                
                    
                        jobname='disp'+'-d%.0e'%d+'-l%.0e'%lpvs+'-'+stage+'-'+'card'+'-v%.0e'%vcard+'-id'+vesselID
                        print(jobname)
                        
                        #get properties
                        rv=sleepdata['Rv'][stage][vesselID]*1e-4 #cm
                        h0=sleepdata['h0'][stage][vesselID]*1e-4*corrfactorh0 #cm
                        rpvs=rv+h0
                        fi=1/sleepdata['period'+'card'][stage][vesselID]
                        
                        
                        wi=2*np.pi*fi
                        
                        ampArea=vcard/wi/lpvs
                
                       
                        ai=ampArea
                        
                        print('rv:',rv)
                        print('h0:',h0)
                        print('f:',fi)
                        print('a:',ai)
                        

                        #estimate best time parameters
                        
                        #the end of the simulation should at least cover 5 periods (to get a linear fit)
                        # and the time to diffuse in the transveral direction
                        #tend=max(5/fi,2*h0**2/d)
                        

                        # Alternatively we decide to impose the final time
                        tend=40
                
                        # we would like 4 output per period
                        toutput=1/fi*4
                        
                        # we set dt as the toutput
                        dt=toutput
                        
                        # we would like at least 1000 time steps
                        while dt>(tend/1000) :
                            dt/=2
                
                        # we would like at max dt=5e-3
                        while dt>(5e-3) :
                            dt/=2
                        
                        #and at least 100 output over the whole simulation
                        while tend/toutput <100 :
                            toutput/=2
                
                        Noutput=tend/toutput
                
                        #estimate of the max velocity
                        Umax=U(0, 0, ai, fi, lpvs, Rv0=rv, h0=h0)
                        #estimate of the Peclet number
                        Pe=h0*Umax/d/2
                        # constraining velocity for the CFL condition
                        uconstrain=max(d/h0,Umax)

                        dx=h0/nr
                        
                        CFL=(uconstrain/(dx/dt))
                        while CFL>0.7 :
                            dt/=2
                            CFL=(uconstrain/(dx/dt))
                            
                        print('# CFL :' ,CFL)
                        print('# Pe : %e'%Pe)
                        
                        if Checkdt :
                            print_checkdt(tend,dt,Noutput,toutput)
 
                        
                        slurmfile=jobname+'.slurm'
                        
                        # #write the slurm file
                        # #gaussian analysis
                        write_state_slurm(jobname, [fi],[ai],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='gaussian', c0valuePVS=1, c0valueSAS=0, sigma=sigma, sasbc='scenarioA', refineleft=False)                              
                        
                        # #update the list of slurm files to be launched in batch
                        slurmfiles.append(slurmfile)
                        
                        
                # then we treat the measure values
            
                # each FB separately
                for FB in ['LF','VLF'] :

                    #get the vessel ID present for all freq
                    vesselid_Rv = sleepdata['Rv'][stage].keys()
                    vesselid_h0 = sleepdata['h0'][stage].keys()

                    vesselid_osc = sleepdata['amp_area'+FB][stage].keys()


                    vesselIDs=intersect(vesselid_Rv,vesselid_h0,vesselid_osc)  
                   
                
                    for vesselID in vesselIDs :
                    
                        jobname='disp'+'-d%.0e'%d+'-l%.0e'%lpvs+'-'+stage+'-'+FB+'-id'+vesselID
                        print(jobname)
                        
                        #get properties
                        rv=sleepdata['Rv'][stage][vesselID]*1e-4 #cm
                        h0=sleepdata['h0'][stage][vesselID]*1e-4*corrfactorh0 #cm
                        rpvs=rv+h0
                        fi=1/sleepdata['period'+FB][stage][vesselID]
                        
                        meanarea=sleepdata['area'][stage][vesselID]*corrfactorarea
                    
                    
                        amp_pvs=sleepdata['amp'+FB][stage][vesselID]*1e-4/2*corrfactorh0
                        amp_endfoot=sleepdata['amp_endfoot'+FB][stage][vesselID]*1e-4/2
                        amp_lumen=sleepdata['amp_lumen'+FB][stage][vesselID]*1e-4/2
                        
                        fi=1/sleepdata['period'+FB][stage][vesselID]
                        
                        
                        
                        # Compute area change from PVS thickness oscillations
                        # A0=np.pi*(rpvs)**2- np.pi*(rv)**2
                        # # compute the change of area
                        # #This doesnt work because we can get negative amplitude. 
                        # #This is because we should not use estimates of radius to compute area
                        # # We would need a model of the area ampliture
                        
                        # #Amin=np.pi*(rpvs+amp_endfoot)**2- np.pi*(rv+amp_lumen)**2
                        # #Amax=np.pi*(rpvs-amp_endfoot)**2- np.pi*(rv-amp_lumen)**2
                        
                        # # So what we do is to fix Rpvs and compute area change using the change of thickness
                        # Amin=np.pi*(rpvs)**2- np.pi*(rv+amp_pvs)**2
                        # Amax=np.pi*(rpvs)**2- np.pi*(rv-amp_pvs)**2        
                        
                        #if Amin>Amax : 
                        #    print('Problem in estimating the area !')
                        #    stop()  
                        # ai=(Amax-Amin)/A0/2
                        
                        
                        
                        # Compute area change directly from area oscilations
                        ai=sleepdata['amp_area'+FB][stage][vesselID]/meanarea*corrfactorarea
                        ai/=2
                        
                        amplituderandom.append(ai)
                        
                          
                        
                        print('stage',stage)
                        print('vessel',vesselID)
                        print('FB',FB)
                        print('rpvs:',rpvs)
                        print('rv:',rv)
                        
                        print('arpvs:',amp_endfoot)
                        print('arv:',amp_lumen)
                        
                        print('h0:',h0)
                        print('f:',fi)
                        print('a:',ai)
    
                                          
    
                        #estimate best time parameters
                        
                        #the end of the simulation should at least cover 3 periods (to get a linear fit)
                        # and the time to diffuse in the transveral direction
                        #tend=max(5/fi,2*h0**2/d)
                        
                        tend=40

                        # we would like 4 output per period
                        toutput=1/fi*4
                        
                        # we set dt as the toutput
                        dt=toutput
                        
                        # we would like at least 1000 time steps
                        while dt>(tend/1000) :
                            dt/=2

                        # we would like at max dt=5e-3
                        while dt>(5e-3) :
                            dt/=2
                        
                        #and at least 100 output over the whole simulation
                        while tend/toutput <100 :
                            toutput/=2
            
                        Noutput=tend/toutput
        
                        #estiamate of the max velocity
                        Umax=U(0, 0, ai, fi, lpvs, Rv0=rv, h0=h0)
                        #estimate of the Peclet number
                        Pe=h0*Umax/d/2
                        # constraining velocity for the CFL condition
                        uconstrain=max(d/h0,Umax)
                        # number of cells in transversal direction
                        nr=8
                        dx=h0/nr
                        
                        CFL=(uconstrain/(dx/dt))
                        while CFL>0.7 :
                            dt/=2
                            CFL=(uconstrain/(dx/dt))
                            
                        print('# CFL :' ,CFL)
                        print('# Pe : %e'%Pe)

                        if  Checkdt:
                            print_checkdt(tend,dt,Noutput,toutput)
                        
                        slurmfile=jobname+'.slurm'

                        
                        #write the slurm file
                        # gaussian analysis
                        write_state_slurm(jobname, [fi],[ai],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='gaussian', c0valuePVS=1, c0valueSAS=0, sigma=sigma, sasbc='scenarioA', refineleft=False)                              

                        
                        # intake
                        #write_state_slurm(jobname, [fi],[ai],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='uniform', c0valuePVS=0, c0valueSAS=1, sigma=sigma, sasbc='scenarioE', refineleft=False)                              

                        #update the list of slurm files to be launched in batch
                        slurmfiles.append(slurmfile)
                        
                    # superposition of LF anf VLF
                    
                    
                    #get the vessel ID present for all freq
                    vesselid_Rv = sleepdata['Rv'][stage].keys()
                    vesselid_h0 = sleepdata['h0'][stage].keys()

                    vesselid_oscLF = sleepdata['amp_area'+'LF'][stage].keys()
                    vesselid_oscVLF = sleepdata['amp_area'+'VLF'][stage].keys()


                    vesselIDs=intersect(vesselid_Rv,vesselid_h0,vesselid_oscLF,vesselid_oscVLF)  
                   
                
                    for vesselID in vesselIDs :
                    
                        jobname='disp'+'-d%.0e'%d+'-l%.0e'%lpvs+'-'+stage+'-'+'LFVLF'+'-id'+vesselID
                        print(jobname)
                        
                        #get properties
                        rv=sleepdata['Rv'][stage][vesselID]*1e-4 #cm
                        h0=sleepdata['h0'][stage][vesselID]*1e-4*corrfactorh0 #cm
                        rpvs=rv+h0
                        fi=1/sleepdata['period'+FB][stage][vesselID]
                        
                        meanarea=sleepdata['area'][stage][vesselID]*corrfactorarea
                    
                    
                       
                        fVLF=1/sleepdata['period'+'VLF'][stage][vesselID]
                        fLF=1/sleepdata['period'+'LF'][stage][vesselID]
                        
                        
                        # Compute area change directly from area oscilations
                        aVLF=sleepdata['amp_area'+'VLF'][stage][vesselID]/meanarea/2*corrfactorarea
                        aLF=sleepdata['amp_area'+'LF'][stage][vesselID]/meanarea/2*corrfactorarea
                        
                        
                          
                        
                        print('stage',stage)
                        print('vessel',vesselID)
                        print('FB',FB)
                        print('rpvs:',rpvs)
                        print('rv:',rv)
                        
                        print('arpvs:',amp_endfoot)
                        print('arv:',amp_lumen)
                        
                        print('h0:',h0)
                        print('f:',[fVLF,fLF])
                        print('a:',[aVLF,aLF])
    
                                          
    
                        #estimate best time parameters
                        
                        #the end of the simulation should at least cover 3 periods (to get a linear fit)
                        # and the time to diffuse in the transveral direction
                        #tend=max(5/fi,2*h0**2/d)
                        
                        tend=40

                        # we would like 4 output per period
                        toutput=1/fLF*4
                        
                        # we set dt as the toutput
                        dt=toutput
                        
                        # we would like at least 1000 time steps
                        while dt>(tend/1000) :
                            dt/=2

                        # we would like at max dt=5e-3
                        while dt>(5e-3) :
                            dt/=2
                        
                        #and at least 100 output over the whole simulation
                        while tend/toutput <100 :
                            toutput/=2
            
                        Noutput=tend/toutput
        
                        #estiamate of the max velocity
                        Umax=U(0, 0, aLF+aVLF, fLF, lpvs, Rv0=rv, h0=h0)
                        #estimate of the Peclet number
                        Pe=h0*Umax/d/2
                        # constraining velocity for the CFL condition
                        uconstrain=max(d/h0,Umax)
                        # number of cells in transversal direction
                        nr=8
                        dx=h0/nr
                        
                        CFL=(uconstrain/(dx/dt))
                        while CFL>0.7 :
                            dt/=2
                            CFL=(uconstrain/(dx/dt))
                            
                        print('# CFL :' ,CFL)
                        print('# Pe : %e'%Pe)
                        
                        if  Checkdt:
                            print_checkdt(tend,dt,Noutput,toutput)

                        slurmfile=jobname+'.slurm'

                        
                        #write the slurm file
                        # gaussian analysis
                        write_state_slurm(jobname, [fi],[ai],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='gaussian', c0valuePVS=1, c0valueSAS=0, sigma=sigma, sasbc='scenarioA', refineleft=False)                              

                        #update the list of slurm files to be launched in batch
                        slurmfiles.append(slurmfile)
                        
                    

            #write the batch file
            with open(seriedir+'batch','w') as f :
                i=1
                for slurmfile in slurmfiles :
                    f.write('sbatch '+slurmfile+' &\n')
                    if not(i%100):
                        f.write('sleep 1 \n')
                    i+=1
            f.close()

   # Dispersion analaysis for veins
    # --------------------
    # for each state we need mean lumen radius, mean PVS thickness + amp period 
    # we compute each FB separatelly
    
    
    # serie name
    seriename='dispersionRandomveinsWT10'
    
    # create a folder for the slurm files and the batch file

    if not os.path.exists('../output/'):
        os.makedirs('../output/')

    if not os.path.exists('../output/supercomputer/'):
        os.makedirs('../output/supercomputer/')


        
        
    seriedirroot='../output/supercomputer/'+seriename+'/'
    
    # sleep data
    folder='../data/statistics/veins_WT10/'
    vessel='Veins'
    mouse='WT10_'
    ftype='randomEffects.txt'
    analysis='sleep_'

    
    sleepdata={}
    sleepdata['Rv']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'medianepisode_'+ftype)
    sleepdata['h0']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'medianepisode_'+ftype)
    
    sleepdata['ampcard']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'cardiac_amp_'+ftype)
    #sleepdata['ampresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'resp_amp_'+ftype)
    sleepdata['ampLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'LF_amp_'+ftype)
    sleepdata['ampVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'VLF_amp_'+ftype)

    sleepdata['amp_endfootcard']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'cardiac_amp_'+ftype)
    #sleepdata['amp_endfootresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'resp_amp_'+ftype)
    sleepdata['amp_endfootLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'LF_amp_'+ftype)
    sleepdata['amp_endfootVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'VLF_amp_'+ftype)

    sleepdata['amp_lumencard']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'cardiac_amp_'+ftype)
    #sleepdata['amp_lumenresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'resp_amp_'+ftype)
    sleepdata['amp_lumenLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'LF_amp_'+ftype)
    sleepdata['amp_lumenVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'VLF_amp_'+ftype)   
    
    sleepdata['amp_areacard']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'cardiac_amp_'+ftype)
    #sleepdata['amp_arearesp']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'resp_amp_'+ftype)
    sleepdata['amp_areaLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'LF_amp_'+ftype)
    sleepdata['amp_areaVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'VLF_amp_'+ftype)

    sleepdata['area']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'medianepisode_'+ftype)

    sleepdata['periodcard']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'cardiac_period_'+ftype)
    #sleepdata['periodresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'resp_period_'+ftype)
    sleepdata['periodLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'LF_period_'+ftype)
    sleepdata['periodVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'VLF_period_'+ftype)
    
    # get the names of the stages
    # for now just sleeping
    spanstages=sleepdata['Rv'] 
    
                        
    # get the number of vessels
    numbervessels=len(sleepdata['Rv']['baseline'])
    
    def intersect(*d):
        result = set(d[0]).intersection(*d[1:])
        return result

    
    amplituderandom=[]

    for lpvs in spanlpvs:
        for d in spandiffusion :
            
            nl=int(lpvs/1e-4)
            nr=8
            sigma=2e-4
            
            # serie name
            serie=seriename+'-d%.0e'%d+'-l%.0e'%lpvs
    
            # create a folder for the slurm files and the batch file
            if not os.path.exists(seriedirroot+serie):
                os.makedirs(seriedirroot+serie)
            
            seriedir=seriedirroot+serie+'/'
    
            
            slurmfiles=[]
            for stage in spanstages:
                
                    
                ### We generale area deformation for the cardiac FB in order to impose the velocity
                    
                # process the cardiac time scale
                # We assume the value of the velocity due to cardiac pulsation
                for vcard in spanVcard:
                    #Umax = a w L . we fix L and take w from measurement. 
                    #It gives us the amplitude of deformation of area to impose.
                    
                    #get the vessel IDs
                    vesselid_Rv = sleepdata['Rv'][stage].keys()
                    vesselid_h0 = sleepdata['h0'][stage].keys()

                    vesselid_osc = sleepdata['amp_areacard'][stage].keys()


                    vesselIDs=intersect(vesselid_Rv,vesselid_h0,vesselid_osc) 
                   

                
                    for vesselID in vesselIDs :
                
                    
                        jobname='disp'+'-d%.0e'%d+'-l%.0e'%lpvs+'-'+stage+'-'+'card'+'-v%.0e'%vcard+'-id'+vesselID
                        print(jobname)
                        
                        #get properties
                        rv=sleepdata['Rv'][stage][vesselID]*1e-4 #cm
                        h0=sleepdata['h0'][stage][vesselID]*1e-4*corrfactorh0 #cm
                        rpvs=rv+h0
                        fi=1/sleepdata['period'+'card'][stage][vesselID]
                        
                        
                        wi=2*np.pi*fi
                        
                        ampArea=vcard/wi/lpvs
                
                       
                        ai=ampArea
                        
                        print('rv:',rv)
                        print('h0:',h0)
                        print('f:',fi)
                        print('a:',ai)
                        

                        #estimate best time parameters
                        
                        #the end of the simulation should at least cover 5 periods (to get a linear fit)
                        # and the time to diffuse in the transveral direction
                        #tend=max(5/fi,2*h0**2/d)
                        

                        # Alternatively we decide to impose the final time
                        tend=40
                
                        # we would like 4 output per period
                        toutput=1/fi*4
                        
                        # we set dt as the toutput
                        dt=toutput
                        
                        # we would like at least 1000 time steps
                        while dt>(tend/1000) :
                            dt/=2
                
                        # we would like at max dt=5e-3
                        while dt>(5e-3) :
                            dt/=2
                        
                        #and at least 100 output over the whole simulation
                        while tend/toutput <100 :
                            toutput/=2
                
                        Noutput=tend/toutput
                
                        #estimate of the max velocity
                        Umax=U(0, 0, ai, fi, lpvs, Rv0=rv, h0=h0)
                        #estimate of the Peclet number
                        Pe=h0*Umax/d/2
                        # constraining velocity for the CFL condition
                        uconstrain=max(d/h0,Umax)

                        dx=h0/nr
                        
                        CFL=(uconstrain/(dx/dt))
                        while CFL>0.7 :
                            dt/=2
                            CFL=(uconstrain/(dx/dt))
                            
                        print('# CFL :' ,CFL)
                        print('# Pe : %e'%Pe)
                        
                        if Checkdt :
                            print_checkdt(tend,dt,Noutput,toutput)
 
                        
                        slurmfile=jobname+'.slurm'
                        
                        # #write the slurm file
                        # #gaussian analysis
                        write_state_slurm(jobname, [fi],[ai],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='gaussian', c0valuePVS=1, c0valueSAS=0, sigma=sigma, sasbc='scenarioA', refineleft=False)                              
                        
                        # #update the list of slurm files to be launched in batch
                        slurmfiles.append(slurmfile)
                        
                        
                # then we treat the measure values
            
                # each FB separately
                for FB in ['LF','VLF'] :

                    #get the vessel ID present for all freq
                    vesselid_Rv = sleepdata['Rv'][stage].keys()
                    vesselid_h0 = sleepdata['h0'][stage].keys()

                    vesselid_osc = sleepdata['amp_area'+FB][stage].keys()


                    vesselIDs=intersect(vesselid_Rv,vesselid_h0,vesselid_osc)  
                   
                
                    for vesselID in vesselIDs :
                    
                        jobname='disp'+'-d%.0e'%d+'-l%.0e'%lpvs+'-'+stage+'-'+FB+'-id'+vesselID
                        print(jobname)
                        
                        #get properties
                        rv=sleepdata['Rv'][stage][vesselID]*1e-4 #cm
                        h0=sleepdata['h0'][stage][vesselID]*1e-4*corrfactorh0 #cm
                        rpvs=rv+h0
                        fi=1/sleepdata['period'+FB][stage][vesselID]
                        
                        meanarea=sleepdata['area'][stage][vesselID]*corrfactorarea
                    
                    
                        amp_pvs=sleepdata['amp'+FB][stage][vesselID]*1e-4/2*corrfactorh0
                        amp_endfoot=sleepdata['amp_endfoot'+FB][stage][vesselID]*1e-4/2
                        amp_lumen=sleepdata['amp_lumen'+FB][stage][vesselID]*1e-4/2
                        
                        fi=1/sleepdata['period'+FB][stage][vesselID]
                        
                        
                        
                        # Compute area change from PVS thickness oscillations
                        # A0=np.pi*(rpvs)**2- np.pi*(rv)**2
                        # # compute the change of area
                        # #This doesnt work because we can get negative amplitude. 
                        # #This is because we should not use estimates of radius to compute area
                        # # We would need a model of the area ampliture
                        
                        # #Amin=np.pi*(rpvs+amp_endfoot)**2- np.pi*(rv+amp_lumen)**2
                        # #Amax=np.pi*(rpvs-amp_endfoot)**2- np.pi*(rv-amp_lumen)**2
                        
                        # # So what we do is to fix Rpvs and compute area change using the change of thickness
                        # Amin=np.pi*(rpvs)**2- np.pi*(rv+amp_pvs)**2
                        # Amax=np.pi*(rpvs)**2- np.pi*(rv-amp_pvs)**2        
                        
                        #if Amin>Amax : 
                        #    print('Problem in estimating the area !')
                        #    stop()  
                        # ai=(Amax-Amin)/A0/2
                        
                        
                        
                        # Compute area change directly from area oscilations
                        ai=sleepdata['amp_area'+FB][stage][vesselID]/meanarea*corrfactorarea
                        ai/=2
                        
                        amplituderandom.append(ai)
                        
                          
                        
                        print('stage',stage)
                        print('vessel',vesselID)
                        print('FB',FB)
                        print('rpvs:',rpvs)
                        print('rv:',rv)
                        
                        print('arpvs:',amp_endfoot)
                        print('arv:',amp_lumen)
                        
                        print('h0:',h0)
                        print('f:',fi)
                        print('a:',ai)
    
                                          
    
                        #estimate best time parameters
                        
                        #the end of the simulation should at least cover 3 periods (to get a linear fit)
                        # and the time to diffuse in the transveral direction
                        #tend=max(5/fi,2*h0**2/d)
                        
                        tend=40

                        # we would like 4 output per period
                        toutput=1/fi*4
                        
                        # we set dt as the toutput
                        dt=toutput
                        
                        # we would like at least 1000 time steps
                        while dt>(tend/1000) :
                            dt/=2

                        # we would like at max dt=5e-3
                        while dt>(5e-3) :
                            dt/=2
                        
                        #and at least 100 output over the whole simulation
                        while tend/toutput <100 :
                            toutput/=2
            
                        Noutput=tend/toutput
        
                        #estiamate of the max velocity
                        Umax=U(0, 0, ai, fi, lpvs, Rv0=rv, h0=h0)
                        #estimate of the Peclet number
                        Pe=h0*Umax/d/2
                        # constraining velocity for the CFL condition
                        uconstrain=max(d/h0,Umax)
                        # number of cells in transversal direction
                        nr=8
                        dx=h0/nr
                        
                        CFL=(uconstrain/(dx/dt))
                        while CFL>0.7 :
                            dt/=2
                            CFL=(uconstrain/(dx/dt))
                            
                        print('# CFL :' ,CFL)
                        print('# Pe : %e'%Pe)

                        if  Checkdt:
                            print_checkdt(tend,dt,Noutput,toutput)
                        
                        slurmfile=jobname+'.slurm'

                        
                        #write the slurm file
                        # gaussian analysis
                        write_state_slurm(jobname, [fi],[ai],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='gaussian', c0valuePVS=1, c0valueSAS=0, sigma=sigma, sasbc='scenarioA', refineleft=False)                              

                        
                        # intake
                        #write_state_slurm(jobname, [fi],[ai],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='uniform', c0valuePVS=0, c0valueSAS=1, sigma=sigma, sasbc='scenarioE', refineleft=False)                              

                        #update the list of slurm files to be launched in batch
                        slurmfiles.append(slurmfile)
                        
                    # superposition of LF anf VLF
                    
                    
                    #get the vessel ID present for all freq
                    vesselid_Rv = sleepdata['Rv'][stage].keys()
                    vesselid_h0 = sleepdata['h0'][stage].keys()

                    vesselid_oscLF = sleepdata['amp_area'+'LF'][stage].keys()
                    vesselid_oscVLF = sleepdata['amp_area'+'VLF'][stage].keys()


                    vesselIDs=intersect(vesselid_Rv,vesselid_h0,vesselid_oscLF,vesselid_oscVLF)  
                   
                
                    for vesselID in vesselIDs :
                    
                        jobname='disp'+'-d%.0e'%d+'-l%.0e'%lpvs+'-'+stage+'-'+'LFVLF'+'-id'+vesselID
                        print(jobname)
                        
                        #get properties
                        rv=sleepdata['Rv'][stage][vesselID]*1e-4 #cm
                        h0=sleepdata['h0'][stage][vesselID]*1e-4*corrfactorh0 #cm
                        rpvs=rv+h0
                        fi=1/sleepdata['period'+FB][stage][vesselID]
                        
                        meanarea=sleepdata['area'][stage][vesselID]*corrfactorarea
                    
                    
                       
                        fVLF=1/sleepdata['period'+'VLF'][stage][vesselID]
                        fLF=1/sleepdata['period'+'LF'][stage][vesselID]
                        
                        
                        # Compute area change directly from area oscilations
                        aVLF=sleepdata['amp_area'+'VLF'][stage][vesselID]/meanarea/2*corrfactorarea
                        aLF=sleepdata['amp_area'+'LF'][stage][vesselID]/meanarea/2*corrfactorarea
                        
                        
                          
                        
                        print('stage',stage)
                        print('vessel',vesselID)
                        print('FB',FB)
                        print('rpvs:',rpvs)
                        print('rv:',rv)
                        
                        print('arpvs:',amp_endfoot)
                        print('arv:',amp_lumen)
                        
                        print('h0:',h0)
                        print('f:',[fVLF,fLF])
                        print('a:',[aVLF,aLF])
    
                                          
    
                        #estimate best time parameters
                        
                        #the end of the simulation should at least cover 3 periods (to get a linear fit)
                        # and the time to diffuse in the transveral direction
                        #tend=max(5/fi,2*h0**2/d)
                        
                        tend=40

                        # we would like 4 output per period
                        toutput=1/fLF*4
                        
                        # we set dt as the toutput
                        dt=toutput
                        
                        # we would like at least 1000 time steps
                        while dt>(tend/1000) :
                            dt/=2

                        # we would like at max dt=5e-3
                        while dt>(5e-3) :
                            dt/=2
                        
                        #and at least 100 output over the whole simulation
                        while tend/toutput <100 :
                            toutput/=2
            
                        Noutput=tend/toutput
        
                        #estiamate of the max velocity
                        Umax=U(0, 0, aLF+aVLF, fLF, lpvs, Rv0=rv, h0=h0)
                        #estimate of the Peclet number
                        Pe=h0*Umax/d/2
                        # constraining velocity for the CFL condition
                        uconstrain=max(d/h0,Umax)
                        # number of cells in transversal direction
                        nr=8
                        dx=h0/nr
                        
                        CFL=(uconstrain/(dx/dt))
                        while CFL>0.7 :
                            dt/=2
                            CFL=(uconstrain/(dx/dt))
                            
                        print('# CFL :' ,CFL)
                        print('# Pe : %e'%Pe)
                        
                        if  Checkdt:
                            print_checkdt(tend,dt,Noutput,toutput)

                        slurmfile=jobname+'.slurm'

                        
                        #write the slurm file
                        # gaussian analysis
                        write_state_slurm(jobname, [fi],[ai],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='gaussian', c0valuePVS=1, c0valueSAS=0, sigma=sigma, sasbc='scenarioA', refineleft=False)                              

                        #update the list of slurm files to be launched in batch
                        slurmfiles.append(slurmfile)
                        
                    

            #write the batch file
            with open(seriedir+'batch','w') as f :
                i=1
                for slurmfile in slurmfiles :
                    f.write('sbatch '+slurmfile+' &\n')
                    if not(i%100):
                        f.write('sleep 1 \n')
                    i+=1
            f.close()
    

    # Dispersion analysis with impermeable SMC
    # --------------------
    # for each state we need mean lumen radius, mean PVS thickness + amp period 
    # we compute each FB separatelly
    
    
    # serie name
    seriename='dispersionSMC50WT10'
    
    # create a folder for the slurm files and the batch file

    if not os.path.exists('../output/'):
        os.makedirs('../output/')

    if not os.path.exists('../output/supercomputer/'):
        os.makedirs('../output/supercomputer/')

        
        
    seriedirroot='../output/supercomputer/'+seriename+'/'
    
    # awake data
    folder='../data/statistics/penetrating_arterioles_WT10/'
    vessel='PenetratingArterioles'
    mouse='WT10_'
    analysis='wake_'
    ftype='randomEffects.txt'
    
    
    awakedata={}
    awakedata['Rv']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'medianepisode_'+ftype)
    awakedata['h0']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'medianepisode_'+ftype)

    
    awakedata['ampcard']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'cardiac_amp_'+ftype)
    awakedata['ampresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'resp_amp_'+ftype)
    awakedata['ampLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'LF_amp_'+ftype)
    awakedata['ampVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'VLF_amp_'+ftype)

    awakedata['amp_endfootcard']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'cardiac_amp_'+ftype)
    awakedata['amp_endfootresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'resp_amp_'+ftype)
    awakedata['amp_endfootLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'LF_amp_'+ftype)
    awakedata['amp_endfootVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'VLF_amp_'+ftype)

    awakedata['amp_lumencard']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'cardiac_amp_'+ftype)
    awakedata['amp_lumenresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'resp_amp_'+ftype)
    awakedata['amp_lumenLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'LF_amp_'+ftype)
    awakedata['amp_lumenVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'VLF_amp_'+ftype)    


    awakedata['area']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'medianepisode_'+ftype)
    
    awakedata['amp_areacard']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'cardiac_amp_'+ftype)
    awakedata['amp_arearesp']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'resp_amp_'+ftype)
    awakedata['amp_areaLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'LF_amp_'+ftype)
    awakedata['amp_areaVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'VLF_amp_'+ftype)

    awakedata['periodcard']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'cardiac_period_'+ftype)
    awakedata['periodresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'resp_period_'+ftype)
    awakedata['periodLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'LF_period_'+ftype)
    awakedata['periodVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'VLF_period_'+ftype)
    
    
    # sleep data
    analysis='sleep_'

    
    sleepdata={}
    sleepdata['Rv']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'medianepisode_'+ftype)
    sleepdata['h0']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'medianepisode_'+ftype)
    
    sleepdata['ampcard']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'cardiac_amp_'+ftype)
    sleepdata['ampresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'resp_amp_'+ftype)
    sleepdata['ampLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'LF_amp_'+ftype)
    sleepdata['ampVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'VLF_amp_'+ftype)

    sleepdata['amp_endfootcard']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'cardiac_amp_'+ftype)
    sleepdata['amp_endfootresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'resp_amp_'+ftype)
    sleepdata['amp_endfootLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'LF_amp_'+ftype)
    sleepdata['amp_endfootVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'endfoot_'+'VLF_amp_'+ftype)

    sleepdata['amp_lumencard']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'cardiac_amp_'+ftype)
    sleepdata['amp_lumenresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'resp_amp_'+ftype)
    sleepdata['amp_lumenLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'LF_amp_'+ftype)
    sleepdata['amp_lumenVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'VLF_amp_'+ftype)   
    
    sleepdata['amp_areacard']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'cardiac_amp_'+ftype)
    sleepdata['amp_arearesp']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'resp_amp_'+ftype)
    sleepdata['amp_areaLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'LF_amp_'+ftype)
    sleepdata['amp_areaVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'VLF_amp_'+ftype)

    sleepdata['area']=ReadRandomEffect(folder+vessel+mouse+analysis+'Area_'+'medianepisode_'+ftype)

    sleepdata['periodcard']=ReadRandomEffect(folder+vessel+mouse+analysis+'lumen_'+'cardiac_period_'+ftype)
    sleepdata['periodresp']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'resp_period_'+ftype)
    sleepdata['periodLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'LF_period_'+ftype)
    sleepdata['periodVLF']=ReadRandomEffect(folder+vessel+mouse+analysis+'PVS_'+'VLF_period_'+ftype)
    
    # get the names of the stages
    # for now just sleeping
    spanstages=sleepdata['Rv'] 
    
                        
    # get the number of vessels
    numbervessels=len(sleepdata['Rv']['baseline'])
    
    def intersect(*d):
        result = set(d[0]).intersection(*d[1:])
        return result

    
    amplituderandom=[]

    for lpvs in spanlpvs:
        for d in spandiffusion :
            
            nl=int(lpvs/1e-4)
            nr=8
            sigma=2e-4
            
            # serie name
            serie=seriename+'-d%.0e'%d+'-l%.0e'%lpvs
    
            # create a folder for the slurm files and the batch file
            if not os.path.exists(seriedirroot+serie):
                os.makedirs(seriedirroot+serie)
            
            seriedir=seriedirroot+serie+'/'
            
            slurmfiles=[]
            for stage in spanstages:
                
                    
                ### We generale area deformation for the cardiac FB in order to impose the velocity
                    
                # process the cardiac time scale
                # We assume the value of the velocity due to cardiac pulsation
                for vcard in spanVcard:
                    #Umax = a w L . we fix L and take w from measurement. 
                    #It gives us the amplitude of deformation of area to impose.
                    
                    #get the vessel IDs
                    vesselid_Rv = sleepdata['Rv'][stage].keys()
                    vesselid_h0 = sleepdata['h0'][stage].keys()

                    vesselid_osc = sleepdata['amp_areacard'][stage].keys()


                    vesselIDs=intersect(vesselid_Rv,vesselid_h0,vesselid_osc) 
                   

                
                    for vesselID in vesselIDs :
                
                    
                        jobname='disp'+'-d%.0e'%d+'-l%.0e'%lpvs+'-'+stage+'-'+'card'+'-v%.0e'%vcard+'-id'+vesselID
                        print(jobname)
                        
                        #get properties
                        rv=sleepdata['Rv'][stage][vesselID]*1e-4 #cm
                        h0=sleepdata['h0'][stage][vesselID]*1e-4*corrfactorh0 #cm
                        rpvs=rv+h0
                        fi=1/sleepdata['period'+'card'][stage][vesselID]
                        
                        
                        wi=2*np.pi*fi
                        
                        ampArea=vcard/wi/lpvs*(1-aSMC)
                
                       
                        ai=ampArea
                        
                        print('rv:',rv)
                        print('h0:',h0)
                        print('f:',fi)
                        print('a:',ai)
                        

                        #estimate best time parameters
                        
                        #the end of the simulation should at least cover 5 periods (to get a linear fit)
                        # and the time to diffuse in the transveral direction
                        #tend=max(5/fi,2*h0**2/d)
                        

                        # Alternatively we decide to impose the final time
                        tend=40
                
                        # we would like 4 output per period
                        toutput=1/fi*4
                        
                        # we set dt as the toutput
                        dt=toutput
                        
                        # we would like at least 1000 time steps
                        while dt>(tend/1000) :
                            dt/=2
                
                        # we would like at max dt=5e-3
                        while dt>(5e-3) :
                            dt/=2
                        
                        #and at least 100 output over the whole simulation
                        while tend/toutput <100 :
                            toutput/=2
                
                        Noutput=tend/toutput
                
                        #estimate of the max velocity
                        Umax=U(0, 0, ai, fi, lpvs, Rv0=rv, h0=h0)
                        #estimate of the Peclet number
                        Pe=h0*Umax/d/2
                        # constraining velocity for the CFL condition
                        uconstrain=max(d/h0,Umax)

                        dx=h0/nr
                        
                        CFL=(uconstrain/(dx/dt))
                        while CFL>0.7 :
                            dt/=2
                            CFL=(uconstrain/(dx/dt))
                            
                        print('# CFL :' ,CFL)
                        print('# Pe : %e'%Pe)
                        
                        if Checkdt :
                            print_checkdt(tend,dt,Noutput,toutput)
 
                        
                        slurmfile=jobname+'.slurm'
                        
                        # #write the slurm file
                        # #gaussian analysis
                        write_state_slurm(jobname, [fi],[ai],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='gaussian', c0valuePVS=1, c0valueSAS=0, sigma=sigma, sasbc='scenarioA', refineleft=False, aSMC=aSMC)                              
                        
                        # #update the list of slurm files to be launched in batch
                        slurmfiles.append(slurmfile)
                        
                        
                # then we treat the measure values
            
                # each FB separately
                for FB in ['LF','VLF'] :

                    #get the vessel ID present for all freq
                    vesselid_Rv = sleepdata['Rv'][stage].keys()
                    vesselid_h0 = sleepdata['h0'][stage].keys()

                    vesselid_osc = sleepdata['amp_area'+FB][stage].keys()


                    vesselIDs=intersect(vesselid_Rv,vesselid_h0,vesselid_osc)  
                   
                
                    for vesselID in vesselIDs :
                    
                        jobname='disp'+'-d%.0e'%d+'-l%.0e'%lpvs+'-'+stage+'-'+FB+'-id'+vesselID
                        print(jobname)
                        
                        #get properties
                        rv=sleepdata['Rv'][stage][vesselID]*1e-4 #cm
                        h0=sleepdata['h0'][stage][vesselID]*1e-4*corrfactorh0 #cm
                        rpvs=rv+h0
                        fi=1/sleepdata['period'+FB][stage][vesselID]
                        
                        meanarea=sleepdata['area'][stage][vesselID]*corrfactorarea
                    
                    
                        amp_pvs=sleepdata['amp'+FB][stage][vesselID]*1e-4/2*corrfactorh0
                        amp_endfoot=sleepdata['amp_endfoot'+FB][stage][vesselID]*1e-4/2
                        amp_lumen=sleepdata['amp_lumen'+FB][stage][vesselID]*1e-4/2
                        
                        fi=1/sleepdata['period'+FB][stage][vesselID]
                        
                        
                        
                        # Compute area change from PVS thickness oscillations
                        # A0=np.pi*(rpvs)**2- np.pi*(rv)**2
                        # # compute the change of area
                        # #This doesnt work because we can get negative amplitude. 
                        # #This is because we should not use estimates of radius to compute area
                        # # We would need a model of the area ampliture
                        
                        # #Amin=np.pi*(rpvs+amp_endfoot)**2- np.pi*(rv+amp_lumen)**2
                        # #Amax=np.pi*(rpvs-amp_endfoot)**2- np.pi*(rv-amp_lumen)**2
                        
                        # # So what we do is to fix Rpvs and compute area change using the change of thickness
                        # Amin=np.pi*(rpvs)**2- np.pi*(rv+amp_pvs)**2
                        # Amax=np.pi*(rpvs)**2- np.pi*(rv-amp_pvs)**2        
                        
                        #if Amin>Amax : 
                        #    print('Problem in estimating the area !')
                        #    stop()  
                        # ai=(Amax-Amin)/A0/2
                        
                        
                        
                        # Compute area change directly from area oscilations
                        ai=sleepdata['amp_area'+FB][stage][vesselID]/meanarea*corrfactorarea
                        ai/=2
                        
                        amplituderandom.append(ai)
                        
                          
                        
                        print('stage',stage)
                        print('vessel',vesselID)
                        print('FB',FB)
                        print('rpvs:',rpvs)
                        print('rv:',rv)
                        
                        print('arpvs:',amp_endfoot)
                        print('arv:',amp_lumen)
                        
                        print('h0:',h0)
                        print('f:',fi)
                        print('a:',ai)
    
                                          
    
                        #estimate best time parameters
                        
                        #the end of the simulation should at least cover 3 periods (to get a linear fit)
                        # and the time to diffuse in the transveral direction
                        #tend=max(5/fi,2*h0**2/d)
                        
                        tend=40

                        # we would like 4 output per period
                        toutput=1/fi*4
                        
                        # we set dt as the toutput
                        dt=toutput
                        
                        # we would like at least 1000 time steps
                        while dt>(tend/1000) :
                            dt/=2

                        # we would like at max dt=5e-3
                        while dt>(5e-3) :
                            dt/=2
                        
                        #and at least 100 output over the whole simulation
                        while tend/toutput <100 :
                            toutput/=2
            
                        Noutput=tend/toutput
        
                        #estiamate of the max velocity
                        Umax=U(0, 0, ai, fi, lpvs, Rv0=rv, h0=h0)
                        #estimate of the Peclet number
                        Pe=h0*Umax/d/2
                        # constraining velocity for the CFL condition
                        uconstrain=max(d/h0,Umax)
                        # number of cells in transversal direction
                        nr=8
                        dx=h0/nr
                        
                        CFL=(uconstrain/(dx/dt))
                        while CFL>0.7 :
                            dt/=2
                            CFL=(uconstrain/(dx/dt))
                            
                        print('# CFL :' ,CFL)
                        print('# Pe : %e'%Pe)

                        if  Checkdt:
                            print_checkdt(tend,dt,Noutput,toutput)
                        
                        slurmfile=jobname+'.slurm'

                        
                        #write the slurm file
                        # gaussian analysis
                        write_state_slurm(jobname, [fi],[ai],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='gaussian', c0valuePVS=1, c0valueSAS=0, sigma=sigma, sasbc='scenarioA', refineleft=False, aSMC=aSMC)                              

                        
                        # intake
                        #write_state_slurm(jobname, [fi],[ai],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='uniform', c0valuePVS=0, c0valueSAS=1, sigma=sigma, sasbc='scenarioE', refineleft=False)                              

                        #update the list of slurm files to be launched in batch
                        slurmfiles.append(slurmfile)
                        
                    # superposition of LF anf VLF
                    
                    
                    #get the vessel ID present for all freq
                    vesselid_Rv = sleepdata['Rv'][stage].keys()
                    vesselid_h0 = sleepdata['h0'][stage].keys()

                    vesselid_oscLF = sleepdata['amp_area'+'LF'][stage].keys()
                    vesselid_oscVLF = sleepdata['amp_area'+'VLF'][stage].keys()


                    vesselIDs=intersect(vesselid_Rv,vesselid_h0,vesselid_oscLF,vesselid_oscVLF)  
                   
                
                    for vesselID in vesselIDs :
                    
                        jobname='disp'+'-d%.0e'%d+'-l%.0e'%lpvs+'-'+stage+'-'+'LFVLF'+'-id'+vesselID
                        print(jobname)
                        
                        #get properties
                        rv=sleepdata['Rv'][stage][vesselID]*1e-4 #cm
                        h0=sleepdata['h0'][stage][vesselID]*1e-4*corrfactorh0 #cm
                        rpvs=rv+h0
                        fi=1/sleepdata['period'+FB][stage][vesselID]
                        
                        meanarea=sleepdata['area'][stage][vesselID]*corrfactorarea
                    
                    
                       
                        fVLF=1/sleepdata['period'+'VLF'][stage][vesselID]
                        fLF=1/sleepdata['period'+'LF'][stage][vesselID]
                        
                        
                        # Compute area change directly from area oscilations
                        aVLF=sleepdata['amp_area'+'VLF'][stage][vesselID]/meanarea/2*corrfactorarea
                        aLF=sleepdata['amp_area'+'LF'][stage][vesselID]/meanarea/2*corrfactorarea
                        
                        
                          
                        
                        print('stage',stage)
                        print('vessel',vesselID)
                        print('FB',FB)
                        print('rpvs:',rpvs)
                        print('rv:',rv)
                        
                        print('arpvs:',amp_endfoot)
                        print('arv:',amp_lumen)
                        
                        print('h0:',h0)
                        print('f:',[fVLF,fLF])
                        print('a:',[aVLF,aLF])
    
                                          
    
                        #estimate best time parameters
                        
                        #the end of the simulation should at least cover 3 periods (to get a linear fit)
                        # and the time to diffuse in the transveral direction
                        #tend=max(5/fi,2*h0**2/d)
                        
                        tend=40

                        # we would like 4 output per period
                        toutput=1/fLF*4
                        
                        # we set dt as the toutput
                        dt=toutput
                        
                        # we would like at least 1000 time steps
                        while dt>(tend/1000) :
                            dt/=2

                        # we would like at max dt=5e-3
                        while dt>(5e-3) :
                            dt/=2
                        
                        #and at least 100 output over the whole simulation
                        while tend/toutput <100 :
                            toutput/=2
            
                        Noutput=tend/toutput
        
                        #estiamate of the max velocity
                        Umax=U(0, 0, aLF+aVLF, fLF, lpvs, Rv0=rv, h0=h0)
                        #estimate of the Peclet number
                        Pe=h0*Umax/d/2
                        # constraining velocity for the CFL condition
                        uconstrain=max(d/h0,Umax)
                        # number of cells in transversal direction
                        nr=8
                        dx=h0/nr
                        
                        CFL=(uconstrain/(dx/dt))
                        while CFL>0.7 :
                            dt/=2
                            CFL=(uconstrain/(dx/dt))
                            
                        print('# CFL :' ,CFL)
                        print('# Pe : %e'%Pe)
                        
                        if  Checkdt:
                            print_checkdt(tend,dt,Noutput,toutput)

                        slurmfile=jobname+'.slurm'

                        
                        #write the slurm file
                        # gaussian analysis
                        write_state_slurm(jobname, [fi],[ai],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='gaussian', c0valuePVS=1, c0valueSAS=0, sigma=sigma, sasbc='scenarioA', refineleft=False, aSMC=aSMC)                              

                        #update the list of slurm files to be launched in batch
                        slurmfiles.append(slurmfile)
                        
                    

            #write the batch file
            with open(seriedir+'batch','w') as f :
                i=1
                for slurmfile in slurmfiles :
                    f.write('sbatch '+slurmfile+' &\n')
                    if not(i%100):
                        f.write('sleep 1 \n')
                    i+=1
            f.close()



    
    #### Transport analysis 
    
    # serie name
    seriename='transportRandomWT10'
    
    # create a folder for the slurm files and the batch file

    if not os.path.exists('../output/'):
        os.makedirs('../output/')

    if not os.path.exists('../output/supercomputer/'):
        os.makedirs('../output/supercomputer/')


        
        
    seriedirroot='../output/supercomputer/'+seriename+'/' 

    for lpvs in spanlpvs:
        for d in spandiffusion :
            
            nl=int(lpvs/1e-4)
            nr=8
            sigma=2e-4
            
            # serie name
            serie=seriename+'-d%.0e'%d+'-l%.0e'%lpvs
    
            # create a folder for the slurm files and the batch file
            if not os.path.exists(seriedirroot+serie):
                os.makedirs(seriedirroot+serie)
            
            seriedir=seriedirroot+serie+'/'
    
            import time
            
            slurmfiles=[]
            # diffusion only
            jobname='diffusion'+'-d%.0e'%d+'-l%.0e'%lpvs
            print(jobname)
                
            #get properties
            rv=np.median(list(sleepdata['Rv']['baseline'].values()))*1e-4 #cm
            h0=np.median(list(sleepdata['h0']['baseline'].values()))*1e-4*corrfactorh0 #cm
            rpvs=rv+h0

            
            
                        
            tend=500

            toutput=0.5
                        
            dt=5e-2
                        
                        
            slurmfile=jobname+'.slurm'

            #write the slurm file
            # gaussian analysis
            #write_state_slurm(jobname, [0],[0],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='gaussian', c0valuePVS=1, c0valueSAS=0, sigma=sigma, sasbc='scenarioA', refineleft=False)                              
                        
            # intake
            write_state_slurm(jobname, [0],[0],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='uniform', c0valuePVS=0, c0valueSAS=1, sigma=sigma, sasbc='scenarioE', refineleft=False, sas=False)                              

            #update the list of slurm files to be launched in batch
            slurmfiles.append(slurmfile)
                        
            # superposition of LF anf VLF
                    
            for stage in spanstages :
                    #get the vessel ID present for all freq
                    vesselid_Rv = sleepdata['Rv'][stage].keys()
                    vesselid_h0 = sleepdata['h0'][stage].keys()

                    vesselid_oscLF = sleepdata['amp_area'+'LF'][stage].keys()
                    vesselid_oscVLF = sleepdata['amp_area'+'VLF'][stage].keys()


                    vesselIDs=intersect(vesselid_Rv,vesselid_h0,vesselid_oscLF,vesselid_oscVLF)  
                   
                
                    for vesselID in vesselIDs :
                        
                       
                          
                        #get properties
                        rv=sleepdata['Rv'][stage][vesselID]*1e-4 #cm
                        h0=sleepdata['h0'][stage][vesselID]*1e-4*corrfactorh0 #cm
                        rpvs=rv+h0
                        
                        meanarea=sleepdata['area'][stage][vesselID]*corrfactorarea
                    
                    
                        fVLF=1/sleepdata['period'+'VLF'][stage][vesselID]
                        fLF=1/sleepdata['period'+'LF'][stage][vesselID]
                        
                        
                        # Compute area change directly from area oscilations
                        aVLF=sleepdata['amp_area'+'VLF'][stage][vesselID]/meanarea/2*corrfactorarea
                        aLF=sleepdata['amp_area'+'LF'][stage][vesselID]/meanarea/2*corrfactorarea
                        
    

                        ### Only slow waves
                        jobname='disp'+'-d%.0e'%d+'-l%.0e'%lpvs+'-'+stage+'-'+'LFVLF'+'-id'+vesselID
                        print(jobname)                          
                        
                        print('stage',stage)
                        print('vessel',vesselID)
                        print('FB',FB)
                        print('rpvs:',rpvs)
                        print('rv:',rv)
                        
                        print('arpvs:',amp_endfoot)
                        print('arv:',amp_lumen)
                        
                        print('h0:',h0)
                        print('f:',[fVLF,fLF])
                        print('a:',[aVLF,aLF])
    
                                          
    
                        #estimate best time parameters
                        
                        #the end of the simulation should at least cover 3 periods (to get a linear fit)
                        # and the time to diffuse in the transveral direction
                        #tend=max(5/fi,2*h0**2/d)
                        
                        tend=200

                        # we would like 4 output per period
                        toutput=0.5
                        
                        # we set dt as the toutput
                        dt=1e-2
                        
            
                        Noutput=tend/toutput
        
                        #estiamate of the max velocity
                        Umax=U(0, 0, aLF+aVLF, fLF, lpvs, Rv0=rv, h0=h0)
                        #estimate of the Peclet number
                        Pe=h0*Umax/d/2
                        # constraining velocity for the CFL condition
                        uconstrain=max(d/h0,Umax)
                        # number of cells in transversal direction
                        dx=h0/nr
                        
                        CFL=(uconstrain/(dx/dt))
                        while CFL>0.7 :
                            dt/=2
                            CFL=(uconstrain/(dx/dt))
                            
                        print('# CFL :' ,CFL)


                        slurmfile=jobname+'.slurm'

                        # intake
                        write_state_slurm(jobname, [fVLF, fLF],[aVLF, aLF],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='uniform', c0valuePVS=0, c0valueSAS=1, sigma=sigma, sasbc='scenarioE', refineleft=False, sas=False)                              

                        #update the list of slurm files to be launched in batch
                        slurmfiles.append(slurmfile)
                        
                        
                    
                    #### with cardiac FB
                    
                    #get the vessel ID present for all freq
                    vesselid_Rv = sleepdata['Rv'][stage].keys()
                    vesselid_h0 = sleepdata['h0'][stage].keys()

                    vesselid_oscLF = sleepdata['amp_area'+'LF'][stage].keys()
                    vesselid_oscVLF = sleepdata['amp_area'+'VLF'][stage].keys()
                    vesselid_oscCard = sleepdata['amp_area'+'card'][stage].keys()


                    vesselIDs=intersect(vesselid_Rv,vesselid_h0,vesselid_oscLF,vesselid_oscVLF, vesselid_oscCard)  
                    
                    for vesselID in vesselIDs :
                        
                        #get properties
                        rv=sleepdata['Rv'][stage][vesselID]*1e-4 #cm
                        h0=sleepdata['h0'][stage][vesselID]*1e-4*corrfactorh0 #cm
                        rpvs=rv+h0
                        
                        meanarea=sleepdata['area'][stage][vesselID]*corrfactorarea
                    
                    
                        fCard=1/sleepdata['period'+'card'][stage][vesselID]
                        fVLF=1/sleepdata['period'+'VLF'][stage][vesselID]
                        fLF=1/sleepdata['period'+'LF'][stage][vesselID]
                        
                        
                        # Compute area change directly from area oscilations
                        aVLF=sleepdata['amp_area'+'VLF'][stage][vesselID]/meanarea/2*corrfactorarea
                        aLF=sleepdata['amp_area'+'LF'][stage][vesselID]/meanarea/2*corrfactorarea
                        
                        vcard=50e-4 #cm/s
                        wi=2*np.pi*fCard
                        aCard=vcard/wi/lpvs
                        
                        jobname='disp'+'-d%.0e'%d+'-l%.0e'%lpvs+'-'+stage+'-'+'LFVLFCard'+'-id'+vesselID
                        print(jobname)                          
                        
                        print('stage',stage)
                        print('vessel',vesselID)
                        print('FB',FB)
                        print('rpvs:',rpvs)
                        print('rv:',rv)
                        
                        print('arpvs:',amp_endfoot)
                        print('arv:',amp_lumen)
                        
                        print('h0:',h0)
                        print('f:',[fVLF,fLF,fCard])
                        print('a:',[aVLF,aLF,aCard])
    
                                          
    
                        #estimate best time parameters
                        
                        
                        tend=200

                        # we would like 4 output per period
                        toutput=0.5
                        
                        # we set dt as the toutput
                        dt=1e-2
                        
            
                        Noutput=tend/toutput
        
                        #estiamate of the max velocity
                        Umax=U(0, 0, aLF+aVLF, fLF, lpvs, Rv0=rv, h0=h0)
                        #estimate of the Peclet number
                        Pe=h0*Umax/d/2
                        # constraining velocity for the CFL condition
                        uconstrain=max(d/h0,Umax)
                        # number of cells in transversal direction
                        dx=h0/nr
                        
                        CFL=(uconstrain/(dx/dt))
                        while CFL>0.7 :
                            dt/=2
                            CFL=(uconstrain/(dx/dt))
                            
                        print('# CFL :' ,CFL)


                        slurmfile=jobname+'.slurm'

                        # intake
                        write_state_slurm(jobname, [fVLF, fLF,fCard],[aVLF, aLF, aCard],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='uniform', c0valuePVS=0, c0valueSAS=1, sigma=sigma, sasbc='scenarioE', refineleft=False, sas=False)                              

                        #update the list of slurm files to be launched in batch
                        slurmfiles.append(slurmfile)
                        

                        ##### Only cardiac
                        
                        jobname='disp'+'-d%.0e'%d+'-l%.0e'%lpvs+'-'+stage+'-'+'Card'+'-id'+vesselID
                        print(jobname)                          
                        
                        print('stage',stage)
                        print('vessel',vesselID)
                        print('FB',FB)
                        print('rpvs:',rpvs)
                        print('rv:',rv)
                        
                        print('arpvs:',amp_endfoot)
                        print('arv:',amp_lumen)
                        
                        print('h0:',h0)
                        print('f:',[fCard])
                        print('a:',[aCard])
    
                                          
    
                        #estimate best time parameters
                        
                        
                        tend=200

                        # we would like 4 output per period
                        toutput=0.5
                        
                        # we set dt as the toutput
                        dt=1e-2
                        
            
                        Noutput=tend/toutput
        
                        #estiamate of the max velocity
                        Umax=U(0, 0, aLF+aVLF, fLF, lpvs, Rv0=rv, h0=h0)
                        #estimate of the Peclet number
                        Pe=h0*Umax/d/2
                        # constraining velocity for the CFL condition
                        uconstrain=max(d/h0,Umax)
                        # number of cells in transversal direction
                        dx=h0/nr
                        
                        CFL=(uconstrain/(dx/dt))
                        while CFL>0.7 :
                            dt/=2
                            CFL=(uconstrain/(dx/dt))
                            
                        print('# CFL :' ,CFL)


                        slurmfile=jobname+'.slurm'

                        # intake
                        write_state_slurm(jobname, [fCard],[aCard],rv,h0, '"${USERWORK}/sleepoutput/'+serie+'"',seriedir+slurmfile, lpvs=lpvs,d=d,dt=dt, toutput=toutput,tend=tend,nl=nl, nr=nr,c0init='uniform', c0valuePVS=0, c0valueSAS=1, sigma=sigma, sasbc='scenarioE', refineleft=False, sas=False)                              

                        #update the list of slurm files to be launched in batch
                        slurmfiles.append(slurmfile)

                    

            #write the batch file
            with open(seriedir+'batch','w') as f :
                i=1
                for slurmfile in slurmfiles :
                    f.write('sbatch '+slurmfile+' &\n')
                    if not(i%100):
                        f.write('sleep 1 \n')
                    i+=1
            f.close()
  

    