################################
### Module for the linescan traces data analysis
##################################

### Needed libraries
import numpy as np
import matplotlib.pyplot as plt 

import matplotlib
#matplotlib.use('TkAgg')
#plt.ioff()

import scipy as sp

from scipy import fftpack
import scipy.signal as sg


import os

from collections import namedtuple


def arc_length(x, y):
    """ Measure the arc length of a curve. Can be an indicator for noise. """
    npts = len(x)
    arc = np.sqrt((x[1] - x[0])**2 + (y[1] - y[0])**2)
    for k in range(1, npts):
        arc = arc + np.sqrt((x[k] - x[k-1])**2 + (y[k] - y[k-1])**2)

    return arc

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]

def remove_nan(y) :
    nans, indexes= nan_helper(y)
    y[nans]=np.interp(indexes(nans), indexes(~nans), y[~nans])
    return (y,nans)


def amp_analysis(signal,signal_smooth,time,cutoff1,cutoff2,view=True,export=[],outputdir='../output/amp-analysis/', outputformat='.pdf'):
    """
    Perform a peak to peak analysis of the amplitude and period of a timeserie in a given frequency band.

    Parameters
    ----------
    signal : TYPE np.array
        DESCRIPTION. the timeserie to analyse
    time : TYPE np.array
        DESCRIPTION. associated time vector
    cutoff1 : TYPE float
        DESCRIPTION.  low frequency cut off (Hertz)
    cutoff2 : TYPE float
        DESCRIPTION. high frequency cut off (Hertz)
    export : TYPE, optional string
        DESCRIPTION. The default is []. name of the output files if wanted
    outpurdir : TYPE, optional string
        DESCRIPTION. The default is '../output/amp-analysis/'. output directory used to save the figures (if export is not empty)

    Returns
    -------
    ampl : TYPE, np.array
        DESCRIPTION. array of the peak to peak amplitudes for the oscillations of the predominant time scale in the given frequency band.
    period: TYPE, np.array
        DESCRIPTION.  array of the peak to peak periods for the oscillations of the predominant time scale in the given frequency band.
    mean : TYPE, np.array
        DESCRIPTION.  array of the peak to peak mean values of the signal for the oscillations of the predominant time scale in the given frequency band.
        

    """
    time_step=time[1]-time[0]
    fs=1/time_step

        
    #spectrum
    frequency, Pxx_spec = sg.periodogram(signal_smooth,fs, 'hanning', scaling='spectrum')

    
    # peak to peak analysis
    
    # if the time windows is too short to cover the band frequency, we return nan values
    
    if time[-1]-time[0] < 1/cutoff2:
         ampl=[] # amplitude array
         period=[]
         t1=[]
         t2=[]
         mean=[]       
         
    else :
        
        # get main frequency
        I=np.where((frequency<cutoff2)&(frequency>cutoff1))[0]
        ifmax=np.argmax(Pxx_spec[I])
        fmax=frequency[I[ifmax]]
        
        if view :
        
            plt.figure()
            plt.plot(frequency[I],Pxx_spec[I])
            plt.plot(frequency[I[ifmax]],Pxx_spec[I[ifmax]],'*r')
            plt.yscale('log')
            plt.xlabel('frequency')
            plt.ylabel('power')
            plt.title('spectrum ' )
            if export :
                plt.savefig(outputdir+export+'-spectrum'+outputformat)
                plt.close()
        
            
        
        #peak detection 
        peaks_indices=sg.argrelextrema(signal_smooth, np.greater, order=int(1/fmax*fs/5))
        pit_indices=sg.argrelextrema(-signal_smooth, np.greater, order=int(1/fmax*fs/5))
        
        # amplitude detection
        #amplitude plot 
        ampl=[] # amplitude array
        period=[]
        t1=[]
        t2=[]
        mean=[]
        t_psys=time[peaks_indices] # time of psys
        t_pdia=time[pit_indices] #time of pdia
        rsys=signal_smooth[peaks_indices] #list of psys 
        rdia=signal_smooth[pit_indices] #list of pdia
        
        if view :
            plt.figure()
            plt.plot(time,signal_smooth,'-k')
        
        for ti in t_psys :
            #index_pdia=np.where((t_pdia>ti-1/fmax) & (t_pdia<ti))[0]
            index_pdia=np.where(((ti-t_pdia)<=1/cutoff1) & ((ti-t_pdia)>=0))[0]
    
            index_psys=np.where(ti==t_psys)[0]
            if(np.size(index_pdia)>0):
                
                a=rsys[index_psys][0]-rdia[index_pdia][-1]
                
                if a > 0 :
                    
                    
                    
                    # We consider this is an oscillation only if the previous valley has a value lower than the current peak.
                    ampl.append(a)
                    
                    t2.append(t_psys[index_psys[0]])
                    t1.append(t_pdia[index_pdia[-1]])
                    
                    if index_pdia[-1]<len(t_pdia)-1:  
                        # period is the time between two valleys
                        period.append((t_pdia[index_pdia[-1]+1]-t_pdia[index_pdia[-1]]))
                    else:
                        # exept the last one that we assume is twice the time between valley and peak
                        period.append((t_psys[index_psys[0]]-t_pdia[index_pdia[-1]])*2)
                    
                    Imean=np.where((time>=t_pdia[index_pdia[-1]])&(time<=t_pdia[index_pdia[-1]]+period[-1]))
                    try:
                        mean.append(np.mean(signal[Imean]))
                    except:
                        mean.append(signal[index_pdia[-1]])
                        continue
                    if view :
                        #plt.plot(t_psys,rsys,'ok')
                        #plt.plot(t_pdia,rdia,'ok')
                        plt.plot([t_psys[index_psys][0],t_pdia[index_pdia][-1]],[rsys[index_psys][0],rdia[index_pdia][-1]],'d:')
            
        
        if view :
                
            plt.xlabel('Time (s)')
            plt.ylabel('Radius(um)')
            plt.title('Peak to peak analysis , f=%.2e Hz'%fmax)  
            lim=[min(time),min(time)+min(8/fmax,max(time))]
            plt.xlim(lim)
            if export :
                plt.savefig(outputdir+export+'-amplitude'+outputformat)
                plt.close()
        
        #period=np.array(period).T[0]
        ampl=np.array(ampl)
        period=np.array(period)
        mean=np.array(mean)
        t1=np.array(t1)
        t2=np.array(t2)
        
    
            
    
    
    return ampl,period,t1,t2,mean

def extract_sequences(data):
    
    """ Get the lists of sequences of REM, NREM and IS sleep stages in a linescan database.
    A sequence objects has the attributes:
        seq.begin : time of beginning
        seq.end : time of ending
        seq-ibegin : index in the database of the beginning
        seq-iend : index in the database of the ending   
    """

    print('Extraction of sleep sequences ... \n')   
    REM_list = []
    NREM_list = []
    IS_list = []
    BASE_list= []
    AWAKENING_list=[]
    LOCO_list=[]
    WHISK_list=[]
    QUIET_list=[]


    REM = namedtuple('REM', ['begin', 'end','ibegin','iend'])
    NREM = namedtuple('NREM', ['begin', 'end','ibegin','iend'])
    IS = namedtuple('IS', ['begin', 'end','ibegin','iend'])
    BASE = namedtuple('BASE', ['begin', 'end','ibegin','iend'])
    AWAKENING = namedtuple('AWAKENING', ['begin', 'end','ibegin','iend'])
    QUIET = namedtuple('QUIET', ['begin', 'end','ibegin','iend'])
    LOCO = namedtuple('LOCO', ['begin', 'end','ibegin','iend'])
    WHISK = namedtuple('WHISK', ['begin', 'end','ibegin','iend'])

    REM_flag = False
    REM_value = 0
    REM_ivalue = 0

    NREM_flag = False
    NREM_value = 0
    NREM_ivalue = 0

    IS_flag = False
    IS_value = 0
    IS_ivalue = 0
    
    BASE_flag = False
    BASE_value = 0
    BASE_ivalue = 0
    
    QUIET_flag = False
    QUIET_value = 0
    QUIET_ivalue = 0
    
    LOCO_flag = False
    LOCO_value = 0
    LOCO_ivalue = 0
    
    WHISK_flag = False
    WHISK_value = 0
    WHISK_ivalue = 0
    
    AWAKENING_flag = False
    AWAKENING_value = 0
    AWAKENING_ivalue = 0

    value=0
    
    #minimum duration of the sequence requiered :
    minduration=1 #seconde(s)

    for i, row in data.iterrows():
        state = row['state']
        valuem1=value
        value = row['t']
        ivalue = i
        
        if (not REM_flag) and ((state == 'REM')|(state == 'Clean REM')):
            REM_flag = True
            REM_value = value 
            REM_ivalue=i
            
        elif (REM_flag) and not ((state == 'REM')|(state == 'Clean REM')):
            if value >= REM_value + minduration:
                new_seq = REM(REM_value, valuem1,REM_ivalue,ivalue-1)
                REM_list.append(new_seq)
            REM_flag = False
        
        if (not NREM_flag) and ((state == 'NREM')|(state == 'Clean NREM')):
            NREM_flag = True
            NREM_value = value 
            NREM_ivalue=i
            
        elif (NREM_flag) and not ((state == 'NREM')|(state == 'Clean NREM')):
            if value >= NREM_value + minduration:
                new_seq = NREM(NREM_value, valuem1,NREM_ivalue,ivalue-1)
                NREM_list.append(new_seq)
            NREM_flag = False
                
        if  (not IS_flag) and ((state == 'IS')|(state == 'Clean IS')):
            IS_flag = True
            IS_value = value 
            IS_ivalue=i
            
        elif (IS_flag) and not ((state == 'IS')|(state == 'Clean IS')):
            if value >= IS_value + minduration:
                new_seq = IS(IS_value, valuem1,IS_ivalue,ivalue-1)
                IS_list.append(new_seq)
            IS_flag = False 
                
        if  (not BASE_flag) and ((state == 'Vessel Baseline')|(state == 'Clean Vessel Baseline')):
            BASE_flag = True
            BASE_value = value 
            BASE_ivalue=i
            
        elif (BASE_flag) and not ((state == 'Vessel Baseline')|(state == 'Clean Vessel Baseline')):
            if value >= BASE_value + minduration:
                new_seq = BASE(BASE_value, valuem1,BASE_ivalue,ivalue-1)
                BASE_list.append(new_seq)
            BASE_flag = False 

        if (not QUIET_flag) and ((state == 'Quiet')|(state == 'Clean Quiet')):
            QUIET_flag = True
            QUIET_value = value 
            QUIET_ivalue=i
            
        elif (QUIET_flag) and not ((state == 'Quiet')|(state == 'Clean Quiet')):
            if value >= QUIET_value + minduration:
                new_seq = QUIET(QUIET_value, valuem1,QUIET_ivalue,ivalue-1)
                QUIET_list.append(new_seq)
            QUIET_flag = False 
                
                
        if (not LOCO_flag) and ((state == 'Locomotion')|(state == 'Clean Locomotion')):
            LOCO_flag = True
            LOCO_value = value 
            LOCO_ivalue=i
            
        elif (LOCO_flag) and not ((state == 'Locomotion')|(state == 'Clean Locomotion')):
            if value >= LOCO_value + minduration:
                new_seq = LOCO(LOCO_value, valuem1,LOCO_ivalue,ivalue-1)
                LOCO_list.append(new_seq)
            LOCO_flag = False
        
        if (not WHISK_flag) and ((state == 'Whisking')|(state == 'Clean Whisking')):
            WHISK_flag = True
            WHISK_value = value 
            WHISK_ivalue=i
            
        elif (WHISK_flag) and not ((state == 'Whisking')|(state == 'Clean Whisking')):
            if value >= WHISK_value + minduration:
                new_seq = WHISK(WHISK_value, valuem1,WHISK_ivalue,ivalue-1)
                WHISK_list.append(new_seq)
            WHISK_flag = False          
            
            
        if (not AWAKENING_flag) and ((state == 'Awakening')|(state == 'Awakening')):
            AWAKENING_flag = True
            AWAKENING_value = value 
            AWAKENING_ivalue=i
            
        elif (AWAKENING_flag) and not ((state == 'Awakening')|(state == 'Clean Awakening')):
            if value >= AWAKENING_value + minduration:
                new_seq = AWAKENING(AWAKENING_value, valuem1,AWAKENING_ivalue,ivalue-1)
                AWAKENING_list.append(new_seq)
            AWAKENING_flag = False         
                        
    print('- %i BASE sequence(s) found'%len(BASE_list))                  
    print('- %i REM sequence(s) found'%len(REM_list))
    print('- %i NREM sequence(s) found'%len(NREM_list))
    print('- %i IS sequence(s) found'%len(IS_list))
    print('- %i AWAKENING sequence(s) found'%len(AWAKENING_list))
    
    
    print('- %i LOCO sequence(s) found'%len(LOCO_list))
    print('- %i WHISK sequence(s) found'%len(WHISK_list))
    print('- %i QUIET sequence(s) found'%len(QUIET_list))
    
    plt.close('all')

                
    return REM_list, NREM_list, IS_list, BASE_list, LOCO_list,WHISK_list,QUIET_list, AWAKENING_list


#search peak
def searchpeak(signal, spantime,fmin,fmax) :
    """ returns the value of the frequency and its index of the peak between the min and max frequencies
        
        Paramters:
            signal (np.array):signal to study 
            spantime (float): time points where to evaluate the signal
            fmin (float): minimum frequency 
            fmax (float): maximum frequency 
        
        Returns:
            peak_freq0 (float): value of peak frequency 
            ipeak_freq0 (int): peak index
    """
    #fs=1/(spantime[1]-spantime[0])
    
    # sample_freq : frequency vector
    # power : corresponding power vector 
    sample_freq,power=periodogram(signal, spantime)
    
    # Find the peak frequency
    pos_mask = np.where((sample_freq > fmin)&(sample_freq <fmax))  # filtre pour les fequences entre min et max
    freqs = sample_freq[pos_mask] # on ne garde que les frequences du filtre
    ipeak_freq0 = pos_mask[0][0]+power[pos_mask].argmax() # indice du pic
    peak_freq0 = freqs[power[pos_mask].argmax()] # frequence du pic

    return peak_freq0, ipeak_freq0

    
#filter functions 

def lowpassfilter( signal, fs, cutoff):
    """ Generate low pass filter
        
        Parameters:
            signal (np.array): signal to filtered
            fs (float): sampling frequency 
            cutoff : choosen value for cutoff
        
        Return:
            filteredsingal (np.array): filtered signal
    """
    nyq = 0.5 *fs 
    normal_cutoff = cutoff / nyq
    b, a = sg.butter(5, normal_cutoff, btype='low', analog=False)  
    filteredsignal = sg.filtfilt(b, a, signal)
    return filteredsignal 


def highpassfilter( signal, fs, cutoff):
    """ Generate high pass filter
        
        Parameters:
            signal (np.array): signal to filtered
            fs (float): sampling frequency 
            cutoff : choosen value for cutoff
            
        Return:
            filteredsingal (np.array): filtered signal
    """
    nyq = 0.5 *fs 
    normal_cutoff = cutoff / nyq
    b, a = sg.butter(5, normal_cutoff, btype='high', analog=False)  
    filteredsignal = sg.filtfilt(b, a, signal)
    return filteredsignal 


def bandpassfilter( signal, fs, cutoff1, cutoff2):
    """ Generate band pass filter
        
        Parameters:
            signal (np.array): signal to filtered
            fs (float): sampling frequency 
            cutoff : choosen value for cutoff
            
        Return:
            filteredsingal (np.array): filtered signal
    """
    #Filtre les hautes frequences
    # cutoff
    cutoff=cutoff2
            
    nyq = 0.5 *fs 
    normal_cutoff = cutoff / nyq
    b, a = sg.butter(3, normal_cutoff, btype='low', analog=False)
    signal_smooth = sg.filtfilt(b, a, signal)

    #Filtre les basses frequences
    
    if cutoff1:
        # cutoff
        cutoff=cutoff1
                
        nyq = 0.5 *fs 
        normal_cutoff = cutoff / nyq
        b, a = sg.butter(3, normal_cutoff, btype='high', analog=False)
        signal_smooth = sg.filtfilt(b, a, signal_smooth)    
    
    return signal_smooth


def bandpassfilter2( signal, fs, cutoff1, cutoff2):
    """ Generate band pass filter
        
        Parameters:
            signal (np.array): signal to filtered
            fs (float): sampling frequency 
            cutoff : choosen value for cutoff
            
        Return:
            filteredsingal (np.array): filtered signal
    """
    nyq = 0.5 *fs 
    normal_cutoff1 = cutoff1 / nyq
    normal_cutoff2 = cutoff2 / nyq
    b, a = sg.butter(5, (normal_cutoff1,normal_cutoff2), btype='bandpass', analog=False)  
    filteredsignal = sg.filtfilt(b, a, signal)
    return filteredsignal


