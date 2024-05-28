#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 15:06:13 2024

@author: kianaburton
"""

from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import math as m
import glob
import uncertainties.unumpy as unumpy  
import uncertainties.umath
from astropy import constants as const
from astropy import units as u
from uncertainties import ufloat


def sort_flare_events(file, upper,lower,integration_time,err,par,name,band_freq):
    flux,sigs,time,timejd,spw1,spw2,XX,YY,scan = [],[],[],[],[],[],[],[],[]
    flux_info = Table.read(file)
    sigs = (np.array(flux_info['Significance']))
    flux = (np.array(flux_info['Flux density [mJy]']))
    time = (np.array(flux_info['Timerange(UTC)']))
    timejd = (np.array(flux_info['JDTime']))
    #spw1 = (np.array(flux_info['spw1'])) #flux in the lower sideband 
    #spw2 = (np.array(flux_info['spw2'])) #flux in the upper sideband 
    #XX=(np.array(flux_info['XX'])) #flux in the XX polarization
    #YY=(np.array(flux_info['YY'])) #flux in the YY polarization
    scan=(np.array(flux_info['# Scan'])) #scan for the flux 
    time_after_start = np.array(flux_info['Time_after_start']) 
    dates = np.array(flux_info['Date']) #uncomment this once the code gets fixed
    
    
    
    # read this in from the file 

    event = []
    
    for i in range(len(time_after_start)-1):
        if abs((time_after_start[i+1]) - time_after_start[i]) <=integration_time:
            num1 = time_after_start[i]
            num2 = time_after_start[i+1]
            event.append(num1)
            event.append(num2) 
    
    num = 1
    event_num = [1]
    event = np.unique(event)
    
    for i in range(len(time_after_start)-1):
        if abs(time_after_start[i+1] -time_after_start[i]) <= integration_time:
            num = num
        else: 
            num = num + 1
            
        event_num.append(num)
    
    
    #this code can be moved to the flare_analysis since the event_num will be recorded in the files
    #the code below can be modified so that np.where(g is max) is used to get the the corresponding times, 
    #spw, and polarization flux
    
    d = np.unique(event_num,return_index=True)[1]
    groups = np.split(flux,d) #the possible flares grouped up by their integration 
    
    peaks = []
    
    for g in groups: 
        try: 
            peaks.append(np.argmax(g))
        except:
            continue       
    
    num_flare_events = len(np.unique(event_num)) #total flaring events - without double counting any
    num_all_flares = len(event_num) #total number of flares
    num_flares_multi_int = len(event) #total number of multi-integration flares 
    
    num_flares_single_int = len(np.unique(event_num))- len(event) #total single integration flares
        
        
    g_sig = np.split(sigs,d)[1:len(d)+1]
    g_flux = np.split(flux,d)[1:len(d)+1]
    g_time = np.split(time,d)[1:len(d)+1]
    g_jd = np.split(timejd,d)[1:len(d)+1]
    g_spw1 = np.split(spw1,d)[1:len(d)+1]
    g_spw2 = np.split(spw2,d)[1:len(d)+1]
    g_XX = np.split(XX,d)[1:len(d)+1]
    g_YY = np.split(YY,d)[1:len(d)+1]
    g_scan = np.split(scan,d)[1:len(d)+1]
    g_time_after = np.split(time_after_start,d)[1:len(d)+1]
    g_dates = np.split(dates,d)[1:len(d)+1]
    
    sigs_peak, flux_peak,time_peak,timejd_peak,spw1_peak = [],[],[],[],[]

    spw2_peak,XX_peak,YY_peak,scan_peak,time_after_start_peak = [],[],[],[],[]
    
    dates_peak  = []
    for i in range(len(g_sig)):
        sigs_peak.append(g_sig[i][peaks[i]])
        flux_peak.append(g_flux[i][peaks[i]])
        time_peak.append(g_time[i][peaks[i]])
        timejd_peak.append(g_jd[i][peaks[i]])
        #spw1_peak.append(g_spw1[i][peaks[i]])
        #spw2_peak.append(g_spw2[i][peaks[i]])
        #XX_peak.append(g_XX[i][peaks[i]])
        #YY_peak.append(g_YY[i][peaks[i]])
        scan_peak.append(g_scan[i][peaks[i]])
        time_after_start_peak.append(g_time_after[i][peaks[i]])
        dates_peak.append(g_dates[i][peaks[i]])
        
    integrations = [len(g) for g in g_sig]#This is an array that tells us how many integrations each flaring event has. 
    #Will be useful to seperate the data by this and pick out the multi integration flares 
    
        
    
    
    #Perform the luminosity calculations 
    distance = par*u.parsec 
    distance  = distance.to(u.centimeter)  #distance to star in cm 
    
    #calculate the peak luminosity#
    ################################
    flux_uncert = np.linspace(0.55,0.55,len(g_flux)) #define some sort of uncertainty since uvmod doesn't give us one
    
    F = flux_peak *u.milliJansky  
    F = F.to(u.jansky) #convert flux peak to jansky 
    F = F.to(u.erg/u.second/u.Hz/u.cm**2) #convert janskys to erg/second/Hz/cm2
    F_uncertainty = flux_uncert* u.milliJansky
    F_uncertainty = F_uncertainty.to(u.jansky)
    F_uncertainty = F_uncertainty.to(u.erg/u.second/u.Hz/u.cm**2)
    
    
    L_peak = F * 4*m.pi*distance**2 #the peak luminosity is the F*pi*distance to star^2
    
    x = unumpy.uarray((F.value,F_uncertainty.value))
    y=x*4*m.pi*distance.value**2
    #y is the array of uncertainities and values for the luminosity
    ###################################
    
    #L_formated = [format(l.value/(10**13),'.1E') for l in L_peak]
    
    #calculate the energy of the flare# 
    ##################################
    #energy is the luminosity of the flare times the seconds 
    
    #step 1- find the total luminosity of each event 
    Energy_erg = []
    freq = band_freq  #frequency in hertz (1.3 mm) this will chnage depending on the band 
    
    for g in g_flux: 
        G = g *u.milliJansky  
        G = G.to(u.jansky) #convert flux peak to jansky 
        G = G.to(u.erg/u.second/u.Hz/u.cm**2) #convert janskys to erg/second/Hz/cm2
        G = G.value
        
        length = len(G)
        Energy_ = ([g * 4*m.pi*distance**2 * length * integration_time * freq for g in G]) 
        Energy_ = [q.value for q in Energy_]
        Energy_sum = np.sum(Energy_)
        Energy_erg.append(Energy_sum)
        #the energy of each event in erg. Multiplied the total luminosity of each event by 
        #total time of the event, and the frequency of the event. 
        #Energy_erg.append(G_sum* 4*m.pi*distance**2 * len(G) * integration_time * freq * u.Hz * u.second)
        
    

    flare_properties = Table([flux_peak,L_peak/10**13,Energy_erg,integrations,sigs_peak,dates_peak,time_after_start_peak,
                         time_peak,timejd_peak,scan_peak],
               names=('Peak Flux mJy', 'Peak Luminosity (10^13 erg/Hz s)','Energy (erg)','num of ints'
                     ,'significance','Date','Time after Start','Time(UTC)','Time jd','Scan'))

    flare_properties.write(name,overwrite=True,format='ascii.csv')
    
    
#######--CHANGE THESE PARAMETERS--##########
#central frequencies can be found in listobs. Use the average of the upper and lower spectral windows
file = 'flare_candidates_scans_all.csv' #the file to read from - created in previous pipeline step 
upper = 240 #245.6 #GHz 
lower = 225 #229.6 #GHz 
#lower 4,6 spw = 225 GHz and upper 8,10 or 16,18 = 240 Ghz spectral index calculation is lower/upper
integration_time = 1
err = 0.55 #choosing a random error in mjy for the alpha and QI calculations
par = 2.94 #distance to star in parsecs 
name = 'flare_properties_Ross154.csv'
band_freq = 225407863158

sort_flare_events(file, upper, lower, integration_time, err, par, name,band_freq)
############################################

