#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 17:32:11 2023

@author: kianaburton
"""

# light curve maker versions 2 

from astropy import units as u
from astropy.time import Time
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import math as m


#should import from the listobs manipulator, which returns and array of values
#the code will produce light curves for the entire .ms. If you only want a few observations, 
#split the .ms, or modify the code to loop over the first few, or last few observations

#---------- INPUT PARAMETERS --------
star = 'GJ191' #name of the star.. this helps with naming the light curve files
vis = 'calibrated_final_avg.ms'
integration_time = 1 #integration time in seconds
listobsfile = 'calibrated_source.ms.listobs.txt'
directory = '/Users/kianaburton/Desktop/MDwarfs/GJ191/' 
#------------------------------------


def make_lc(star,vis,integration_time,listobsfile,directory):
    
    import Listobs_manipulator
    titles,dates_formated,times = Listobs_manipulator.get_info(star,listobsfile)
    
    for p in range(len(titles)):
        title = titles[p] #title for light curve and flux info. This determines the name that the light curve saves under
        vis = vis #visibility file for UVModelfit
        date = dates_formated[p][0] #date to read from
        integration_time = integration_time 
        obs_times = times[p] #scan times
        
        # don't forget to change the path variable down below to your directory!
        
        #--------------------------------------------------------------
    
        rms = 5.2 #rms in mJy
        
        #-----------converts the obs_time array into an array of time ranges that correspond to each integration, in the format that UVmodelfit requires--------
        
        holder,start,end = [],[],[]
        
        for time in obs_times:
          holder.append(time.split('-'))
        for i in range(len(holder)):    # ATTN!! previously len(holder)-1, but w 12m obs, this doesnt work. might need to switch back w aca obs
          start.append(date+'T'+holder[i][0])
          end.append(date+'T'+holder[i][1])
        
        
        #get time intervals to feed into UVModelfit....
        scan_values,field_values,array,new_array,timeranges,flux,good_times,sourcepar = [],[],[],[],[],[],[],[]
        flux.clear()
        good_times.clear()
        
        no_data = []
        empty_field = True 
        
        print(start)
        
        for j in range(len(start)):
            timebin = integration_time*u.second
            time = Time(start[j],scale='utc')
            end_time = Time(end[j],scale='utc')
        
            sourcepar = [1,0,0]
            print(sourcepar)
          
          
            while time <= end_time:
                new_time = time + timebin
                array.append(time)
                time = new_time
        
            for time in array:
                new_array.append(time.to_value(format='isot').replace('-','/').replace('T','/'))
        
            #get array in proper form to feed into UVModelfit
            for i in range(1,len(new_array)):
                timeranges.append(new_array[i-1] + '~' + new_array[i])
            
            print(len(timeranges))
        
            #run UVModelfit
            for i in range(len(timeranges)):
                counter = 0
                empty_field = False
        
                try:
                    print('----------------------------------------------')
                    print("------------------RUNNING UVMODELFIT ON INTEGRATION {0}/{1}--------------------".format(i,len(timeranges)))
                    uvmodelfit(vis=vis,sourcepar=sourcepar,varypar=[True,False,False]
                                ,comptype='P',timerange=timeranges[i],outfile='trash{0}.cl'.format(i))
        
                    path = directory +'trash{0}.cl'.format(i) #Change this to your directory!
                    cl.open(path)
                    cl.getcomponent(0)
                    r = flux.append((cl.getfluxvalue(0)[0]* u.jansky).to(u.millijansky).value)
                    good_times.append(timeranges[i])
                    cl.close()
        
        
                except RuntimeError:
                    print("Oops, no data for field. Saving flux as 0mJy")
                    flux.append(0)
                    good_times.append(timeranges[i])
                    no_data.append(timeranges[i])
                    continue
        
            timeranges.clear()
            array.clear()
            new_array.clear()
        
        print("!!COMPLETE!!")
        
        
        #---------- converts the time ranges to JD time for each integration, for easy light curve plotting --------
        
        #print(good_times)
        jd_times = []
        for n in good_times:  
            jd_times.append(n.split("~")[0]) #the first time in the timerange gets converted to the JD time
          
        values = []
        for n in jd_times:
            n = n.replace('/','-',2).replace('/','T',3)
            j = Time(n,scale='utc')
            values.append(j.jd)  #not generalized
        
        
        fluxes = []
        for i in range(len(flux)):
            fluxes.append(flux[i])
        
        
        
        #------------- sort values, and get the new indices of the those values to rearrange the flux densities! -------------
        
        sorted_values = sorted(values)
        s = np.array(values) 
        sort_index = np.argsort(s) 
        sorted_flux = []
        for i in range(len(fluxes)):
            sorted_flux.append(fluxes[sort_index[i]])
        print(len(sorted_flux))
        print(np.max(sorted_flux))
        
        sorted_good_times = []
        
        for i in range(len(fluxes)):
            sorted_good_times.append(good_times[sort_index[i]])
        
        scale = sorted_values[0]
        sorted_values_scaled = sorted_values - scale
        
        
        
        
        
        #--------- Save light curve information ------------- 
        
        np.savetxt("lightcurve_"+title+".csv",np.transpose([sorted_values,sorted_flux,sorted_good_times]),
                 delimiter=',',fmt='%s',header='JDTime, Flux density [mJy], Timerange(UTC)')
        
        
        
        #---------- Plots light curve----------
        
        
        plt.figure(figsize=(8,8))
        plt.rcParams["font.family"] = "serif"
        plt.rcParams["font.size"] = 18
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hlines(3*rms,np.min(sorted_values_scaled),np.max(sorted_values_scaled),color='gray',ls='--') 
        
        ax.tick_params(axis='x', which='major',pad=5)
        ax.errorbar(sorted_values_scaled,sorted_flux,yerr=rms,ecolor='black')
        ax.scatter(sorted_values_scaled,sorted_flux)
        ax.plot(sorted_values_scaled,sorted_flux)
        
        plt.ylabel("Flux density [mJy]")
        plt.xlabel("JD + {0}".format(round(scale,3)))
        plt.title(title) 

make_lc(star,vis,integration_time,listobsfile,directory)
    
