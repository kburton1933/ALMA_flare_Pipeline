#!/usr/bin/env python
# coding: utf-8

# In[3]:


import glob
from astropy import units as u
from astropy.time import Time
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
#import math as m
from astropy.io import fits
import pandas as pd


# In[3]:


### Create lc files with the appropriate scan information ###
# selects all of the light curve files in the directory #
def create_lc_with_scans():

    file1 = 'listobs_scan_info.csv' #file information
    #file2 = 'listobs_example.csv'
    light_curves = glob.glob("lightcurve*.csv") #grabs the current light curve files in the directory
    #light_curves = ['lightcurve_GJ581_07-Sep-2018_Obs5.csv']
    #'lightcurve_GJ581_01-Sep-2018_Obs3.csv',
    #'lightcurve_GJ581_04-Sep-2018_Obs4.csv',
    #'lightcurve_GJ581_07-Sep-2018_Obs5.csv',
    #'lightcurve_GJ581_07-Sep-2018_Obs6.csv',
    #'lightcurve_GJ581_08-Sep-2018_Obs7.csv',
    #'lightcurve_GJ581_09-Sep-2018_Obs8.csv',
    #'lightcurve_GJ581_10-Sep-2018_Obs9.csv',
   # 'lightcurve_GJ581_10-Sep-2018_Obs10.csv',
   # 'lightcurve_GJ581_14-Sep-2018_Obs11.csv',
    #'lightcurve_GJ581_21-Sep-2018_Obs12.csv']

    scan_info = Table.read(file1) #gets the info produces by scan_finder
    begin_time = np.array(scan_info['Begin Time JD'])
    end_time = np.array(scan_info['End Time JD'])


    #gets the observation number column for the file produced by listobs_manipulator.
    #obs_info = np.array(Table.read(file2)['Observation']) 
    scan_IDs = np.array(scan_info['ScanID'])


    #Loop through all of the light curve files in the directory, and add in the appropriate scan info 

    scan_ = []
    print(len(light_curves))
    for k in range(len(light_curves)):
        flux_info = np.array(Table.read(light_curves[k])) 
        scan = []
        flux_info = Table.read(light_curves[k])
        time_values = flux_info['# JDTime']
        flux_values = flux_info['Flux density [mJy]']
        flare_timeranges = flux_info['Timerange(UTC)'] 
        
        try: 
            val = np.array(Table.read(light_curves[k])['# JDTime'])
            print(len(val))
            for i in range(len(val)):
                for j in range(len(begin_time)):
                    if round(begin_time[j],5) <= round(val[i],5) < round(end_time[j],5): #may need to adjust this number if you get an inconsistent data column error 
                        scan.append(np.array(scan_info['ScanID'])[j]) 
        #create a new files w scan information
            if len(scan)> len(time_values):
                scan = scan[:len(time_values)]
            if len(scan)<len(time_values):
                scan.append(scan[-1]) #might get an error if the scan array is much shorter than the others
            print(len(time_values))
            print(len(scan))
            print(len(flux_values))
            print(len(flare_timeranges)) 
            master = Table([time_values, flux_values, flare_timeranges,scan],
                       names=('# JDTime', 'Flux density [mJy]', 'Timerange(UTC)','ScanID'))
            master.write(light_curves[k].replace('.csv','').replace('lightcurve','lc')+'_scan_info.csv',overwrite=True) 
        except:
           print('ugh')
           #df = pd.DataFrame([time_values, flux_values, flare_timeranges,scan])
           #difference = abs(len(flux_values) - len(scan))
           #print(difference)
           #for i in range(difference):
               #scan.append(scan[-1])
           #master = Table([df.fillna(0.0).values],
                       #names=('# JDTime', 'Flux density [mJy]', 'Timerange(UTC)','ScanID'))
           #master.write(light_curves[k].replace('.csv','').replace('lightcurve','lc')+'_scan_info.csv',overwrite=True)
    return scan_IDs

# In[5]:


#finds the 3 lowest flux integrations (excluding negatives) from each scan, calculates the 
#rms for each of these, and records

def select_choices():
    #open the newly created lc files, and concatenate all of the light curve info ###
    #including the flare timeranges, flux, and scans# 
    lc = glob.glob("lc*.csv")
    all_times = []
    all_scans = []
    all_flux = []

    for i in range(len(lc)):
        all_times.append(np.array(Table.read(lc[i])['Timerange(UTC)']))
        all_scans.append(np.array(Table.read(lc[i])['ScanID']))
        all_flux.append(np.array(Table.read(lc[i])['Flux density [mJy]']))

    all_times = np.concatenate(all_times)
    all_scans = np.concatenate(all_scans)
    all_flux = np.concatenate(all_flux)

    #loop through the different scans in scan_ID and have the code pick the 3 lowest flux values to 
    #compute an rms for. 

    #record the scan, the flux, and the timerange
  
    try:
        file1 = 'listobs_scan_info.csv' #file information
        scan_info = Table.read(file1) #gets the info produces by scan_finder
        scan_IDs = np.array(scan_info['ScanID'])
    except(FileNotFoundError):
        print('Listobs_scan_info.csv not found. Please run the Scan_finder.py part of the pipeline')


    indices = []
    for s in scan_IDs: 
        indices.append(np.where(all_scans==s)[0])
    #print(f'These are the scan_IDs {scan_IDs}')
    #print(indices)
    #print(len(indices))

    choice_index = []
    choices = []
    #print(len(all_flux))
    for i in range(len(indices)):
        m = np.where((all_flux[indices[i]] > 0))[0] #ignore where flux is 0, bc 0 indicates no data.
        choices.append(np.sort(all_flux[indices[i]][m])[:3]) #take the lowest 3 of the flux values
        choice_index.append(np.argsort(all_flux[indices[i]][m])[:3]) #gets the indices of the lowest 3 of the flux values
        #print(f'the choices are: {choices[i]} for scan {scan_IDs[i]}')
    

    timerange_choice = []
    for i in range(len(indices)):
        timerange_choice.append(all_times[indices[i]][choice_index[i]])
        
    no_scan_info = []
    for i in range(len(choices)): 
        if len(choices[i]) == 0: 
            print(f'No information for scan {scan_IDs[i]}. Removing the scan from further analysis')
            no_scan_info.append(i)

    std = np.array([np.std(s) for s in choices])
    
    
    scan_IDs = np.delete(scan_IDs, no_scan_info)
    timerange_choice = np.delete(timerange_choice, no_scan_info)
    choice_index = np.delete(choice_index, no_scan_info)
    choices = np.delete(choices, no_scan_info)
    std = np.delete(std, no_scan_info)
    
    print(scan_IDs)
    
    return scan_IDs, timerange_choice, choice_index, choices, std 

def find_rms_scan_candidates(vis,imsize,cell,region):
    scan_IDs, timerange_choice, choice_index, choices, std = select_choices()
    #now create a represenative rms file, record the stdev between the flux
    rep_rms = Table([timerange_choice, choices, scan_IDs,std],
                   names=('Timerange (UTC)','3 Lowest flux Per Scan','Scans','Standard Deviation in Flux'));
    #now pick out the flare candidates from the represenative rms file 

    rms_mean = []
    #print(len(scan_IDs))
    print(choices)
    for i in range(len(scan_IDs)):
       
       
        first= timerange_choice[i][0]
        second = timerange_choice[i][1]
        third = timerange_choice[i][2]
        
        try:

            #choose imagenames so that the files don't get rewritten 
            tclean(vis=vis, imagename='file1{0}'.format(i),
                   imsize=imsize, cell=cell,selectdata=True,
                   weighting='natural', niter=100,
                   interactive=False, timerange=first)
            tclean(vis=vis, imagename='file2{0}'.format(i),
                   imsize=imsize, cell=cell,selectdata=True,
                   weighting='natural', niter=100,
                   interactive=False, timerange=second)
            tclean(vis=vis, imagename='file3{0}'.format(i),
                   imsize=imsize, cell=cell,selectdata=True,
                   weighting='natural', niter=100,
                   interactive=False, timerange=third)
            
        except:
            print('tclean is not working')
            

        try:
            name1 = 'file1{0}.image'.format(i) 
            name2 = 'file2{0}.image'.format(i)
            name3 = 'file3{0}.image'.format(i)

            names = np.array([name1,name2,name3])
            rms= []

            #exportfits(imagename=name +'.image',fitsimage=imagename+'.fits')
            for n in names: 
            #compute an rms for all three of the time-ranges, then take the average of the rms
                rms.append(imstat(n,region=region)['rms'][0])


            rms_mean.append(rms)
        
        except: 
            print('something else has gone wrong')
            rms_mean.append('-')
        
    
 

    timerange_choice = [str(g) for g in timerange_choice]
    choices = [str(g) for g in choices]
    rms_mean = [str(g) for g in rms_mean]
    rep_rms = Table([timerange_choice, choices, scan_IDs,std,rms_mean],
                    names=('Timerange (UTC)','3 Lowest flux Per Scan','Scans','Standard Deviation in Flux','rms [mJy]'))
    rep_rms.write('rep_rms.csv',overwrite=True,format='ascii.csv')


# In[ ]:


def find_flare_candidates(integration,vis,vis_xx,vis_yy,directory):

    #create a dictionary to lookup the spectral window ID for each scan. The spectral window IDs can be found in the 
    #listobs_scan_info file that gets created from the scan finder script 

    #create a representative rm
    #create a representative rms from the file 
    rms_ = (Table.read('rep_rms.csv',format='ascii.csv')['rms [mJy]'])
    #rms_ = np.array([s for s in rms_])
    #rms_ = np.array([np.array(s) for s in rms_])
    
    m = [str(g) for g in rms_]
    m_ = [g.replace(' ',',').replace(',,',',') for g in m]


    rms_ = []
    for g in m_:
        try: 
            rms_.append(eval(g))
        except:
            rms_.append('-')

    rms_avg = []
    for s in rms_:
        try:
            rms_avg.append(np.mean(s)) #uses the average rms as the representative rms 
        except:
            rms_avg.append(0) 


    table = Table.read('rep_rms.csv' ,format= 'ascii.csv')

    #create a dicitonary to match up the scans with the rms. This will help decide which rms should be used
    #per integration in the light curve files 

    dictionary = {}

    keys = table['Scans']
    values = rms_avg

    for key, value in zip(keys, values):
        dictionary[key] = value
    print(dictionary)


    lc = glob.glob("lc*.csv")
    integration_time = integration

    all_times = []
    all_scans = []
    all_flux = []
    all_times_jd = []
    time_after_start = []

    for i in range(len(lc)):
        all_times.append(np.array(Table.read(lc[i])['Timerange(UTC)']))
        all_times_jd.append(np.array(Table.read(lc[i])['# JDTime']))
        all_scans.append(np.array(Table.read(lc[i])['ScanID']))
        all_flux.append(np.array(Table.read(lc[i])['Flux density [mJy]']))
        #time_after_start.append(np.arange(0,len(Table.read(lc[i])),integration_time)) 
        #this will work so long as the date format for the lcs stay same


    all_times = np.concatenate(all_times)
    all_scans = np.concatenate(all_scans)
    all_flux = np.concatenate(all_flux)
    all_times_jd = np.concatenate(all_times_jd)
    time_after_start = np.arange(0,len(all_times)*integration_time,integration_time)

    dates = []
    lc = glob.glob("lc*.csv")
    for l in lc: 
        for i in range(len(Table.read(l))):
            dates.append(l.split('_')[2] + "_"+ l.split('_')[3])
    dates = np.array(dates)


    #create an array with the corresponding representative rms information
    rms_array = []
    for i in range(len(all_flux)): 
        rms_array.append(dictionary[all_scans[i]])
    rms_array = 1000* np.array(rms_array) #convert the rms from Janskys to millijanskys 


    #grab all of the 3 sigma and above integrations as flare candidates. Save everything to file
    flare_indices = np.where(all_flux>= 3*rms_array)[0]
    rms_values = np.take(rms_array,flare_indices)
    possible_flares = np.take(all_flux,flare_indices)
    flare_scans = np.take(all_scans,flare_indices)
    flare_times = np.take(all_times,flare_indices)
    flare_times_jd = np.take(all_times_jd,flare_indices)
    dates2 = np.take(dates,flare_indices)
    time_after_start2 = np.take(time_after_start,flare_indices)
    significance = possible_flares/rms_values

    np.savetxt("flare_candidates_scans_all.csv".format(np.min(flare_scans),np.max(flare_scans)),np.transpose([flare_scans,possible_flares,
significance,rms_values,time_after_start2,dates2,flare_times,flare_times_jd]),delimiter=',',fmt='%s',header= 'Scan, Flux density [mJy],Significance,rms[mJy],Time_after_start,Date,Timerange(UTC), JDTime')

    flare_indices = np.where(all_flux>= 8*rms_array)[0]
    rms_values = np.take(rms_array,flare_indices)
    possible_flares = np.take(all_flux,flare_indices)
    flare_scans = np.take(all_scans,flare_indices)
    flare_times = np.take(all_times,flare_indices)
    flare_times_jd = np.take(all_times_jd,flare_indices)
    dates3 = np.take(dates,flare_indices)
    time_after_start2 = np.take(time_after_start,flare_indices)
    significance = possible_flares/rms_values

    print('possible_flares')
    print(len(possible_flares))
    print('flare_times')
    print(len(flare_times))
    print('flare_times_jd')
    print(len(flare_times_jd))
    print('dates:')
    print(len(dates))
    print('signficance:')
    print(len(significance))
    print('time_after_start:')
    print(len(time_after_start2))
    print('rms_values:')
    print(len(rms_values))
    print('flare_scans:')
    print(len(flare_scans))

    directory  = directory
    spw = Table.read('listobs_scan_info.csv')['Spwid']
    spw = np.array(spw)
    dict_spw = {}

    keys_ = (Table.read('listobs_scan_info.csv')['ScanID'])
    values_ = spw
    for key, value in zip(keys_, values_):
        dict_spw[key] = value
    print(dict_spw)

    vis = vis
    vis_xx = vis_xx
    vis_yy = vis_yy

    spw1_flux,spw2_flux, Exx, Eyy = [],[],[],[]
    print(len(possible_flares))

    for i in range(len(possible_flares)):
        #lookup which scan the integration is part of and use it to get the corresponding spw ID
        #records the spw ID of the upper and lower sidebands as strings for easy input into uvmod
        piece = dict_spw[flare_scans[i]].replace(' ',',').replace(',,',',').split(',')
        spw1 = piece[0] + ',' + piece[1]
        spw2 = piece[2] +','+ piece[3]
       
        print(spw1)
        try:
            uvmodelfit(vis=vis,sourcepar=[1,0,0],varypar=[True,False,False]
                       ,comptype='P',timerange=str(flare_times[i]),spw = spw1, outfile='spw1{0}.cl'.format(i))


            path = directory +'spw1{0}.cl'.format(i) #Change this to your directory!
            cl.open(path)
            cl.getcomponent(0)
            r = spw1_flux.append((cl.getfluxvalue(0)[0]* u.jansky).to(u.millijansky).value)
            #print("spw1 {0}/50 done". format(i))
            cl.close()
        except:
            r = spw1_flux.append(0)
 

        try: 
            uvmodelfit(vis=vis,sourcepar=[1,0,0],varypar=[True,False,False],
                       comptype='P',timerange=str(flare_times[i]),spw = spw2, outfile='spw2{0}.cl'.format(i)) 


            path =directory +'spw2{0}.cl'.format(i) #Change this to your directory!
            cl.open(path)
            cl.getcomponent(0)
            r = spw2_flux.append((cl.getfluxvalue(0)[0]* u.jansky).to(u.millijansky).value)
            #print("spw2 {0}/50 done". format(i))
            cl.close()
        except:
            r = spw2_flux.append(0)
            
        #get the XX polarization flux
        try:
            uvmodelfit(vis=vis_xx,sourcepar=[1,0,0],varypar=[True,False,False],
                       comptype='P',timerange=str(flare_times[i]), outfile='XX{0}.cl'.format(i)) 

            path =directory +'XX{0}.cl'.format(i) #Change this to your directory!
            cl.open(path)
            cl.getcomponent(0)
            r = Exx.append((cl.getfluxvalue(0)[0]* u.jansky).to(u.millijansky).value)
            cl.close()
        except:
            r = Exx.append(0)
           

        #record the YY polarization flux 
        try:
            uvmodelfit(vis=vis_yy,sourcepar=[1,0,0],varypar=[True,False,False],
                       comptype='P',timerange=str(flare_times[i]), outfile='YY{0}.cl'.format(i)) 

            path =directory+'YY{0}.cl'.format(i) #Change this to your directory!
            cl.open(path)
            cl.getcomponent(0)
            r = Eyy.append((cl.getfluxvalue(0)[0]* u.jansky).to(u.millijansky).value)
            cl.close()
        except:
            r = Eyy.append(0)
      


        #save all of this information to a file 
    print((Exx))
    print(len(Exx))
    print((Eyy))
    print(len(Eyy))
    print((spw1_flux))
    print(len(spw1_flux))
    print((spw2_flux))
    print(len(spw2_flux))
    
    np.savetxt("flare_candidates_scans_{0}-{1}.csv".format(np.min(flare_scans),np.max(flare_scans)),
               np.transpose([flare_scans,possible_flares,significance,rms_values,
                             time_after_start2,dates3,flare_times,flare_times_jd,spw1_flux,spw2_flux,
                             Exx,Eyy]), delimiter=',',fmt='%s',header= 
               'Scan, Flux density [mJy], Significance,rms[mJy],Time_after_start,Date,Timerange(UTC), JDTime, spw1,spw2,XX,YY')





integration = 1 #in seconds
vis = 'calibrated_final_avg.ms'
vis_xx = 'calibrated_final_avg_XX.ms'
vis_yy = 'calibrated_final_avg_YY.ms'
directory = '/Users/kianaburton/Desktop/MDwarfs/GJ191/' #Don't forget the forward slash in the directory!
imsize=[720,720]
cell=['0.1arcsec']
region = 'Ross_region.crtf'


#create_lc_with_scans()
#select_choices()
#find_rms_scan_candidates(vis, imsize,cell,region) #uses natural weighting 
find_flare_candidates(integration,vis,vis_xx,vis_yy,directory) #change varypar and sourcepar in the code if necessary


