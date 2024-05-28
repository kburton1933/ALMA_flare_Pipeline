#!/usr/bin/env python
# coding: utf-8

# In[233]:


from astropy.table import Table
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np
import math as m
import glob



def make_file(vis):
    vis = vis
    metadata = listobs(vis)

    string = str(metadata)
    string = string.replace('array','np.array') 

    m = open('dict_listobs.txt','w') #should have this save to a different name other than dict.txt
    m.write(string)
    m.close()

    #open the file back up
    with open('dict_listobs.txt') as f: #change this to dict_listobs.txt
        d = f.read()

    #turn string from file back into a dictionary

    #data = d.replace('array','np.array')
    data = eval(d)


    #from the list of scans, populate the array
    scans  = list(data.keys())

    header  = ['scanID','beginTime','EndTime','FieldName','SpwIds','IntegrationTime']
    begin,end, field,integration,spwid,scanid = [],[],[],[],[],[]


    for s in scans: #where scans is a string of the list of scans 
        if 'scan_' in s: #only works if the key is scan_
            begin.append(data[s]['0']['BeginTime'])
            end.append(data[s]['0']['EndTime'])
            field.append(data[s]['0']['FieldId'])
            integration.append(data[s]['0']['IntegrationTime'])
            spwid.append(data[s]['0']['SpwIds'])
            scanid.append(data[s]['0']['scanId'])

    #and of course listobs would give the scan information out of order... 

    sorted_index = np.argsort(scanid) 
    scanid = np.array(scanid)[sorted_index]
    begin = np.array(begin)[sorted_index]
    end = np.array(end)[sorted_index]
    field = np.array(field)[sorted_index]
    integration = np.array(integration)[sorted_index]
    spwid = (np.array(spwid)[sorted_index])
    spwid = [str(s).replace('[','').replace(']','') for s in spwid]


    # convert MJD to JD, so that those times are also in the file 
    # also create arrays for the date and utc time 

    bjd = Time(np.array(begin),format='mjd')
    ejd = Time(np.array(end),format ='mjd')
    begin_jd = bjd.jd
    end_jd = ejd.jd


    times_begin = bjd.to_value('datetime')
    times_end = ejd.to_value('datetime')

    date_begin = [str(s).split(' ')[0].replace('-','/') for s in times_begin]
    date_end = [str(s).split(' ')[0].replace('-','/') for s in times_end]

    t_begin = [str(s).split(' ')[1] for s in times_begin]
    t_end = [str(s).split(' ')[1] for s in times_end]


    #now find a way to record the corresponding observation ID and date, maybe just use this in conjunction with the 
    #listobs maker file that finds the appropriate listobs files

    #save all of this to a file
    file = Table([scanid, field, integration,spwid,begin,
           end, begin_jd, end_jd, date_begin,
           date_end, t_begin, t_end],names=('ScanID','field','integration_time'
                                            ,'Spwid','Begin Time MJD','End Time MJD ',
                                           'Begin Time JD','End Time JD','Begin Date',
                                           'End Date','Begin Time UTC','End Time UTC'))
    file.write('listobs_scan_info.csv',overwrite=True,format='ascii.csv') 



make_file('calibrated_final_avg.ms')

