#!/usr/bin/env python
# coding: utf-8

# In[643]:


from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import math as m
import glob
import re


# In[644]:


def monthToNum(shortMonth):
    return {
            'Jan': '01',
            'Feb': '02',
            'Mar': '03',
            'Apr': '04',
            'May': '05',
            'Jun': '06',
            'Jul': '07',
            'Aug': '08',
            'Sep': '09', 
            'Oct': '10',
            'Nov': '11',
            'Dec': '12'
    }[shortMonth]


# In[653]:


#read in the listobs file.

def get_info(star,listobsfile):

    #change this to check if a certain file is already there, if not, run listobs, 
    #and create the necessary file

    r = listobsfile
    with open(r) as f:
        lines = [line for line in f]


    import re

    #takes the part of the listobs file from the first observationID, to the line '(nRows = Total number...)'
    #since the necessary information should be between these lines 
    new_line = []
    arr = []
    for i in range(len(lines)):
        if 'ObservationID = 0' in lines[i]:
            start = i
        if '(nRows = Total number of rows per scan)' in lines[i]:
            end = i
    data = lines[start:end]

    o = []
    for i in range(len(data)):
        if data[i] == '\n':
            o.append(i)

    #grabs the timeranges by looking for data of the format '00:00:00.0 - 00:00:00.0'
    result = []
    num = []
    lines = np.array(lines)
    for i in range(len(lines)):
        result.append(re.findall("\w+\:\w+\:\w+\.\w - \w+\:\w+\:\w+\.\w", lines[i]))
        num.append(i)


    #records where the timerange data is, and filters out the empty arrays
    location = []
    for i in range(len(result)):
        if len(result[i]) != 0:
            location.append(num[i])


    good_data = lines[location];

    #records the dates, observations

    #get the date information 
    dates = []
    for i in range(len(good_data)):
        dates.append(re.findall("\w+\-\w+\-\w+", good_data[i]))

    new_dates = []
    observations = []

    dates_info = []
    for i in range(len(dates)):
        if len(dates[i]) != 0:
            dates_info.append(i)

    #to split arrays based on index, can just use np.split. 
    n = 0
    for i in range(len(dates)):
        if len(dates[i]) != 0:
            new_dates.append(dates[i])
        else:
            new_dates.append(new_dates[i-1])
        if i in dates_info:
            n = n+1
            observations.append('Obs{0}'.format(n))
        else:
            observations.append('Obs{0}'.format(n))

    #now record the timeranges, etc.
    timeranges = np.concatenate(result)
    timeranges = list(timeranges)

    new_timeranges =[]
    for n in timeranges:
        n = n.replace(' ','',1).replace(' ','',2)
        new_timeranges.append(n)

    #now record the scan, field, and spwIds info
    #as of now, I'm not able to record the fieldID and spwIds since the file structures change 
    scan, fieldID, spwIds = [],[],[]
    for i in range(len(good_data)):
            p = good_data[i].split('  ')
            if len(p[9]) == 0:
                try:
                    scan.append(p[3])
                except IndexError:
                    scan.append('-')
                    continue    
                try:
                    fieldID.append(p[6])
                except IndexError:
                    fieldID.append('-')
                    continue   
            else:
                scan.append('-')
                fieldID.append('-')

    master = Table([new_dates, observations, new_timeranges,
                    scan,fieldID],
                   names=('Date', 'Observation', 'Timerange (UTC)','scan','field ID'))

    master.write('listobs_example.ecsv',overwrite=True) 
    
    d = Table.read('listobs_example.ecsv')

    ob = np.array(d['Observation'])

    if dates_info[0] == 0:
        new_dates_info = np.delete(dates_info,0)
    else: 
        new_dates_info = dates_info

    obs_times = np.split(ob,new_dates_info)
    date = np.split(np.concatenate(np.array(d['Date'])),new_dates_info)

    titles = []
    star = star
    for i in range(len(obs_times)):
        titles.append(star + "_" + date[i][0] +"_" +obs_times[i][0])
    titles = np.array(titles)

    times = np.split(np.array(d['Timerange (UTC)']), new_dates_info)

    #changes the dates into a readable format )(YYYY-MM-DD)
    dates_formated = []
    for n in np.concatenate(date): 
        split_ = n.split('-')
        month_num = monthToNum(split_[1])
        dates_formated.append(split_[2]+'-'+month_num+'-'+split_[0])
    dates_formated = np.split(dates_formated,new_dates_info)

    return titles,dates_formated,times



#star = 'Prox'
#titles,dates_formated,times = get_info(star,'july_16_17_listobs.txt')


# In[ ]:




