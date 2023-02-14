#!/usr/bin/env python
# Create a GeoDataFrame out of TCCON files
#Author: E.-M. Schoemann
#Date 20.07.2020

import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob, os
import geopandas
import xarray as xr
import math

datapath = "."
savepath = "."

def CreateDF(DS):
    d = {'CO2': DS.xco2_ppm.values, 'CO2_error': DS.xco2_ppm_error.values,'CO': DS.xco_ppb.values,'CO_error':DS.xco_ppb_error.values,'Year':DS.year.values,'day_num':DS.day.values, 'hour_num':DS.hour.values, 'time':DS.time.values}
    df = pd.DataFrame(data=d)

    return df

def CreateDataFrameTCCON(station_id):
    print("Reading data:")
    filepath = "/home/atmo/users/mhaun/farewell_package/data/TCCON/data/gosat/netcdf/"+station_id+"*.public.nc"
    DS = xr.open_mfdataset(filepath,combine='by_coords',concat_dim='None',use_cftime=None,decode_times=False)
    #create Dataframe
    df= CreateDF(DS)

    #create date variable
    print("create timestamp")
    date = []
    datet = []
    Minute = []
    Month = []
    Day = []
    Hour = []
    Sec = []

    Monthdays = [0,31,28,31,30,31,30,31,31,30,31,30,31]
    Monthdays_ly = [0,31,29,31,30,31,30,31,31,30,31,30,31]
    if station_id in 'wg_db_ll':
        print('Caution: Assume every year to have 366 days but put all data of 29.2 on 28.2 for non gap year')

    for index, row in df.iterrows():
        if (row.Year % 4 == 0 and row.Year not in [1900,2100,2200,2300]) or (station_id in 'wg_db_ll'):
            mo = Monthdays_ly
        else:
            mo = Monthdays
        for l in range(1,len(mo)+1):
            if row.day_num <= sum(mo[0:l+1]):
                Day.append(int(row.day_num - sum(mo[0:l])))
                Month.append(l)
                break
            else:
                if l == len(mo):
                    print(row.Year)
                    print(row.day_num)
                    print(l)
                    print(sum(mo[0:l+1]))
                pass
        
        if row.hour_num < 0:
            Day[-1] = Day[-1] - 1
            if Day[-1] <= 0:
                 Month[-1] = Month[-1]-1
                 Day[-1] = mo[Month[-1]]
                 if Month == 0:
                     row.Year = row.Year -1
                     Month = 12
                     Day = 31
            
            Hour.append((math.floor(row.hour_num)))
            Minute.append((math.floor((row.hour_num-Hour[-1])*60)))
            Sec.append((math.floor((row.hour_num-Hour[-1]-(Minute[-1]/60))*3600)))
            Hour[-1] = 24 + Hour[-1]
        else:
            Hour.append(math.floor(row.hour_num))
            Minute.append(math.floor((row.hour_num-Hour[-1])*60))
            Sec.append(math.floor((row.hour_num-Hour[-1]-(Minute[-1]/60))*3600))

   
        try:
            date.append(datetime.date(int(row.Year),int(Month[-1]),int(Day[-1])))
            if Hour[-1] not in range(0,24):
                print(row.hour_num)
                print(Hour[-1])
            datet.append(datetime.datetime(int(row.Year),int(Month[-1]),int(Day[-1]),int(Hour[-1]),int(Minute[-1]),int(Sec[-1])))    
        except:
            #print(Day[-1])
            #print(Month[-1])
            date.append(datetime.date(int(row.Year),int(Month[-1]),int(28)))
            datet.append(datetime.datetime(int(row.Year),int(Month[-1]),28,int(Hour[-1]),int(Minute[-1]),int(Sec[-1])))   
    df.insert(loc=1,column='Date',value=date)
    df.insert(loc=1,column='DateTime',value=datet)
    df.insert(loc=1,column='Month',value=Month)
    df.insert(loc=1,column='Day',value=Day)
    df.insert(loc=1,column='Hour',value=Hour)
    df.insert(loc=1,column='Min',value=Minute)
    df.insert(loc=1,column='Sec',value=Sec)    

    print('save data')
    df.to_pickle(datapath + "/TCCON/DataFrames/DF2_"+station_id+".pkl")


