#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 08:41:27 2021
Create pandas dataframe out of OBsPack data
@author: eschoema
"""

import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob, os
import geopandas
import xarray as xr
import math
import pandas as pd

datapath = "."
savepath = "."

def CreateDF(DS):
    #print(DS.C_ppm.values)
    d = {'CO2': DS.value.values*1000000, 
         #'CO2_quality': DS.qcflag.values,
         'Obs_flag':DS.obs_flag.values,
         #'CO2_error': DS.value_std_dev.values,
         'time':DS.time.values,
         'Lat':DS.latitude.values,
         'Long':DS.longitude.values}
    nanvec = np.ones(len(DS.value.values))
    nanvec[:] = np.nan
    try:
        d.update({'CO2_error': DS.value_std_dev.values*1000000}) 
    except:
        d.update({'CO2_error': nanvec}) 
    try:
        d.update({'CO2_quality': DS.qcflag.values}) 
    except:
        d.update({'CO2_quality': nanvec}) 
    
    df = pd.DataFrame(data=d)
    

    return df

def CreateDataFrameObsPack(station_idlist):
    for station_id in station_idlist:
        print(station_id)
        folderPath = datapath + "/ObsPack_CO2_GLOBALVIEWplusv7.0/obspack_co2_1_GLOBALVIEWplus_v7.0_2021-08-18/data/nc/"
        print("Reading data:")
        filepath1 = folderPath+station_id+"*.nc"#_L6_CPartitioned.nc"
        
        
        for num, filepath in enumerate(glob.glob(filepath1)):    
            print(filepath)
            DS = xr.open_mfdataset(filepath,combine='by_coords',concat_dim='None',decode_times = False)
            if num == 0:        
            #create Dataframe
                df= CreateDF(DS)
            else:      
            #create Dataframe
                df3= CreateDF(DS)
                df = df.append(df3, ignore_index=True)
    
        df.to_pickle(datapath + "/ObsPack_CO2_GLOBALVIEWplusv7.0/DataFrames/DF1_"+station_id+".pkl")
        #create date variable
        print("create timestamp")
        #date = []
        
        year = df.apply(lambda x: (datetime.date(1970,1,1)+datetime.timedelta(seconds = int(x.time))).year,axis=1)
        df.insert(loc=1,column='Year',value=year)
        month = df.apply(lambda x: (datetime.date(1970,1,1)+datetime.timedelta(seconds = int(x.time))).month,axis=1)
        df.insert(loc=1,column='Month',value=month)
        day = df.apply(lambda x: (datetime.date(1970,1,1)+datetime.timedelta(seconds = int(x.time))).day,axis=1)
        df.insert(loc=1,column='Day',value=day)
        date = df.apply(lambda x: datetime.date(int(x.Year),int(x.Month),int(x.Day)),axis=1)
        df.insert(loc=1,column='Date',value=date)
    
    
        print('save data')
        df.to_pickle(datapath + "/ObsPack_CO2_GLOBALVIEWplusv7.0/DataFrames/DF1_"+station_id+".pkl")

#main
#co2_aia_aircraft-flask_2_representative.nc
CreateDataFrameObsPack(["co2_aia_aircraft-flask_2_representative",
                      "co2_ara_surface-flask_2_representative",
                      "co2_bhd_surface-flask_1_representative",
                      "co2_bhd_surface-flask_426_representative",
                      "co2_bhd_surface-insitu_15_baseline",
                      "co2_cfa_surface-flask_2_representative",
                      "co2_cgo_surface-flask_1_representative",
                      "co2_cgo_surface-flask_2_representative",
                      "co2_cgo_surface-flask_4_representative",
                      "co2_gpa_surface-flask_2_representative",
                      "co2_ota_surface-flask_2_representative",
                      ])
