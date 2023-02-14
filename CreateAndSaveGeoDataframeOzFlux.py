#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 15:51:32 2021
Create pandas dataframes out of OzFlux data files
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
    try:
        Lat = DS.latitude.values
        Long = DS.longitude.values
    except:
        try:
            Lat = float(DS.latitude)
            Long = float(DS.longitude)
        except:
            Lat = np.nan
            Long = np.nan
    DS = DS.squeeze() # drops latitude and longitude columne which only contains one value

    try:
        d = {'Year':DS.Year.values[:],
            'Day':DS.Day.values[:], 
             'Month':DS.Month.values[:], 
             'Minute':DS.Minute.values[:],
             'Hour':DS.Hour.values[:],
             'Second':DS.Second.values[:],
             'Lat':np.ones(len(DS.Second.values[:]))*Lat,#DS.latitude.values,
             'Long':np.ones(len(DS.Second.values[:]))*Long}#DS.longitude.values}
        nanvec = np.ones(len(DS.Second.values))*np.nan
    except:
        time = DS.time.values
        years, months, days, hours, mins, secs = [],[],[],[],[],[]
        for times in time:
            years.append(int(str(times)[0:4]))
            months.append(int(str(times)[5:7]))
            days.append(int(str(times)[8:10]))
            hours.append(int(str(times)[11:13]))
            mins.append(int(str(times)[14:16]))
            secs.append(int(str(times)[17:19]))
        d = {'Year':years,
             'Day':days,
             'Month':months, 
             'Minute':mins,
             'Hour':hours,
             'Second':secs,
             'Lat':np.ones(len(DS.time.values))*Lat,#DS.latitude.values,
             'Long':np.ones(len(DS.time.values))*Long}#DS.longitude.values}    
        nanvec = np.ones(len(DS.time.values))*np.nan
    
    try:
         d.update({'Fc':DS.Fc.values[:]})
         d.update({'Fc_quality': DS.Fc_QCFlag.values[:]})
    except:
        try:
            d.update({'Fc':DS.Fco2.values[:]})
            d.update({'Fc_quality': DS.Fco2_QCFlag.values[:]})
        except:
            d.update({'Fc': nanvec})
            d.update({'Fc_quality': nanvec})
    try:
         d.update({'Ts':DS.Ts.values[:]})
         d.update({'Ts_quality': DS.Ts_QCFlag.values[:]})
    except:
        d.update({'Ts': nanvec})
        d.update({'Ts_quality': nanvec})
    try:
         d.update({'Ta':DS.Ta.values[:]})
         d.update({'Ta_quality': DS.Ta_QCFlag.values[:]})
    except:
        d.update({'Ta': nanvec})
        d.update({'Ta_quality': nanvec})
    try:
        d.update({'CO2':DS.Cc.values[:]})
        d.update({'CO2_quality': DS.Cc_QCFlag.values[:]})
    except:
        d.update({'CO2': nanvec})
        d.update({'CO2_quality': nanvec})
    try:
         d.update({'Precipitation': DS.Precip.values[:]})
    except:
        try:
            d.update({'Precipitation': DS.Rain_W2K.values[:]})
        except:
            d.update({'Precipitation': nanvec})
    try:
         d.update({'NEE':DS.NEE.values[:]})
         d.update({'NEE_quality': DS.NEE_QCFlag.values[:]})
    except:
        d.update({'NEE': nanvec})
        d.update({'NEE_quality': nanvec})
    try:
         d.update({'NEP':DS.NEP.values[:]})
         d.update({'NEP_quality': DS.NEP_QCFlag.values[:]})
    except:
        d.update({'NEP': nanvec})
        d.update({'NEP_quality': nanvec})    
    try:
        d.update({'GPP': DS.GPP.values[:]})
        d.update({'GPP_quality': DS.GPP_QCFlag.values[:]})
    except:
        d.update({'GPP': nanvec})
        d.update({'GPP_quality': nanvec})
    try:
        d.update({'Sws': DS.Sws.values[:]})
        d.update({'Sws_quality': DS.Sws_QCFlag.values[:]})
    except:
        d.update({'Sws': nanvec})
        d.update({'Sws_quality': nanvec})
    try:
        d.update({'ER_dark': DS.ER_dark.values[:]})
        d.update({'ER_dark_quality': DS.ER_dark_QCFlag.values[:]})
    except:
        d.update({'ER_dark': nanvec})
        d.update({'ER_dark_quality': nanvec})
    try:
        d.update({'ER_night': DS.ER_night.values[:]})
        d.update({'ER_night_quality': DS.ER_night_QCFlag.values[:]})
    except:
        d.update({'ER_night': nanvec})
        d.update({'ER_night_quality': nanvec})
    df = pd.DataFrame(data=d)

    return df

def CreateDataFrameOzFlux(station_id,level):
    folderPath = datapath + "/OzFlux/"
    if level == 'L6':
        folderPath = datapath + "/OzFlux/ASM_L6/"
    print("Reading data:" +station_id)
    filepath1 = folderPath+station_id+"*.nc"#_L6_CPartitioned.nc"
    
    for num, filepath in enumerate(glob.glob(filepath1)):    
        print(filepath)
        try:
            DS = xr.open_mfdataset(filepath,combine='by_coords',concat_dim='None')
        except:
            DS = xr.open_mfdataset(filepath,combine='by_coords',concat_dim='None', drop_variables = 'time')
        if num == 0:        
        #create Dataframe
            df= CreateDF(DS)
        else:      
        #create Dataframe
            df3= CreateDF(DS)
            df = df.append(df3, ignore_index=True)

    df.to_pickle(datapath + "/OzFlux/DataFrames/DF5_noDate_"+level+station_id+".pkl")
    #create date variable
    print("create timestamp")
    #date = []
    
    date = df.apply(lambda x: datetime.date(int(x.Year),int(x.Month),int(x.Day)),axis=1)
        
    df.insert(loc=1,column='Date',value=date)


    print('save data')
    df.to_pickle(datapath + "/OzFlux/DataFrames/DF5_"+level+station_id+".pkl")

#main
filenamelist = ['ADdry','AdelaideRiver','ADwet','ASM','BeaconFarm_IFR',
                'BeaconFarm_UUW','Boyagin1','Boyagin_understory','Calperum','CapeTribulation',
                'Collie', 'CowBay','CumberlandPlain','CumberlandMelaleuca','DalyPasture',
                'DalyRegrowth','DalyUncleared','Dargo','Digby','DryRiver',
                'Emerald','GatumPasture','Gingin','GWW','HowardSprings',
                'HowardUnderstory','Kopuatai','Litchfield','Longreach','Nimmo',
                'Otway','RDMF','Ridgefield','Riggs','Robson','Samsford',
                'SturtPlains','TTE','Tumbarumba','Wallaby','Warra',
                'Whroo','WombatStateForest','Yanco']
            #before 2023["WombatStateForest","Whroo","Tumbarumba","SturtPlains",
            #"Samford","Robson","Ridgefield","Litchfield",
            #"GWW","HowardSprings","HowardUnderstory",
            #"CowBay","CapeTribulation",
            #"Boyagin","BeaconFarm_IFR","BeaconFarm_UUW",
            #"Yanco","ASM","DalyUncleared","DryRiver","Gingin","TTE",
            #"Calperum","Warra","CumberlandPlain"]
            #["CumberlandPlain","Tumbarumba","SturtPlains",
            #"Samford","GWW","HowardSprings","HowardUnderstory",
            #"CowBay","CapeTribulation",
            #"Boyagin","BeaconFarm_IFR","BeaconFarm_UUW",
            #"Yanco","ASM","DalyUncleared","DryRiver","Gingin","TTE"]#,"CumberlandPlain","Calperum"]#"Warra","Robson","Ridgefield","Litchfield","WombatStateForest","Whroo"
for name in filenamelist:
    CreateDataFrameOzFlux(name,'') #'L6' for L6, '' for L3
