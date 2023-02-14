#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Create a GeoDataFrame out of CT2019B CO2 data


"""
Created on Fri Apr  9 15:07:32 2021

@author: eschoema
"""

import datetime
#from datetime import timedelta
import read_remotec_out
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb
from functions import createPandasDateFrame 
import glob, os
import geopandas
import xarray as xr
import time
from RegionParam import getRegion

datapath = "."
savepath = "."

def CreateDF(DS):
    for l in range(len(DS.variables['co2'].values[0])):
        if l == 0:
            la = np.array(DS.variables['co2'].values[0][l]*DS.variables['air_mass'].values[0][l])
            la = la[np.newaxis,...]
        else:
            la2 = np.array(DS.variables['co2'].values[0][l]*DS.variables['air_mass'].values[0][l])
            la2 = la2[np.newaxis,...]
            la = np.append(la,la2, axis = 0)
    co2a = np.sum(la,axis = 0)/np.sum(DS.variables['air_mass'].values[0],axis = 0)
                    
    dCO2 = pd.DataFrame(data=co2a,index = DS.variables['latitude'].values,
                    columns =  DS.variables['longitude'].values, dtype='float' )
    dCO2 = pd.melt(dCO2.reset_index(), id_vars='index',
                   value_name ='CO2',var_name = 'Long')
    dCO2["Lat"] = dCO2['index']
    dCO2.insert(loc=1,column='Year',value= np.ones(len(dCO2["Lat"]))*DS.time_components.values[0][0]) 
    dCO2.insert(loc=1,column='Month',value= np.ones(len(dCO2["Lat"]))*DS.time_components.values[0][1])
    dCO2.insert(loc=1,column='Day',value= np.ones(len(dCO2["Lat"]))*15)

    return dCO2

def CreateDataFrameCTnoncs(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num):
    savepath = datapath + "/CT2019/dataframes/"
    print("Start reading data:")
    for num, filepath in enumerate(glob.glob(datapath + "/CT2019/MonthlyCT2019B/*.nc")):      
        DS = xr.open_mfdataset(filepath,combine='by_coords',concat_dim='None')
        if num == 0:        
        #create Dataframe
            df3= CreateDF(DS)
            df = df3[(df3.Long >= Long_min)
                    & (df3.Long <= Long_max)
                    & (df3.Lat >= Lat_min)
                    & (df3.Lat <= Lat_max)]
        else:      
        #create Dataframe
            df3= CreateDF(DS)
            df2 = df3[(df3.Long >= Long_min)
                    & (df3.Long <= Long_max)
                    & (df3.Lat >= Lat_min)
                    & (df3.Lat <= Lat_max)]
            df = df.append(df2, ignore_index=True)
    
    df.drop_duplicates(keep = 'first')   
    
    
    print("finished reading data") 
    #create date variable
    print("create timestamp")
    date = []
    for i in range(len(df.Lat)):
        date.append(datetime.date(int(df.Year.values[i]),int(df.Month.values[i]),int(df.Day.values[i]))) 
    df.insert(loc=1,column='Date',value=date)
    
    df.to_pickle(savepath+"DF_Monthlynoncs2_"+RegionName+str(Num)+".pkl")
    #create GeoDataFrame
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Long, df.Lat))
    gdf.crs = {'init' :'epsg:4326'}
 
    if Num >= 900:
        Transcom = pd.read_pickle("savepath + "/Transcom_Regions.pkl")
        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]
   
    gdf.to_pickle(savepath+"GDF_Monthlynoncs2_"+RegionName+str(Num)+".pkl")

if __name__=='__main__':
    Numm = 949
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)
    CreateDataFrameCTnoncs(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm)  
    