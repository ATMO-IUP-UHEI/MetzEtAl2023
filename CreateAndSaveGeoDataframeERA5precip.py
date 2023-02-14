#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 14:03:05 2021
Create pandas dataframes from ERA5 precipitation data
@author: eschoema
"""

import datetime
#from datetime import timedelta
import read_remotec_out
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb
from functions import getNumDayOfMonth, createPandasDateFrame, getReferenceDate, getReferencesDateDay 
import glob, os
import geopandas
import xarray as xr
import time
from shapely.geometry import Polygon
import math
from skimage.measure import block_reduce
from RegionParam import getRegion

#IMPORTANT see: https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Monthlymeans

datapath = "."
savepath = "."

def CreateDF(DS):
    #for each timestep
    for d in range(len(DS['tp'])):
        nam = 'tp'
        dfyCO2 = pd.DataFrame(data=block_reduce(
                                    DS.variables[nam][d].values,
                                    block_size=(4,4),
                                    func=np.mean),
                                index = list(np.array(range(-95,-505,-10))/10), 
                                columns = list(np.array(range(1125,1805,10))/10),
                                dtype='float' )
        dfyCO2 = pd.melt(dfyCO2.reset_index(), id_vars='index',
                                    value_name =nam,var_name = 'Long')
        dfyCO2["Lat"] = dfyCO2['index']
        dfyCO2.insert(loc=1,column='Year',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[0:4]))        
        dfyCO2.insert(loc=1,column='Month',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[5:7]))
        dfyCO2.insert(loc=1,column='Day',value= np.ones(len(dfyCO2["Lat"]))*int(str(DS.time.values[d])[8:10]))
        if d == 0:
            df = dfyCO2.copy()
        else:
            df = df.append(dfyCO2)                         
    return df

def CreateDataFrameERA5P(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num):
       
    if not os.path.isfile(datapath + "/ERA5/Australia/DataFrames/DF1_"+RegionName+str(Num)+".pkl"):
        # check if dataframe without time variable already has been created 
        if not os.path.isfile(datapath + "/ERA5/Australia/DataFrames/DFD0_"+RegionName+str(Num)+".pkl"):
            print("Start reading data:")
            filepath = datapath + "/ERA5/Australia/MonthlyPrecipitation.nc"               
            DS = xr.open_mfdataset(filepath, combine='by_coords',concat_dim='None')#decode_times=False)
   
            #create Dataframe
            df3= CreateDF(DS)
            # preselect region as rectangle
            df = df3[(df3.Long >= Long_min)
             & (df3.Long <= Long_max)
             & (df3.Lat >= Lat_min)
             & (df3.Lat <= Lat_max)]
        
            print("finished reading data") 
            df.to_pickle(datapath + "/ERA5/Australia/DataFrames/DFD0_"+RegionName+str(Num)+".pkl")
        else:
            df = pd.read_pickle(datapath + "/ERA5/Australia/DataFrames/DFD0_"+RegionName+str(Num)+".pkl")
        #create date variable
        print("create timestamp")
        date = df.apply(lambda x: datetime.date(int(x.Year),int(x.Month),int(x.Day)),axis=1)
        df.insert(loc=1,column='Date',value=date)
        DayPMonth = df.apply(lambda x: getNumDayOfMonth(int(x.Year),int(x.Month)),axis=1) 
        df.insert(loc=1,column='NumberDays',value=DayPMonth) 
        df.insert(loc=1,column='monthlytp',value=df.tp*df.NumberDays) 
        df.to_pickle(datapath + "/ERA5/Australia/DataFrames/DF1_"+RegionName+str(Num)+".pkl")
    #create GeoDataFrame
    else:
        df = pd.read_pickle(datapath + "/ERA5/Australia/DataFrames/DF1_"+RegionName+str(Num)+".pkl") 


    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Long, df.Lat),crs = 4326)
 
    if Num >= 900:
        Transcom = pd.read_pickle(savepath + "/Transcom_Regions.pkl")
        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]
   
    gdf.to_pickle(datapath + "/ERA5/Australia/DataFrames/GDF1_"+RegionName+str(Num)+".pkl")

Numm = 949
RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)
CreateDataFrameERA5P(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm)
