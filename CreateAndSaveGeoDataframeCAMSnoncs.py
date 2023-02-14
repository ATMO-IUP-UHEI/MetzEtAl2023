#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Create a GeoDataFrame out of CAMS CO2  files

"""
Created on Tue May  4 18:33:52 2021

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

def CreateDF(DS,year,month):
    for d in range(len(DS)):                
        dCO2 = pd.DataFrame(data=DS.variables['XCO2'].values[d],index = DS.variables['latitude'].values,
                        columns =  DS.variables['longitude'].values, dtype='float' )
        dCO2 = pd.melt(dCO2.reset_index(), id_vars='index',
                       value_name ='CO2',var_name = 'Long')
        dCO2["Lat"] = dCO2['index']
        dCO2.insert(loc=1,column='Year',value= np.ones(len(dCO2["Lat"]))*year) 
        dCO2.insert(loc=1,column='Month',value= np.ones(len(dCO2["Lat"]))*month)
        #dCO2.insert(loc=1,column='Day',value= np.ones(len(dCO2["Lat"]))*np.floor(DS.variables['time'].values[d]/24))
        dCO2.insert(loc=1,column='Day',value= np.ones(len(dCO2["Lat"]))*15)

        if d == 0:
            df = dCO2.copy()
        else:
            df = df.append(dCO2) 

    return df

def CreateDataFrameCAMSnoncs(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num):
    savepath = datapath + "/CAMS/dataframes/"
    if not os.path.isfile(savepath+"DFD0_Monthlynoncs3_"+RegionName+str(Num)+".pkl"):
        print("Start reading data:")
        
        for num, filepath in enumerate(glob.glob(datapath + "/CAMS/MeanColumn/*surface*.nc")):      
            #print(filepath)
            DS = xr.open_mfdataset(filepath,combine='by_coords',concat_dim='None',decode_times =False)
            year = int(filepath[-9:-5])
            month = int(filepath[-5:-3])
            if num == 0:        
            #create Dataframe
                df3= CreateDF(DS,year,month)
                df = df3[(df3.Long >= Long_min)
                        & (df3.Long <= Long_max)
                        & (df3.Lat >= Lat_min)
                        & (df3.Lat <= Lat_max)]
            else:      
            #create Dataframe
                df3= CreateDF(DS,year,month)
                df2 = df3[(df3.Long >= Long_min)
                        & (df3.Long <= Long_max)
                        & (df3.Lat >= Lat_min)
                        & (df3.Lat <= Lat_max)]
                df = df.append(df2, ignore_index=True)
        
        df.drop_duplicates(keep = 'first')   
        
        df.to_pickle(savepath+"DFD0_Monthlynoncs3_"+RegionName+str(Num)+".pkl")
        print("finished reading data") 
    else:
        df = pd.read_pickle(savepath+"DFD0_Monthlynoncs3_"+RegionName+str(Num)+".pkl")
    #create date variable
    print("create timestamp")
    date = []
    for i in range(len(df.Lat)):
        
        #date.append(datetime.date(int(df.Year.values[i]),int(df.Month.values[i]),int(df.Day.values[i]))) 
        date.append(datetime.date(int(df.Year.values[i]),int(df.Month.values[i]),15)) 
    df.insert(loc=1,column='Date',value=date)
    
    df.to_pickle(savepath+"DF_Monthlynoncs3_"+RegionName+str(Num)+".pkl")
    #create GeoDataFrame
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Long, df.Lat))
    gdf.crs = {'init' :'epsg:4326'}
 
    if Num >= 900:
        Transcom = pd.read_pickle(savepath + "/Transcom_Regions.pkl")
        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]
   
    gdf.to_pickle(savepath+"GDF_Monthlynoncs3_"+RegionName+str(Num)+".pkl")

#main
if __name__=='__main__':
    Numm = 901
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)
    CreateDataFrameCAMSnoncs(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm) 
      