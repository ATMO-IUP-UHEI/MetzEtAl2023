#!/usr/bin/env python
# Create a GeoDataFrame out of OCO-2 CO2 files
#Author: E.-M. Schoemann
#Date 20.07.2020

import datetime
import read_remotec_out
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb
from functions import createPandasDateFrame 
import glob, os
import geopandas
import xarray as xr
from RegionParam import getRegion

datapath = "."
savepath = "."

def CreateDF(DS,DS2,DS3):
    d = {'CO2': DS.xco2.values, 
         'CO2_error': DS.xco2_uncertainty.values, 
         'Lat':DS.latitude.values,
         'Long':DS.longitude.values,
         'Year':DS.date.values[:,0],
         'Month':DS.date.values[:,1], 
         'Day':DS.date.values[:,2],
         'Hour':DS.date.values[:,3],
         'quality':DS.xco2_quality_flag.values, 
         'glint':DS3.glint_angle.values,
         'land':DS3.land_water_indicator.values}
         #'surface':np.abs(DS2.surface_type_flipped.values * DS3.land_water_indicator.values)}#for old 9r version
    df = pd.DataFrame(data=d)

    return df

def CreateDataFrameOCO(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num):
    print("Start reading data:")
    for k in range(2014,2022):
        print(k)
        for num, filepath in enumerate(glob.glob(datapath + "/OCO2_L2_Lite_FP10r_CO2/"+str(k)+"/*.nc4")):      
            DS = xr.open_mfdataset(filepath,combine='by_coords',concat_dim='None',use_cftime=None)
            #DS2 = xr.open_mfdataset(filepath,group = 'Auxiliary', combine='by_coords',concat_dim='None')
            DS3 = xr.open_mfdataset(filepath,group = 'Sounding', combine='by_coords',concat_dim='None')
            DS2 = 'empty'
            if num == 0 and k == 2014:        
            #create Dataframe
                df3= CreateDF(DS,DS2,DS3)
                df = df3[(df3.Long >= Long_min)
                        & (df3.Long <= Long_max)
                        & (df3.Lat >= Lat_min)
                        & (df3.Lat <= Lat_max)]
            else:      
            #create Dataframe
                df3= CreateDF(DS,DS2,DS3)
                df2 = df3[(df3.Long >= Long_min)
                        & (df3.Long <= Long_max)
                        & (df3.Lat >= Lat_min)
                        & (df3.Lat <= Lat_max)]
                df = df.append(df2, ignore_index=True)
        
    print("finished reading data") 

    #create date variable
    print("create timestamp")
    date = []
    for i in range(len(df.Year)):
        date.append(datetime.date(df.Year[i],df.Month[i],df.Day[i]))    
    df["Date"] = date
    

    df.to_pickle(datapath + "/OCO-2/DataFrames/DF2_"+RegionName+str(Num)+".pkl")
    #create GeoDataFrame
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Long, df.Lat))
    gdf.crs = {'init' :'epsg:4326'}

    if Num >= 900:
        Transcom = pd.read_pickle(savepath + "/Transcom_Regions.pkl")
        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]

    gdf.to_pickle(datapath + "/OCO-2/DataFrames/GDF2_"+RegionName+str(Num)+".pkl")


if __name__=='__main__':
    Numm  = 949
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)
    CreateDataFrameOCO(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm)
