#!/usr/bin/env python
# Create a GeoDataFrame out of ACOS/GOSAT files
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

datapath = "."

gdf_R = pd.read_pickle(datapath + "/GOSAT_Markus/dataframes/GDF09_19_6.pkl")

def CreateDF(DS,DS2,DS3):
    d = {'CO2': DS.xco2.values, 
         'CO2_error': DS.xco2_uncertainty.values, 
         'Lat':DS.latitude.values,
         'Long':DS.longitude.values,
         'Year':DS.date.values[:,0],
         'Month':DS.date.values[:,1], 
         'Day':DS.date.values[:,2],
         'Hour':DS.date.values[:,3], 
         'Min':DS.date.values[:,4],
         'Sec':DS.date.values[:,5],
         'CO2_uncorr': DS2.xco2_raw.values,
         'quality': DS.xco2_quality_flag.values,
         'glint': DS3.glint_angle.values,
         'land_frac': DS3.land_fraction,
         'DWS':DS2.dws.values,
         'delGrad':DS2.co2_grad_del.values,
         'psurf':DS2.psurf.values,
         'dpfrac':DS2.dpfrac.values,
         'albedo_sco2':DS2.albedo_sco2.values, 
         'gain':DS3.gain.values}
    df = pd.DataFrame(data=d)

    return df

def CreateDataFrameACOS(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num):
    print("Start reading data:")
    for num, filepath in enumerate(glob.glob(datapath + "/ACOS/ACOS_L2_Lite_FP.9r/*.nc4")):      
        DS = xr.open_mfdataset(filepath,combine='by_coords',concat_dim='None',use_cftime=None)
        DS2 = xr.open_mfdataset(filepath,group = 'Retrieval',combine='by_coords',concat_dim='None',use_cftime=None)
        DS3 = xr.open_mfdataset(filepath,group = 'Sounding',combine='by_coords',concat_dim='None',use_cftime=None)
        if num == 0:        
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
    inRemotec = []
    for i in range(len(df.Year)):
        if int(df.Sec[i]+2) in round(gdf_R[(gdf_R.Year==df.Year[i])&(gdf_R.Month == df.Month[i]) & (gdf_R.Day == df.Day[i]) & (gdf_R.Hour == df.Hour[i]) & (gdf_R.Min == df.Min[i])].Sec).values:
            inRemotec.append(True)
        else:
            inRemotec.append(False) 
        date.append(datetime.date(df.Year[i],df.Month[i],df.Day[i]))    
    df["Date"] = date
    df["in_Remotec_data"] = inRemotec

    df.to_pickle(datapath + "/ACOS/DataFrames/DF6_"+RegionName+str(Num)+".pkl")
    #create GeoDataFrame
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Long, df.Lat))
    gdf.crs = {'init' :'epsg:4326'}

    if Num >= 900:
        Transcom = pd.read_pickle("/home/eschoema/Transcom_Regions.pkl")
        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]

    gdf.to_pickle(datapath + "/ACOS/DataFrames/GDF6_"+RegionName+str(Num)+".pkl")


