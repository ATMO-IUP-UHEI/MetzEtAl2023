#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 16:20:05 2021

@author: eschoema
"""
# Script to create dataframe for a month containing TM5-4DVAR XCO2 values created from 
# the 3 hourly mixingratios combined with pressurelevels

datapath = "."
savepath = "."

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
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Script create monthly mean values for unsampled tm5 data")
    parser.add_argument("startdate", type=str, help="enter startdate as string in form 'yyyy-mm'")
    parser.add_argument("enddate", type=str, help="enter enddate as string in form 'yyyy-mm'")    
    parser.add_argument("Num", type=int, help="Region Number")
    parser.add_argument("Dtype", type=str, help="Enter Assimilation type")
    
    args = parser.parse_args()
    return args


def CreateDF(DS,DS2,yearv,monthv,dayv):
    bt = DS.variables['bt'].values
    at = DS.variables['at'].values
    # for all the 8 time steps per file
    for t in range(len(DS2.variables['mix'].values[0])):
        surfp = DS2.variables['pressure'].values[t]
        # for all 25 levels per timestep
        for l in range(len(DS2.variables['mix'].values[0][t])):
            if l == 0:
                la = np.array(DS2.variables['mix'].values[0][t][l]*((surfp*(bt[l]-bt[l+1]))+(np.ones((90,120))*(at[l]-at[l+1]))))
                la = la[np.newaxis,...]
            else:
                la2 = np.array(DS2.variables['mix'].values[0][t][l]*((surfp*(bt[l]-bt[l+1]))+(np.ones((90,120))*(at[l]-at[l+1]))))
                la2 = la2[np.newaxis,...]
                la = np.append(la,la2, axis = 0)
        # now la contains for every level a 2D array with mixing ratio * dPressure dim = [level, lat, long]
        co2a = np.sum(la,axis = 0)/((surfp*(bt[0]-bt[-1]))+(np.ones((90,120))*(at[0]-at[-1])))
        # now co2a contains a 90x120 array with the XCO2 values for this timestep 
        if t == 0:
            co2 = co2a[np.newaxis,...]
        else:
            co2 = np.append(co2,co2a[np.newaxis,...], axis = 0)
        # co2 now contains a 90x120 array with the XCO2 values for each timestep
    co2mean = np.nanmean(co2, axis = 0)
    # co2mean now contains a 90x120 array with the daily mean XCO2 values for one day           
    dCO2 = pd.DataFrame(data=co2mean,index = DS2.variables['latitude'].values,
                    columns =  DS2.variables['longitude'].values, dtype='float' )
    dCO2 = pd.melt(dCO2.reset_index(), id_vars='index',
                   value_name ='CO2',var_name = 'Long')
    dCO2["Lat"] = dCO2['index']

    dCO2.insert(loc=1,column='Year',value= np.ones(len(dCO2["Lat"]))*yearv) 
    dCO2.insert(loc=1,column='Month',value= np.ones(len(dCO2["Lat"]))*monthv)
    dCO2.insert(loc=1,column='Day',value= np.ones(len(dCO2["Lat"]))*dayv)
    dCO2.insert(loc=1,column='Date',value= np.repeat(datetime.date(yearv,monthv,dayv),len(dCO2["Lat"])))

    return dCO2

def CreateDataFrameTM5noncsY(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num, yearv, monthv,Dtype):
    savepath = datapath + "/TK5_4DVAR/dataframes/Unsampled/"
    print("Start reading data:")
    for num, filepath in enumerate(glob.glob(datapath + "/TK5_4DVAR/Molefractions/"+Dtype+"/output/2009010100-2019070100/mix/"+str(yearv)+"/"+str(monthv).zfill(2)+"/*.nc4")):  
        dayv = int(filepath[-6:-4])
        DS = xr.open_mfdataset(filepath,combine='by_coords',concat_dim='None')
        DS2 = xr.open_mfdataset(filepath,group='glb300x200',combine='by_coords',concat_dim='None')
        if num == 0:        
        #create Dataframe
            df3= CreateDF(DS,DS2,yearv,monthv,dayv)
            df = df3[(df3.Long >= Long_min)
                    & (df3.Long <= Long_max)
                    & (df3.Lat >= Lat_min)
                    & (df3.Lat <= Lat_max)]
        else:      
        #create Dataframe
            df3= CreateDF(DS,DS2,yearv,monthv,dayv)
            df2 = df3[(df3.Long >= Long_min)
                    & (df3.Long <= Long_max)
                    & (df3.Lat >= Lat_min)
                    & (df3.Lat <= Lat_max)]
            df = df.append(df2, ignore_index=True)
    
    df.drop_duplicates(keep = 'first')   
    df = df.groupby(['Year','Month','Lat','Long'])['CO2'].mean().reset_index()
    df.insert(loc=1,column='Date',value= np.repeat(datetime.date(yearv,monthv,15),len(df["Lat"])))
    
    print("finished reading data") 
   
    df.to_pickle(savepath+"DF_Monthlynoncs1_"+Dtype+"_"+RegionName+str(Num)+"_"+str(yearv)+str(monthv).zfill(2)+".pkl")
    #create GeoDataFrame
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Long, df.Lat))
    gdf.crs = {'init' :'epsg:4326'}
 
    if Num >= 900:
        Transcom = pd.read_pickle(savepath + "/Transcom_Regions.pkl")
        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]
   
    gdf.to_pickle(savepath+"GDF_Monthlynoncs1_"+Dtype+"_"+RegionName+str(Num)+"_"+str(yearv)+str(monthv).zfill(2)+".pkl")
    objects = dir()
    for obj in objects:
      if not obj.startswith("__") and not obj in [DS,DS2,Lat_min,Lat_max,Long_min,Long_max,RegionName,Num, yearv, monthv, Dtype]:
        del globals()[obj]
    
  
def CreateDataFrameTM5noncs(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm, year_min,year_max, month_min,month_max, Dtype):
    savepath = datapath + "/TK5_4DVAR/dataframes/Unsampled/"
    
    for y in range(year_min,year_max+1,1):
        mstart = 1
        mend = 12
        print(y)
        if y == year_min:
            mstart = month_min
        if y == year_max:
            mend == month_max
        for m in range(mstart,mend+1):
            print(m)
            CreateDataFrameTM5noncsY(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm, y, m, Dtype)
            
    print("Merging dataFrames:")
    #for num, file in enumerate(glob.glob("CoSampledCAMS_surf*.pkl")):
    for num, file in enumerate(glob.glob(savepath+"GDF_Monthlynoncs1_"+Dtype+"_"+RegionName+str(Numm)+"_*.pkl")):
        print("read: "+file)
        data = pd.read_pickle(file)
        if num == 0:
            #create Dataframe
            df= data
        else:
            #add data to the Dataframe
            df = df.append(data, ignore_index=True)  #append one dataframe to another dataframe
    
    dfm = df.groupby(['Date','Year','Month'])['CO2'].mean().reset_index()
    dfm.insert(loc=1,column='MonthDate',value=dfm.Date)
    
    df.to_pickle(savepath+"GDF_MonthlynoncsGrid_"+Dtype+"_"+RegionName+str(Numm)+"_merged.pkl")
    dfm.to_pickle(savepath+"GDF_MonthlynoncsMean_"+Dtype+"_"+RegionName+str(Numm)+"_merged.pkl")
    
#main
# SETTINGS:
if __name__=='__main__':
    args = parse_arguments()
    year_min = int(args.startdate[0:4])#2009
    month_min = int(args.startdate[5:7])#4
    year_max = int(args.enddate[0:4])#2019
    month_max = int(args.enddate[5:7])#12
    Numm = args.Num
    Dtype = args.Dtype
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)
    CreateDataFrameTM5noncs(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm, year_min,year_max, month_min,month_max, Dtype)
