#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 15.12.2022

Skript to create geopandas dataframe based on my own transcom region definition
out of the monthly 1x1° datasets of MIP

@author: eschoema
"""

import datetime
#from datetime import timedelta
import numpy as np
import pandas as pd
from functions import getAreaOfGrid
import geopandas
import xarray as xr
import time
import math
from RegionParam import getRegion

datapath = "."
savepath = "."

def CreateDF(DS):
    for d in range(len(DS['land'])): #over time
           
        dfy = pd.DataFrame(data=DS.variables['land'][d].values,
                                index = DS.variables['latitude'].values, 
                                columns = DS.variables['longitude'].values, 
                                dtype='float' )
        dfy = pd.melt(dfy.reset_index(), id_vars='index',
                                    value_name ='land',var_name = 'Long')
        dfy["Lat"] = dfy['index']
        dfy.insert(loc=1,column='Year',value= np.ones(len(dfy["Lat"]))*int(DS.start_date[d][0]))        
        dfy.insert(loc=1,column='Month',value= np.ones(len(dfy["Lat"]))*int(DS.start_date[d][1]))
        dfy.insert(loc=1,column='Day',value= np.ones(len(dfy["Lat"]))*15)
        
            
        if d == 0:
            df = dfy.copy()
        else:
            df = df.append(dfy)
        del dfy         
    return df

def CreateDataFrameMIPflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num,version):
    print("Start reading data")
    if version == 'regular':
        filepath = datapath + "/OCO-2_v10_MIP/LNLGIS/EnsMean_gridded_fluxes_LNLGIS.nc4"
        filepathStd = datapath + "/OCO-2_v10_MIP/LNLGIS/EnsStd_gridded_fluxes_LNLGIS.nc4"
    elif version == 'comp':
        filepath = datapath + "/OCO-2_v10_MIP/MIP_v10_Sourishnew/gridded_fluxes_original_20221106/flux_mip/gridded_fluxes/EnsMean_gridded_fluxes_LNLGIS.nc4"
        filepathStd = datapath + "/OCO-2_v10_MIP/MIP_v10_Sourishnew/gridded_fluxes_original_20221106/flux_mip/gridded_fluxes/EnsStd_gridded_fluxes_LNLGIS.nc4"
    elif version == 'regularIS':
        filepath = datapath + "/OCO-2_v10_MIP/IS/EnsMean_gridded_fluxes_IS.nc4"
        filepathStd = datapath + "/OCO-2_v10_MIP/IS/EnsStd_gridded_fluxes_IS.nc4"
    else:
        print('version does not exist')

    DSMean = xr.open_dataset(filepath)
    DSStd = xr.open_dataset(filepathStd)
   
    #create Dataframe
    df3Mean= CreateDF(DSMean)
    df3Std= CreateDF(DSStd)
    df3Mean.drop(columns='index',inplace = True)
    df3Std.drop(columns='index',inplace = True)
    dfMean = df3Mean[(df3Mean.Long >= Long_min) & (df3Mean.Long <= Long_max) & (df3Mean.Lat >= Lat_min) & (df3Mean.Lat <= Lat_max)]
    dfStd = df3Std[(df3Std.Long >= Long_min) & (df3Std.Long <= Long_max) & (df3Std.Lat >= Lat_min) & (df3Std.Lat <= Lat_max)]
    dfStd.rename(columns = {'land':'landStd'},inplace=True)
    df = pd.merge(dfMean,dfStd,on =['Year','Month','Day','Lat','Long'])


    print("create timestamp")
    date = df.apply(lambda x: datetime.date(int(x.Year),
                                                     int(x.Month),
                                                     int(x.Day)),axis=1)
    df.insert(loc=1,column='Date',value=date)
    df.insert(loc=1,column='MonthDate',value=date)

    #gdf with latitude dependent area of a 1x1° cell
    Area = getAreaOfGrid()
    df = pd.merge(df, Area, on=['Lat'])     
    # former unit gC/(m^2*year)  -> TgC/(gridcell*month) by multiplying with * m²/gridcell * 1/12 *1//10**12
    landt = df.apply(lambda x: (x.land * 1/12 * x.Area/10**12),axis=1)
    print('works')
    landStdt = df.apply(lambda x: (x.landStd * x.Area/10**12),axis=1)
    
    df.insert(loc=1, column = 'Landtot', value = landt)
    df.insert(loc=1, column = 'LandStdtot', value = landStdt)
    
    if version == 'regular':
        df.to_pickle(datapath + "/OCO-2_v10_MIP/dataframes/DFMIP_"+RegionName+str(Num)+".pkl")
    elif version == 'regularIS':
        df.to_pickle(datapath + "/OCO-2_v10_MIP/dataframes/DFMIP_IS_"+RegionName+str(Num)+".pkl")
    else:
        df.to_pickle(datapath + "/OCO-2_v10_MIP/MIP_v10_Sourishnew/dataframes/DFMIP_"+RegionName+str(Num)+".pkl")
    
    #create GeoDataFrame
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(df.Long, df.Lat),crs = 4326)
 
    if Num >= 900:
        Transcom = pd.read_pickle(savepath + "/Transcom_Regions.pkl")
        igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
        gdf = gdf.loc[igdf]
   
    if version == 'regular':
        gdf.to_pickle(datapath + "/OCO-2_v10_MIP/dataframes/GDFMIP_"+RegionName+str(Num)+".pkl")
    elif version == 'regularIS':
        gdf.to_pickle(datapath + "/OCO-2_v10_MIP/dataframes/GDFMIP_IS_"+RegionName+str(Num)+".pkl")
    else:
        df.to_pickle(datapath + "/OCO-2_v10_MIP/MIP_v10_Sourishnew/dataframes/GDFMIP_"+RegionName+str(Num)+".pkl")
    
    
    del gdf,df
    
#main
if __name__=='__main__':
    Numm  = 949
    RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)

    CreateDataFrameMIPflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,'regularIS')    
