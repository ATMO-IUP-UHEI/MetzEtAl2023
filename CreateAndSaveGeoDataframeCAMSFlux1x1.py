#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# Create a GeoDataFrame out of CAMS 1x1 CO2 flux data
#Author: E.-M. Schoemann
#Date 26.11.2020

import datetime
#from datetime import timedelta
import read_remotec_out
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb
from functions import createPandasDateFrame, getGdfOfGrid 
import glob, os
import geopandas
import xarray as xr
import time
import math
from skimage.measure import block_reduce
from RegionParam import getRegion
from shapely.geometry import Polygon

datapath = "."
savepath = "."

def CreateDF(DS2,DS,Datatype,year,month):
     
    dArea = pd.DataFrame(data=DS2.variables['area'].values,index = DS.variables['latitude'].values,
                                  columns =  DS.variables['longitude'].values, dtype='float' )
    dArea = pd.melt(dArea.reset_index(), id_vars='index',
                                 value_name ='area',var_name = 'Long')
    dArea["Lat"] = dArea['index']
    
    dArea2 = pd.DataFrame(data=DS2.variables['lsf'].values,index = DS.variables['latitude'].values,
                                  columns =  DS.variables['longitude'].values, dtype='float' )
    dArea2 = pd.melt(dArea2.reset_index(), id_vars='index',
                                 value_name ='lsf',var_name = 'Long')
    dArea2["Lat"] = dArea2['index']
    dArea = pd.merge(dArea,dArea2,on=['Lat','Long'],how='left')
    #print(dfy2)
        
    for d in [0]:#range(len(DS['flux_apos_bio'])):
        for nam in ['flux_apos_bio','flux_apri_bio','flux_apos_ocean','flux_apri_ocean','flux_foss']:        
            
            dfyCO2 = pd.DataFrame(data=DS.variables[nam].values*1000*(44/12),index = DS.variables['latitude'].values, columns = DS.variables['longitude'].values, dtype='float' )
            dfyCO2 = pd.melt(dfyCO2.reset_index(), id_vars='index',
                                      value_name =nam+'0',var_name = 'Long')
            dfyCO2["Lat"] = dfyCO2['index']
            dfyCO2 = pd.merge(dfyCO2, dArea, on=['Lat','Long'], how = 'left')
            dfyCO2.insert(loc=1,column='Year',value= np.ones(len(dfyCO2["Lat"]))*year)        
            dfyCO2.insert(loc=1,column='Month',value= np.ones(len(dfyCO2["Lat"]))*month)
            
            if nam in 'flux_apos_bio':
                dfy = dfyCO2.copy()
            else:
                dfy = pd.merge(dfyCO2, dfy, on=['Lat','Long'])        
                dfy = dfy.drop(columns=['Month_y','Year_y'])
                dfy = dfy.rename(columns={"Month_x": "Month","Year_x":"Year"})
            del dfyCO2
        if d == 0:
            df = dfy.copy()
        else:
            df = df.append(dfy)
                     
    return df

def CreateDataFrameCAMSflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num,Datatype=''):
    if Datatype in 'Sat':
        starty = 2014
        paths = datapath + "/CAMS/Flux/cams73_latest_co2_flux_satellite_mm_"
        DS2 = xr.open_mfdataset(datapath + "/CAMS/Flux/cams73_latest_co2_flux_satellite_mm_201905.nc",combine='by_coords',concat_dim='None',decode_times=False)
    
    else:
        starty = 2009
        paths = datapath + "/CAMS/Flux/cams73_latest_co2_flux_surface_mm_"
        DS2 = xr.open_mfdataset(datapath + "/CAMS/Flux/cams73_latest_co2_flux_surface_mm_201905.nc",combine='by_coords',concat_dim='None',decode_times=False)
    
    for y in range(starty,2020): #not os.path.isfile(datapath + "/GFAS/dataframes/DF1_"+RegionName+str(Num)+".pkl"):
        print("Start reading data for : " +str(y))
        dp = paths+str(y)+"*.nc"
        for num, filepath in enumerate(glob.glob(dp)):
            #print(filepath)
            DS = xr.open_mfdataset(filepath, combine='by_coords',concat_dim='None',drop_variables = 'time_components')#decode_times=False)
            mon = int(filepath[-5:-3])
            #create Dataframe
            df3= CreateDF(DS2,DS,Datatype,y,mon)
            del DS
            df2 = df3[(df3.Long >= Long_min)
             & (df3.Long <= Long_max)
             & (df3.Lat >= Lat_min)
             & (df3.Lat <= Lat_max)]

            if num == 0:
                df = df2.copy()
            else:           
                df = df.append(df2, ignore_index=True)
        
        print("finished reading data") 

        #create date variable
        print("create timestamp")
        date = df.apply(lambda x: datetime.date(int(x.Year),
                                                  int(x.Month),
                                                  15),axis=1)
        df.insert(loc = 1, value = date, column = 'Date') 
        df.to_pickle(datapath + "/CAMS/dataframes/DF1FLUX"+Datatype+str(y)+"_"+RegionName+str(Num)+".pkl")
        df.drop(columns = ['index_x_x','index_y_x','index_y_y','index_x_y','lsf_x','lsf_y'],inplace=True)
        #create GeoDataFrame
        gdf = geopandas.GeoDataFrame(
            df, geometry=geopandas.points_from_xy(df.Long, df.Lat),crs = 4326)

        print('Interpolate on 1°x1° grid')
        #get reference 1x1° grid geodataframe of the region
        try:
            gdfGrid = pd.read_pickle(savepath + "/Grid1_1"+RegionName+".pkl")
        except:
            gdfGrid = getGdfOfGrid(Num)
        #get grid cell size
        dLong = (DS2.variables['longitude'].values[2]
                -DS2.variables['longitude'].values[1])
        dLat = (DS2.variables['latitude'].values[2]
                -DS2.variables['latitude'].values[1])
        #create polygon geomatry for geodataframe
        geom = gdf.apply(lambda x: Polygon(zip(
                [x.Long-dLong/2,x.Long-dLong/2,x.Long+dLong/2,x.Long+dLong/2],
                [x.Lat-dLat/2,x.Lat+dLat/2,x.Lat+dLat/2,x.Lat-dLat/2])),axis = 1)
        gdf.insert(loc = 1, value = geom, column = 'geomPoly')
        gdf = gdf.set_geometry(gdf.geomPoly)
        #cut the gdf grid with the 1x1° grid and keep all subgridcells of the intersection
        inter = geopandas.overlay(gdfGrid,gdf,how='intersection')
        # to m² insteadt of degree to get area of the subgridcells
        inter = inter.to_crs({'proj':'cea'}) 
        inter.insert(loc = 1, value = inter.area,column = 'subarea')
        inter = inter.to_crs(crs=4326) # back to degrees
        firstVar = True
        for DataVar in ['flux_apos_bio0','flux_apri_bio0','flux_apos_ocean0','flux_apri_ocean0','flux_foss0']:
            nonanindex = ~np.isnan(inter[DataVar])
            #calculate weighted mean on the 1x1° grid (using GridID)
            newdf2 = inter.loc[inter.index[nonanindex]].groupby(
                        ['GridID','Month','Year'])[DataVar].agg(
                        lambda x: np.average(
                        x, weights=inter.loc[x.index, "subarea"])).reset_index()
            if firstVar:
                newdf = newdf2.copy()
            else:
                newdf = newdf.merge(newdf2,on=['GridID','Month','Year'],how='outer')
            firstVar = False
        gdfGrid = gdfGrid.drop(columns='geometry')
        gdfGrid = gdfGrid.rename(columns={'geomPoint':'geometry'})
        gdfGrid = gdfGrid.set_geometry(gdfGrid.geometry)
        del gdf
        #get  Lat, Lon and geometry in the new gdf 'newdf'
        gdf = pd.merge(gdfGrid[['GridID','geometry','Lat','Long','AreaGrid']],newdf,on=['GridID'])
        gdf = gdf.drop(columns = 'GridID')
        for DataVar in ['flux_apos_bio','flux_apri_bio','flux_apos_ocean','flux_apri_ocean','flux_foss']:
            gdf.insert(loc = 1, value = gdf[DataVar+'0']*gdf['AreaGrid'],column = DataVar+'_tot')
        date = gdf.apply(lambda x: datetime.date(int(x.Year),
                                                     int(x.Month),
                                                     int(15)),axis=1)
        gdf.insert(loc=1,column='Date',value=date)
        Monthdate = gdf.apply(lambda x: datetime.date(int(x.Year),
                                                        int(x.Month),
                                                        15),axis=1)
        gdf.insert(loc=1,column='MonthDate',value=Monthdate)
        #this should change nothing as the 1x1° gdf is already on Australian region
        if Num >= 900:
            if y == starty:
                Transcom = pd.read_pickle(savepath + "/Transcom_Regions.pkl")
            igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
            gdf = gdf.loc[igdf]

        gdf.to_pickle(datapath + "/CAMS/dataframes/GDF1x1FLUX"+Datatype+str(y)+"_"+RegionName+str(Num)+".pkl")

        del gdf,df
        
def CombineGDFCAMS(RegionName, Num,Datatype):
    if 'Sat' in Datatype:
        starty = 2014
    else:
        starty = 2009
    for y in range(starty,2020):
        df2 = pd.read_pickle(datapath + "/CAMS/dataframes/GDF1x1FLUX"+Datatype+str(y)+"_"+RegionName+str(Num)+".pkl")
        if y == starty:
            df = df2.copy()
        else:           
            df = df.append(df2, ignore_index=True)
            
    df.to_pickle(datapath + "/CAMS/dataframes/GDF1x1FLUX"+Datatype+"_"+RegionName+str(Num)+".pkl")
'''
#main
Numm  = 949
RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)

CreateDataFrameCAMSflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,'Sur')   
CombineGDFCAMS(RegionName, Numm,'Sur') 
'''