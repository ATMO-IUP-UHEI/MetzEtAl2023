#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 17:38:52 2021
Create pandas dataframes out of TRENDY flux data
@author: eschoema
"""

import datetime
#from datetime import timedelta
import numpy as np
import pandas as pd
from functions import getAreaOfGrid, getGdfOfGrid 
import glob, os
import geopandas
import xarray as xr
import math
from shapely.geometry import Polygon
from skimage.measure import block_reduce
from RegionParam import getRegion
from functools import reduce
import cftime
from netCDF4 import Dataset

datapath = "."
savepath = "."

#https://gis.stackexchange.com/questions/326408/how-aggregate-data-in-a-geodataframe-by-the-geometry-in-a-geoseries

def CreateDF(DS,Datatype,ModelInfo,ncTime):
    timeVar = ModelInfo.tVar.values[0]
       
    # get the index of the first timevalue after 2008/1/1 and before 2008/1/2
    #cftime can not deal with month calendars
    if ModelInfo.tformat.values[0] == 'm':
        itstart = int(math.ceil((2008-ModelInfo.starttime.values[0].year)*12
                      +1-ModelInfo.starttime.values[0].month
                      -DS.variables[timeVar].values[0]))
    else: 
        itstart = cftime.date2index(datetime.datetime(2008,1,1),
                                    ncTime,
                                    calendar = ncTime.calendar,
                                    select = 'after')
        #check whether itstart before 2008/2/1
        if itstart == cftime.date2index(datetime.datetime(2008,2,1),
                                        ncTime,
                                        calendar = ncTime.calendar,
                                        select = 'after'):
            print('Error in determining startindex at 2008/1/1')
    print(itstart)
    print(len(DS.variables[Datatype]))
    print('start reading the data')
    for t in range(itstart,len(DS.variables[Datatype])):
        if len(DS.variables['lat'+ModelInfo.LatName.values[0]].values) == 180: 
            #the data is on a 1°x1° grid    
            dfFlux = pd.DataFrame(
                            data=DS.variables[Datatype][t].values,#data=DS.variables[Datatype][t][:][0].values,
                            index = DS.variables['lat'+ModelInfo.LatName.values[0]].values, 
                            columns = DS.variables['lon'+ModelInfo.LonName.values[0]].values, 
                            dtype='float')
        elif len(DS.variables['lat'+ModelInfo.LatName.values[0]].values) == 360: 
            #the data is on a 0.5°x0.5° grid 
            dfFlux = pd.DataFrame(
                            data=block_reduce(DS.variables[Datatype][t].values,#block_reduce(DS.variables[Datatype][t][:][0].values,
                                              block_size=(2,2),
                                              func=np.nanmean),
                            index = block_reduce(
                                       DS.variables['lat'+ModelInfo.LatName.values[0]].values, 
                                       block_size=((2,)),
                                       func=np.mean),
                            columns = block_reduce(
                                       DS.variables['lon'+ModelInfo.LonName.values[0]].values,
                                       block_size=((2,)),
                                       func=np.mean),
                            dtype='float' )
        else:
            #the regridding is done when creating the geodataframe   
            dfFlux = pd.DataFrame(
                            data=DS.variables[Datatype][t].values,#data=DS.variables[Datatype][t][0].values,
                            index = DS.variables['lat'+ModelInfo.LatName.values[0]].values, 
                            columns = DS.variables['lon'+ModelInfo.LonName.values[0]].values, 
                            dtype='float')
        #Unpivot dataframe (df)     
        dfFlux = pd.melt(dfFlux.reset_index(), id_vars='index',
                                      value_name =Datatype,var_name = 'Long')
        dfFlux["Lat"] = dfFlux['index']
        #### Latitude problem in CABLE-POP
        # The CABLE-POP caontains latitudes from -85.3 to 93.5 
        # Maps of the data show that there is no shift in the latitude. 
        # The coordinates are right and unphysical coordinates >90 get deleted
        dfFlux = dfFlux.drop(dfFlux[(dfFlux.Lat >90)].index,axis=0).reset_index()
        #get date
        if ModelInfo.tformat.values[0] == 'm':
            date = datetime.date(
                     int(ModelInfo.starttime.values[0].year + math.floor(DS.variables[timeVar].values[t]/12)),
                     int(ModelInfo.starttime.values[0].month +DS.variables[timeVar].values[t]%12),
                     15)
        else:
            cfdate = cftime.num2date(ncTime[t].data,
                                   ncTime.units,
                                   ncTime.calendar)
            #the 'only_use_cftime_datetimes=False' option to cobert directely 
            #into datetime.date does not wirk for noleap calendar
            date = datetime.date(cfdate.year,cfdate.month,cfdate.day)
        year = date.year
        month = date.month
        day = date.day
        dfFlux.insert(loc=1,column='Year',
                      value= np.ones(len(dfFlux["Lat"]))*int(year))        
        dfFlux.insert(loc=1,column='Month',
                      value= np.ones(len(dfFlux["Lat"]))*int(month))
        dfFlux.insert(loc=1,column='Day',
                      value= np.ones(len(dfFlux["Lat"]))*int(day))
        dfFlux.insert(loc=1,column='Date',
                      value= [date]*len(dfFlux["Lat"]))
        dfFlux.insert(loc=1,column='MonthDate',
                      value= [datetime.date(year,month,15)]*len(dfFlux["Lat"]))
        if t == itstart:
            print('the first time step i on '+str(year)+str(month).zfill(2)+str(day).zfill(2))
            df = dfFlux.copy()
        else:
            df = df.append(dfFlux)
                     
    return df

def CreateDataFrameTRENDYflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Num,
                            VegModel,ModelInfo):
    #geodataframe (gdf) with borders of the Transcom regions
    Transcom = pd.read_pickle(savepath + "/Transcom_Regions.pkl")
    #gdf with latitude dependent area of a 1x1° cell
    Area = getAreaOfGrid()
        
    datapath = datapath + "/TRENDY/"+VegModel+"/"  
    for num, Datatype in enumerate(ModelInfo.Variables.values[0].split('-')):
        try:
            gdfGrid = pd.read_pickle(savepath + "/Grid1_1"+RegionName+".pkl")
        except:
            gdfGrid = getGdfOfGrid(Num)
        filepath = datapath + VegModel + "_S3_" + Datatype + ".nc"
        print(filepath)
        #decode time false as 365day/noleap calender can not be decoded
        DS = xr.open_mfdataset(filepath,combine='by_coords',concat_dim='None'
                               ,decode_times=False)
        print(DS.keys())
        print(DS)
        #read time variable as netCDF4._netCDF4.Variable to use cftime module
        ncTime = Dataset(filepath).variables[ModelInfo.tVar.values[0]]
        #check that ncTime and DS.time are in the same order
        if not (ncTime== DS[ModelInfo.tVar.values[0]].values).all():
            print('Error in time varaible comparison of netcdf and xarray')
        #create Dataframe
        df3= CreateDF(DS,Datatype,ModelInfo,ncTime)
        
        df = df3[(df3.Long >= Long_min)
         & (df3.Long <= Long_max)
         & (df3.Lat >= Lat_min)
         & (df3.Lat <= Lat_max)]
       
        print("finished reading data")
        df.to_pickle(datapath + "/TRENDY/dataframes/DF1"+
                     ModelInfo.Model.values[0]+Datatype.split('_')[0]+"_"+RegionName+str(Num)+".pkl")
        
        #create GeoDataFrame
        gdf = geopandas.GeoDataFrame(
            df, geometry=geopandas.points_from_xy(df.Long, df.Lat),crs=4326)

        #interpolate grid on 1x1° grid     
        if len(DS.variables['lat'+ModelInfo.LatName.values[0]].values) not in [180,360]:
            print('Interpolate on 1°x1° grid')
            #get grid cell size
            dLong = (DS.variables['lon'+ModelInfo.LonName.values[0]].values[2]
                    -DS.variables['lon'+ModelInfo.LonName.values[0]].values[1])
            dLat = (DS.variables['lat'+ModelInfo.LatName.values[0]].values[2]
                    -DS.variables['lat'+ModelInfo.LatName.values[0]].values[1])
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
            # back to degrees
            inter = inter.to_crs(crs=4326)
            #get indices of all not nan values 
            nonanindex = ~np.isnan(inter[Datatype])
            #calculate weighted mean on the 1x1° gris (using GridID)
            newdf = inter.loc[inter.index[nonanindex]].groupby(
                        ['GridID','Month','Year','Day'])[Datatype].agg(
                        lambda x: np.average(
                        x, weights=inter.loc[x.index, "subarea"])).reset_index()
            # change geometry column in reference grid to prepare for merging
            gdfGrid = gdfGrid.drop(columns='geometry')
            gdfGrid = gdfGrid.rename(columns={'geomPoint':'geometry'})
            gdfGrid = gdfGrid.set_geometry(gdfGrid.geometry)
            del gdf
            #get  Lat, Lon and geometry in the new gdf 'newdf'->'gdf'
            gdf = pd.merge(gdfGrid[['GridID','geometry','Lat','Long']],newdf,on=['GridID'])
            gdf = gdf.drop(columns = 'GridID')

            #insert date columns
            date = gdf.apply(lambda x: datetime.date(int(x.Year),
                                                     int(x.Month),
                                                     int(x.Day)),axis=1)
            gdf.insert(loc=1,column='Date',value=date)
            Monthdate = gdf.apply(lambda x: datetime.date(int(x.Year),
                                                          int(x.Month),
                                                          15),axis=1)
            gdf.insert(loc=1,column='MonthDate',value=Monthdate)
        else:
            if Num >= 900: 
                #not needed for the interpolated gdf as gdfGrid already cut by Transcom border
                Transcom = pd.read_pickle(savepath + "/Transcom_Regions.pkl")
                igdf = gdf.within(Transcom[(Transcom.transcom == RegionName)].geometry.iloc[0])
                gdf = gdf.loc[igdf]
        if Num == 950: 
            SemiArid = pd.read_pickle(datapath + "/ERA5/Australia/DataFrames/MaskDryPoly0.02_4.pkl")
            iSemiAridgdf = gdf.within(SemiArid.geometry.iloc[0])
            gdf = gdf.loc[iSemiAridgdf]
        elif Num == 951: 
            Wet = pd.read_pickle(datapath + "/ERA5/Australia/DataFrames/MaskWetPoly0.02_4.pkl")
            iWetgdf = gdf.within(Wet.geometry.iloc[0])
            gdf = gdf.loc[iWetgdf]
        elif Num == 952: 
            SemiAridLC = pd.read_pickle(datapath + "/ERA5/Australia/DataFrames/MaskDryPoly_LC2014_Alstroem.pkl")
            iSemiAridgdfLC = gdf.within(SemiAridLC.geometry.iloc[0])
            gdf = gdf.loc[iSemiAridgdfLC]
        elif Num == 953: 
            WetLC = pd.read_pickle(datapath + "/ERA5/Australia/DataFrames/MaskWetPoly_LC2014_Alstroem.pkl")
            iWetgdfLC = gdf.within(WetLC.geometry.iloc[0])
            gdf = gdf.loc[iWetgdfLC]
        
        
        print("calculate total fluxes") 
        gdf = pd.merge(gdf,Area,on=['Lat'],how='left')
        gdf.insert(loc=1,column=Datatype.split('_')[0] + 'tot',
                      value= gdf[Datatype]*gdf['Area']) 
        gdf = gdf.rename(columns={Datatype:Datatype.split('_')[0]})
        
        gdf.to_pickle(datapath + "/TRENDY/dataframes/GDF1"
                      +ModelInfo.Model.values[0]+Datatype.split('_')[0]+"_"+RegionName+str(Num)+".pkl")
        
        del gdf,df,DS

TrendyModels = ['CABLE-POP', 'CLASSIC', 'DLEM', 'IBIS', 'ISAM', 'JSBACH', 
            'LPX-Bern', 'OCN', 'ORCHIDEE', 'ORCHIDEE-CNP', 'ORCHIDEEv3', 
            'LPJ-GUESS','CLM5.0','ISBA-CTRIP','JULES-ES-1p0','LPJ','SDGVM',
            'VISIT','YIBs']

tCalendar = ['standard','365_day','noleap','','noleap','proleptic_gregorian',
        'noleap','365_day','noleap','noleap','noleap',
        '','noleap','gregorian','365_day','365_day','360_day',
        '','']

LatName = ['itude','itude','','itude','','itude',
             'itude','itude','itude','','',
             'itude','','_FULL','itude','','itude',
             '','itude']

LonName = ['gitude','gitude','','gitude','','gitude',
             'gitude','gitude','gitude','','',
             'gitude','','_FULL','gitude','','gitude',
             '','gitude']

ModelVars = ['gpp-npp-ra-rh',
             'gpp-nbp-npp-ra-rh-fFire',
             'gpp-npp-ra-rh',
             'gpp-nbp-npp-ra-rh',
             'gpp-nbp-npp-ra-rh',
             'gpp-nbp-npp-ra-rh-fFire-fLuc',#'rh', #
             'gpp-nbp-npp-ra-rh-fFire',
             'gpp-nbp-npp-ra-rh',
             'gpp-nbp-npp-ra-rh-fFire',
             'gpp-nbp-npp-ra-rh-fFire',
             'gpp-nbp-npp-ra-rh',
             'gpp-nbp-npp-ra-fFire',
             'gpp-nbp-npp-ra-rh-fFire',
             'gpp-nbp-npp-ra-rh-fFire',
             'gpp-nbp-rh-npp_nlim',
             'gpp-nbp-npp-ra-rh-fFire',
             'gpp-nbp-npp-ra-rh-fFire',
             'gpp-nbp-ra-rh-fFire',
             'gpp-nbp-npp-ra-rh']
#msl contains multiple depth. When treating this varialbe in the func CreateDF muss ein [0] hinter [t]
#ModelVars = ['msl','msl','msl','msl',
#            'msl','msl','msl','msl',
#            'msl','msl','msl','msl',
#            'msl','msl','msl','msl',
#            'msl','msl','msl']

#the time format is only relevant for models with mothly values as for those the cftime function can not be applied
tformat = ['s','d','m','m','d','h',
           'd','d','d','d','d',
           'd','d','m','s','d','h',
           'm','m']

tVar = ['time','time','time','time','time','time',
        'time','time','time','time','time_counter',
        'time','time','time_counter','time','time','time',
        'time','time']

tstart = [datetime.date(1700,1,1),
          datetime.date(1700,12,31),
          datetime.date(1700,1,1),
          datetime.date(1700,1,1),
          datetime.date(1700,1,1),
          datetime.date(1700,1,16),
          datetime.date(1700,1,1),
          datetime.date(1700,1,1),
          datetime.date(1700,1,1),
          datetime.date(1701,1,1),
          datetime.date(1700,1,1),
          datetime.date(1700,1,1),
          datetime.date(1700,1,1),
          datetime.date(1700,1,16),
          datetime.date(2010,1,1),
          datetime.date(1700,1,1),
          datetime.date(1900,1,1),
          datetime.date(1860,1,1),
          datetime.date(1700,1,1)]

dfModelInfo = pd.DataFrame({'Model':TrendyModels,
                            'tcalendar':tCalendar,
                          'Variables': ModelVars,
                          'tformat':tformat,
                          'tVar': tVar,
                          'starttime':tstart,
                          'LatName':LatName,
                          'LonName':LonName} )

Numm = 901
RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)
    
    
for model in ['JSBACH']:#['IBIS', 'ISAM', 'JSBACH' 
            #'LPX-Bern', 'OCN', 'ORCHIDEE', 'ORCHIDEE-CNP', 'ORCHIDEEv3', 
            #'CLM5.0','JULES-ES-1p0','SDGVM',
            #'VISIT','CABLE-POP', 'CLASSIC', ]:#YiBs, LPJ-GUESS, LPJ,'DLEM',,'ISBA-CTRIP'
    ModelInfo = dfModelInfo[(dfModelInfo.Model == model)]
    CreateDataFrameTRENDYflux(Lat_min,Lat_max,Long_min,Long_max,RegionName,Numm,
                            model,ModelInfo)

