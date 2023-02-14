#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-
# Script to plot TRANSCOM Map and FLUXCOM and ObsPack Stations
# author: E.-M. Schoemann, 7.1.2022

import pandas as pd
import geopandas
import xarray as xr
import os
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import matplotlib as mpl

datapath = "."
savepath = "."

font = {'family' : 'Arial'}
mpl.rc('font', **font)

#TCCON
dfTCCON = pd.DataFrame(data={'Name':['Wollogong','Darwin'],
                              'Lat':[-34.406,-12.4559],
                              'Lon':[150.8798,130.9266]})
gdfTCCON = geopandas.GeoDataFrame(dfTCCON,
                                  geometry = geopandas.points_from_xy(dfTCCON.Lon,dfTCCON.Lat),
                                  crs = 4326)

#OzFlux
dfOzFlux = pd.DataFrame(data={'Name':['ASM','DalyUncleared','DryRiver'],
                              'Lat':[-22.2828,-14.1592,-15.2588],
                              'Lon':[133.2493,131.3881,132.3706]})
gdfOzFlux = geopandas.GeoDataFrame(dfOzFlux,
                                  geometry = geopandas.points_from_xy(dfOzFlux.Lon,dfOzFlux.Lat),
                                  crs = 4326)

# create Geodataframe with polygons and aggregate TRANSCOM Regions
if not os.path.isfile(datapath + "/CT2019/TranscomRegions.pkl"):
    filepath = datapath + "/CT2019/regions.nc"
    Transcom = xr.open_mfdataset(filepath, combine='by_coords',concat_dim='None',drop_variables='time_components')#,decode_times=False)
    dfTranscom = pd.DataFrame(data=Transcom.variables['transcom_regions'].values,
                                index = Transcom.variables['latitude'].values,
                                columns =  Transcom.variables['longitude'].values, dtype='float' )
    dfTranscom = pd.melt(dfTranscom.reset_index(), id_vars='index',
                                     value_name ='transcom_regions',var_name = 'Long')
    dfTranscom["Lat"] = dfTranscom['index']
    print('creating Polygons')
    polyList = dfTranscom.apply(lambda x: Polygon(zip([x.Long-0.5,x.Long-0.5,x.Long+0.5,x.Long+0.5],[x.Lat-0.5,x.Lat+0.5,x.Lat+0.5,x.Lat-0.5])),axis=1)
    dfTranscom.insert(loc=1,column='geomPoly',value=polyList)
    print('creating gdf')
    gdfTranscom = geopandas.GeoDataFrame(dfTranscom, crs = 'epsg:4326',geometry =dfTranscom.geomPoly)             
    print('dissolving')
    gdfTranscomRegions = gdfTranscom.dissolve(by='transcom_regions') 
    gdfTranscomRegions.to_pickle(datapath + "/CT2019/TranscomRegions.pkl")
else:
    gdfTranscomRegions = pd.read_pickle(datapath + "/CT2019/TranscomRegions.pkl")

# get Australian Transcom Region
TranscomPoly = pd.read_pickle(savepath + "/Transcom_Regions.pkl")
TranscomAU = TranscomPoly[(TranscomPoly.transcom == 'AU')]
world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
# additional POC station
dfPOC = pd.DataFrame(data={'Name':['POC'],
                              'Lat':[-20],
                              'Lon':[175]})
gdfPOC = geopandas.GeoDataFrame(dfPOC,
                                  geometry = geopandas.points_from_xy(dfPOC.Lon,dfPOC.Lat),
                                  crs = 4326)

# get ObsPack Station data
dfObsPack = pd.read_pickle(datapath + "/obspack_co2_1_CARBONTRACKER_CT2019B_2020-11-03/data/nc/SationsList.pkl")
gdfObsPack = geopandas.GeoDataFrame(dfObsPack
                ,geometry = geopandas.points_from_xy(dfObsPack.Longitude, dfObsPack.Latitude)
                ,crs = 4326)
kindObs = gdfObsPack.apply(lambda x: x.Type.split('-')[0], axis = 1)
gdfObsPack.insert(loc = 1, value = kindObs,column = 'MajorType')
gdfObsPack = gdfObsPack[(gdfObsPack.YrStart <=2018)&(gdfObsPack.YrEnd >=2009)]
gdfAirCraft = gdfObsPack[(gdfObsPack.MajorType == 'aircore')|(gdfObsPack.MajorType == 'aircraft')]
gdfShip = gdfObsPack[(gdfObsPack.MajorType == 'shipboard')]
gdfSurface = gdfObsPack[(gdfObsPack.MajorType == 'surface')|(gdfObsPack.MajorType == 'tower')]

#get FLUXCOM stations
dfFLUXCOM = pd.read_csv(datapath + "/FLUXCOM/Tramonta2016FCStations.txt"
                        ,delimiter = ' '
                        ,header= 0)
gdfFLUXCOM = geopandas.GeoDataFrame(dfFLUXCOM
                ,geometry = geopandas.points_from_xy(dfFLUXCOM.Long, dfFLUXCOM.Lat)
                ,crs = 4326)

MaskDry = pd.read_pickle(datapath + "/ERA5/Australia/DataFrames/MaskDry0.02_4.pkl")
MaskWet = pd.read_pickle(datapath + "/ERA5/Australia/DataFrames/MaskWet0.02_4.pkl")

#get CarboScope Stations
CS06 = pd.read_csv(datapath + "/CarboScope/coords.s06v2021.txt", delimiter = ' ',skipinitialspace=True,skiprows = 0)
CS10 = pd.read_csv(datapath + "/CarboScope/coords.s10v2021.txt", delimiter = ' ',skipinitialspace=True,skiprows = 0)
CSEXT = pd.read_csv(datapath + "/CarboScope/coords.sEXTv2021.txt", delimiter = ' ',skipinitialspace=True,skiprows = 0)
gdfCS06 = geopandas.GeoDataFrame(CS06
                ,geometry = geopandas.points_from_xy(CS06.Lon, CS06.Lat)
                ,crs = 4326)
gdfCS10 = geopandas.GeoDataFrame(CS10
                ,geometry = geopandas.points_from_xy(CS10.Lon, CS10.Lat)
                ,crs = 4326)
gdfCSEXT = geopandas.GeoDataFrame(CSEXT
                ,geometry = geopandas.points_from_xy(CSEXT.Lon, CSEXT.Lat)
                ,crs = 4326)
print('CS06:'+str(len(gdfCS06[(gdfCS06.Lat>= -50)&(gdfCS06.Lat <= -5)&(gdfCS06.Lon>= 110)&(gdfCS06.Lon <= 180)])))
print('CS10:'+str(len(gdfCS10[(gdfCS10.Lat>= -50)&(gdfCS10.Lat <= -5)&(gdfCS10.Lon>= 110)&(gdfCS10.Lon <= 180)])))
print('CSEXT:'+str(len(gdfCSEXT[(gdfCSEXT.Lat>= -50)&(gdfCSEXT.Lat <= -5)&(gdfCSEXT.Lon>= 110)&(gdfCSEXT.Lon <= 180)])))

# figure S3
fig, ax = plt.subplots(figsize = (12,5))

world.plot(color = 'lightgrey',ax=ax)
TranscomAU.plot(color = 'grey',ax=ax)#'nipy_spectral')
gdfObsPack.plot(color = 'purple',  markersize= 120,marker = 'x',ax=ax, label = 'ObsPack surface and tower',zorder=2)
gdfPOC.plot(ax=ax,color = 'purple',label = 'ObsPack ship POC',zorder=2)
gdfFLUXCOM.plot(color = 'red', markersize= 120, marker = '+',ax=ax, label = 'FLUXNET tower',zorder=2)
MaskWet.boundary.plot(color = 'black',ax=ax,zorder=1)
MaskDry.plot(color = 'blue',ax=ax, label = 'semi-arid area',zorder=1)
gdfOzFlux.plot(ax=ax,color = 'red',zorder=2,label = 'OzFlux tower')
gdfTCCON.plot(ax=ax,color = 'orange',marker = "^",markersize= 100,zorder=2,label = 'TCCON stations')
for x, y, label in zip(gdfOzFlux.geometry.x, gdfOzFlux.geometry.y, gdfOzFlux.Name):
    if label == 'DryRiver':
        ax.annotate(label, xy=(x, y), xytext=(4, -10), textcoords="offset points",color = 'red',fontsize=12)
    else:
        ax.annotate(label, xy=(x, y), xytext=(4, 4), textcoords="offset points",color = 'red',fontsize=12)
for x, y, label in zip(gdfTCCON.geometry.x, gdfTCCON.geometry.y, gdfTCCON.Name):
    ax.annotate(label, xy=(x, y), xytext=(4, 4), textcoords="offset points",color = 'orange',fontsize=12)

ax.legend(loc = 3)
ax.set_xlim([110,180])
ax.set_ylim([-50,-10])
ax.set_ylabel('Latitude', fontsize = 12)
ax.set_xlabel('Longitude', fontsize = 12)
plt.savefig(savepath + "/Results/Plots/Paper2021/final/Supp_TranscomMapStations_SemiArid_AU7.png", dpi=300,bbox_inches='tight')
plt.savefig(savepath + "/Results/Plots/Paper2021/final/Supp_TranscomMapStations_SemiArid_AU7.pdf", dpi=300,bbox_inches='tight',format = 'pdf')
