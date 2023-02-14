#!/usr/bin/env python

#Skript to create geopandas dataframe with transcom rgeions
#Created on October 20th 2020

#@author: E.-M. Schoemann


import numpy as np
import xarray as xr
#import rasterio
import pandas as pd
import geopandas
import glob
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import contextily as ctx
import mplleaflet
import geoplot as gplt
import geoplot.crs as gcrs
import argparse
import h5py
from netCDF4 import Dataset

savepath = "."

# Asia 8
df = pd.DataFrame(data = {'transcom':'Eurasia Temperate'},index=[0])

#create dataframe with border between eurasia temperate and tropical asia
lat_coord = [0  ,  34 , 36, 40, 43, 43.5, 43.5, 44.5, 46, 47, 47, 51, 56.5, 56.5, 52, 53, 53, 50, 50, 51, 51, 50, 48, 50, 47, 41, 49, 51, 16.44, 28.19, 27.1, 27.1, 22.73, 24.83, 26.34, 25.03, 21.95, 20.98, 23.03, 21.57, 18.95, 16.44  , 0]
long_coord= [-10, -10, 25, 26, 29, 40  , 45  , 44  , 53, 54, 59, 62, 64.5, 67  , 71, 73, 81, 86, 99, 99,102,102,122,125,128,132,149 ,179, 179,  122.59, 119.27, 111.34, 101.69, 97.61, 95.1, 92.02, 90.47, 86.84, 82.52, 79.67, 79.92, 81.95, 81.95]

geometry = [Polygon(zip(long_coord,lat_coord))]
gdf = geopandas.GeoDataFrame(df, crs="epsg:4326",geometry=geometry)
# dataset containing all continents and their countries
world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

ETs = world[(world.continent=='Asia')|(world.continent=='Europe')] #select only south america
ETs = ETs[['continent','geometry']]
ETs.insert(loc=1,column='dis',value = np.ones(len(ETs.continent)))
ETs = ETs.dissolve(by = 'dis').reset_index() #make one polygon out off all countries

ET = geopandas.overlay(ETs, gdf, how='intersection') #take intersection of the two dataframes
ET = ET[['continent','geometry']]
ET.insert(loc=2,column='transcom',value= 'ET')



# north Asia 9
df = pd.DataFrame(data = {'transcom':'Eurasia Boreal'},index=[0])

#create dataframe with border between eurasia temperate and tropical asia
lat_coord = [89, 89, 51, 56.5, 56.5, 52, 53, 53, 50, 50, 51, 51, 50, 48, 50, 47, 41, 49, 51]
long_coord= [179,62, 62, 64.5, 67  , 71, 73, 81, 86, 99, 99,102,102,122,125,128,132,149,179]

geometry = [Polygon(zip(long_coord,lat_coord))]
gdf = geopandas.GeoDataFrame(df, crs="epsg:4326",geometry=geometry)
# dataset containing all continents and their countries
world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

EBs = world[(world.continent=='Asia')|(world.continent=='Europe')] #select only south america
EBs = EBs[['continent','geometry']]
EBs.insert(loc=1,column='dis',value = np.ones(len(EBs.continent)))
EBs = EBs.dissolve(by = 'dis').reset_index() #make one polygon out off all countries

EB = geopandas.overlay(EBs, gdf, how='intersection') #take intersection of the two dataframes
EB = EB[['continent','geometry']]
EB.insert(loc=2,column='transcom',value= 'EB')


# Europe Transcom
df = pd.DataFrame(data = {'transcom':'Europe'},index=[0])

#create dataframe with border between eurasia temperate and tropical asia
lat_coord = [34 , 36, 40, 43, 43.5, 43.5, 44.5, 46, 47, 47, 51, 72 , 72]
long_coord= [-11, 25, 26, 29, 40  , 45  , 44  , 53, 54, 59, 62, 62, -11]

geometry = [Polygon(zip(long_coord,lat_coord))]
gdf = geopandas.GeoDataFrame(df, crs="epsg:4326",geometry=geometry)
# dataset containing all continents and their countries
world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

Eurs = world[(world.continent=='Asia')|(world.continent=='Europe')] #select only south america
Eurs = Eurs[['continent','geometry']]
Eurs.insert(loc=1,column='dis',value = np.ones(len(Eurs.continent)))
Eurs = Eurs.dissolve(by = 'dis').reset_index() #make one polygon out off all countries

Eur = geopandas.overlay(Eurs, gdf, how='intersection') #take intersection of the two dataframes
Eur = Eur[['continent','geometry']]
Eur.insert(loc=2,column='transcom',value= 'Eur')

#South America
df = pd.DataFrame(data = {'transcom':'Tropical Asia'},index=[0])

#create dataframe with border between eurasia temperate and tropical asia
lat_coord = [16.44, 18.95, 21.57, 23.03, 20.98, 21.95, 25.03, 26.34, 24.83, 22.73, 27.10, 27.10, 28.19, -19.12, -19.12, 16.44]
long_coord= [81.95, 79.92, 79.67, 82.52, 86.84, 90.47, 92.02, 95.10, 97.61,101.69,111.34,119.27,122.59, 172.82, 100,    81.95]

geometry = [Polygon(zip(long_coord,lat_coord))]
gdf = geopandas.GeoDataFrame(df, crs="epsg:4326",geometry=geometry)
# dataset containing all continents and their countries
world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

TAs = world[(world.continent=='Asia')|(world.name=='Papua New Guinea')] #select only south america
TAs = TAs[['continent','geometry']]
TAs.insert(loc=1,column='dis',value = np.ones(len(TAs.continent)))
TAs = TAs.dissolve(by = 'dis').reset_index() #make one polygon out off all countries

TA = geopandas.overlay(TAs, gdf, how='intersection') #take intersection of the two dataframes
TA = TA[['continent','geometry']]
TA.insert(loc=2,column='transcom',value= 'TA')


#South America
df = pd.DataFrame(data = {'transcom':'southern american temperate'},index=[0])

#create dataframe with border between southern american tropical and south american temperate in the north and a rectangle for the east, west and south border
lat_coord = [-7.7  ,-7.7  ,-11  ,-11  ,-9   ,-9   ,-3   ,-2   ,-2  ,-58,-58,-7.7]
long_coord= [-80.5 ,-77   ,-67  ,-63.5,-61  ,-53  ,-46.5,-39.5,-30 ,-30,-84,-84 ]

geometry = [Polygon(zip(long_coord,lat_coord))]
gdf = geopandas.GeoDataFrame(df, crs="epsg:4326",geometry=geometry)
# dataset containing all continents and their countries
world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

SA = world[world.continent == 'South America'] #select only south america
SA = SA[['continent','geometry']]
SA = SA.dissolve(by = 'continent').reset_index() #make one polygon out off all countries

SAT = geopandas.overlay(SA, gdf, how='intersection') #take intersection of the two dataframes
SAT = SAT[['continent','geometry']]
SAT.insert(loc=2,column='transcom',value= 'SAT')

#South America Tropical southwards the equator
df = pd.DataFrame(data = {'transcom':'southern american tropical'},index=[0])

#create dataframe with border between southern american tropical and south american temperate in the north and a rectangle for the east, west and south border
lat_coord = [-7.7  ,-7.7  ,-11  ,-11  ,-9   ,-9   ,-3   ,-2   ,-2  ,0,0,-7.7]
long_coord= [-80.5 ,-77   ,-67  ,-63.5,-61  ,-53  ,-46.5,-39.5,-30 ,-30,-84,-84 ]

geometry = [Polygon(zip(long_coord,lat_coord))]
gdf = geopandas.GeoDataFrame(df, crs="epsg:4326", geometry=geometry)

# dataset containing all continents and their countries
world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

SArSO = world[world.continent == 'South America'] #select only south america
SArSO = SArSO[['continent','geometry']]
SArSO = SArSO.dissolve(by = 'continent').reset_index() #make one polygon out off all countries

SATrSO = geopandas.overlay(SArSO, gdf, how='intersection') #take intersection of the two dataframes
SATrSO = SATrSO[['continent','geometry']]
SATrSO.insert(loc=2,column='transcom',value= 'SATr_SE')


#South America Tropical
df = pd.DataFrame(data = {'transcom':'southern american tropical'},index=[0])

#create dataframe with border between southern american tropical and south american temperate in the north and a rectangle for the east, west and south border
lat_coord = [-7.7  ,-7.7  ,-11  ,-11  ,-9   ,-9   ,-3   ,-2   ,-2  ,15,15,-7.7]
long_coord= [-80.5 ,-77   ,-67  ,-63.5,-61  ,-53  ,-46.5,-39.5,-30 ,-30,-84,-84 ]

geometry = [Polygon(zip(long_coord,lat_coord))]
gdf = geopandas.GeoDataFrame(df, crs="epsg:4326", geometry=geometry)

# dataset containing all continents and their countries
world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

SAr = world[world.continent == 'South America'] #select only south america
SAr = SAr[['continent','geometry']]
SAr = SAr.dissolve(by = 'continent').reset_index() #make one polygon out off all countries

SATr = geopandas.overlay(SAr, gdf, how='intersection') #take intersection of the two dataframes
SATr = SATr[['continent','geometry']]
SATr.insert(loc=2,column='transcom',value= 'SATr')

#whole South America
SAwhole = world[world.continent == 'South America'] #select only south america
SAwhole = SAwhole[['continent','geometry']]
SAwhole = SAwhole.dissolve(by = 'continent').reset_index()

SAwhole.insert(loc=2,column='transcom',value= 'SAT_SATr')

#Africa
df2 = pd.DataFrame(data = {'transcom':'southern Africa'},index=[0])

lat_coord2 = [0  ,0  ,-37  ,-37  ,0]
long_coord2= [7  ,62 ,62   ,7    ,7]

geometry2 = [Polygon(zip(long_coord2,lat_coord2))]
gdf2 = geopandas.GeoDataFrame(df2, crs="epsg:4326", geometry=geometry2)


world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

Af = world[world.continent == 'Africa']
Af = Af[['continent','geometry']]
Af = Af.dissolve(by = 'continent').reset_index()

SAf = geopandas.overlay(Af, gdf2, how='intersection')
SAf = SAf[['continent','geometry']]
SAf.insert(loc=2,column='transcom',value= 'SA')

#whole Africa
Afwhole = world[world.continent == 'Africa'] #select only south america
Afwhole = Afwhole[['continent','geometry']]
Afwhole = Afwhole.dissolve(by = 'continent').reset_index()

Afwhole.insert(loc=2,column='transcom',value= 'Af')


#Australia
AU = world[(world.name == 'New Zealand')|(world.name == 'Australia')]
AU = AU[['continent','geometry']]
AU = AU.dissolve(by = 'continent').reset_index()

AU.insert(loc=2,column='transcom',value= 'AU')


#Indonesia and Papua New Giunea
IP = world[(world.name == 'Indonesia')|(world.name == 'Papua New Guinea')]
IP = IP[['continent','geometry']]

IP.insert(loc=1,column='dis',value = ['1','1'])
IP = IP.dissolve(by = 'dis').reset_index()
IP = IP.drop('dis', 1)

IP.insert(loc=2,column='transcom',value= 'I_P')


#Europe with Russia
Eu = world[(world.continent == 'Europe')]
Eu = Eu[['continent','geometry']]
Eu = Eu.dissolve(by = 'continent').reset_index()
Eu.insert(loc=2,column='transcom',value= 'EU')


#Asia without Russia ad without I_P
As = world[(world.continent == 'Asia')]
As = As[['continent','geometry']]
As = As.dissolve(by = 'continent').reset_index()
As = geopandas.overlay(As, TA, how='difference')
As = As[['continent','geometry']]
As.insert(loc=2,column='transcom',value= 'As')

#North America South
df3 = pd.DataFrame(data = {'transcom':'southern America'},index=[0])

lat_coord2 = [49  ,49  ,0,0]
long_coord2= [-160  ,-60 ,-60 ,-160]

geometry3 = [Polygon(zip(long_coord2,lat_coord2))]
gdf3 = geopandas.GeoDataFrame(df3, crs="epsg:4326", geometry=geometry3)


NAS = world[(world.continent == 'North America')&~(world.name == 'Canada')&~(world.name == 'Greenland')]
NAS = NAS[['continent','geometry']]
NAS = NAS.dissolve(by = 'continent').reset_index()
NAS = geopandas.overlay(NAS, gdf3, how='intersection')
NAS = NAS[['continent','geometry']]
NAS.insert(loc=2,column='transcom',value= 'NAS')

#North America North
NAN = geopandas.overlay(world[(world.continent == 'North America')],NAS,how='difference')
NAN = NAN[['continent','geometry']]
NAN = NAN.dissolve(by = 'continent').reset_index()
NAN.insert(loc=2,column='transcom',value= 'NAN')

Transcom = AU.copy()
Transcom = Transcom.append(SAT)
Transcom = Transcom.append(SAwhole)
Transcom = Transcom.append(SATr)
Transcom = Transcom.append(SATrSO)
Transcom = Transcom.append(SAf)
Transcom = Transcom.append(Afwhole)
Transcom = Transcom.append(IP)
Transcom = Transcom.append(TA)
Transcom = Transcom.append(Eu)
Transcom = Transcom.append(As)
Transcom = Transcom.append(NAN)
Transcom = Transcom.append(NAS)
Transcom = Transcom.append(Eur)
Transcom = Transcom.append(ET)
Transcom = Transcom.append(EB)
Transcom.to_pickle(savepath + "/Transcom_Regions.pkl")


