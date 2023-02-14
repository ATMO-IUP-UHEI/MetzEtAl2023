#!/usr/bin/env python
# regridd the 0.1x0.1° GFAS data to 1°x1° resolution
# Author: E.-M. Schoemann


import numpy as np
import pandas as pd
import sys
import datetime
import geopandas
import xarray as xr
import h5py 
from netCDF4 import Dataset
import read_remotec_out
import time
import calendar
import argparse
from skimage.measure import block_reduce

datapath = "."

#tutorial on https://unidata.github.io/netcdf4-python/netCDF4/index.html#section6

# Argument Parser
def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to regrid CAMS GFAS data from 0.1°x0.1° to 1°x1°")
    parser.add_argument("date", type=str, help="enter date as string in form 'yyyymm'")
    args = parser.parse_args()
    return args


#settings:
args = parse_arguments()
date = args.date
inputfile = datapath + "/GFAS/cams_gfas_daily_0_1x0_1_"+date+".nc"


#read out file
DS = xr.open_mfdataset(inputfile, combine='by_coords',concat_dim='None',decode_times=False)

#settings for the first file
#newday = True
filename = datapath + "/GFAS/Daily_1x1_regridded/cams_gfas_daily_regrid1x1_"+date #output file


    
rootgrp = Dataset(filename+".nc", "w", format="NETCDF4_CLASSIC")

# set dimensions
time = rootgrp.createDimension("time", None) #unlimited size
latitude = rootgrp.createDimension("latitude", 180)
longitude = rootgrp.createDimension("longitude", 360)

# set varaibles
times = rootgrp.createVariable(varname = "time", datatype = "i4", dimensions=("time"))
times.units = "hours since 1900-01-01 00:00:00.0"
times.calendas = "gregorian"
latitude = rootgrp.createVariable(varname = "latitude", datatype = "f4", dimensions=("latitude"))
latitude.units = "degrees_north"
latitude.comments = "Center latitude of the 1°x1° grid cell"
latitude.long_name = "latitude"
longitude = rootgrp.createVariable(varname = "longitude", datatype = "f4", dimensions=("longitude"))
longitude.units = "degrees_east"
longitude.comments = "Center longitude of the 1°x1° grid cell"
longitude.long_name = "longitude"
co2fire = rootgrp.createVariable(varname = "co2fire", datatype = "f4", dimensions=("time","latitude","longitude"),fill_value = -32767)
co2fire.units = "kg m**-2 s**-1"
co2fire.scale_factor = 1.0276405974730713E-9
co2fire.add_offset = 3.367167181680266E-5
co2fire.missing_value = -32767
co2fire.long_name = "Wildfire flux of Carbon Dioxide"
cofire = rootgrp.createVariable(varname = "cofire", datatype = "f4", dimensions=("time","latitude","longitude"),fill_value = -32767)
cofire.units = "kg m**-2 s**-1"
cofire.scale_factor = 6.006043256737052E-11
cofire.add_offset = 1.9679401335024625E-6
cofire.missing_value = -32767
cofire.long_name = "Wildfire flux of Carbon Monoxide"
cfire = rootgrp.createVariable(varname = "cfire", datatype = "f4", dimensions=("time","latitude","longitude"),fill_value = -32767)
cfire.units = "kg m**-2 s**-1"
cfire.scale_factor = 3.0006472231709516E-10
cfire.add_offset = 9.83192069144194E-6
cfire.missing_value = -32767
cfire.long_name = "Wildfire overall flux of burnt Carbon"
ch4fire = rootgrp.createVariable(varname = "ch4fire", datatype = "f4", dimensions=("time","latitude","longitude"),fill_value = -32767)
ch4fire.units = "kg m**-2 s**-1"
ch4fire.scale_factor = 2.726054636698209E-12
ch4fire.add_offset = 8.93219062260535E-8
ch4fire.missing_value = -32767
ch4fire.long_name = "Wildfire flux of Methane"
mami = rootgrp.createVariable(varname = "mami", datatype = "f4", dimensions=("time","latitude","longitude"),fill_value = -32767)
mami.units = "m"
mami.scale_factor = 0.10616025513863245
mami.add_offset = 3478.4469198724305
mami.missing_value = -32767
mami.long_name = "Mean altitude of maximum injection"
apt = rootgrp.createVariable(varname = "apt", datatype = "f4", dimensions=("time","latitude","longitude"),fill_value = -32767)
apt.units = "m"
apt.scale_factor = 0.15259487586406847
apt.add_offset = 4999.923702562068
apt.missing_value = -32767
apt.long_name = "Altitude of plume top"
# ---------------------------------------------------------------------------------------

# empty the parameter lists
#aapb= []
aapt = []
acfire = []
ach4fire = []
aco2fire = []
acofire= []
#ainjh = []
amami = []


times[:] = DS.time.values
longitude[:] = np.array(range(5,3600,10))/10
latitude[:] = np.array(range(1795,-5,-10))/10-90
# for each timestep
for i in range(len(DS.time)):
    print(DS.time.values[i])
    print(datetime.datetime.now())

    #aapb.append((block_reduce(DS.apb.values[i],block_size=(10,10),func=np.mean)).tolist())
    aapt.append((block_reduce(DS.apt.values[i],block_size=(10,10),func=np.mean)).tolist())
    acfire.append((block_reduce(DS.cfire.values[i],block_size=(10,10),func=np.mean)).tolist())
    acofire.append((block_reduce(DS.cofire.values[i],block_size=(10,10),func=np.mean)).tolist())
    aco2fire.append((block_reduce(DS.co2fire.values[i],block_size=(10,10),func=np.mean)).tolist())
    ach4fire.append((block_reduce(DS.ch4fire.values[i],block_size=(10,10),func=np.mean)).tolist())
    #ainjh.append((block_reduce(DS.injh.values[i],block_size=(10,10),func=np.mean)).tolist())
    amami.append((block_reduce(DS.mami.values[i],block_size=(10,10),func=np.mean)).tolist())

    
#write in variables
apt[:] = aapt
cfire[:] = acfire
cofire[:] = acofire
co2fire[:] = aco2fire
ch4fire[:] = ach4fire
mami[:] = amami


