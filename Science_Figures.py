#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 09:22:29 2023

@author: eschoema

Skript for Figures in Science Paper next to the skripts Plot_Carbon_Fluxes ans Plot_Timeseries_Detmers
"""

import pandas as pd
import geopandas 
import numpy as np
from RegionParam import getRegion
import matplotlib.pyplot as plt
import matplotlib as mpl
from functions import getMeanAmplitude

font = {'family' : 'Arial'}
mpl.rc('font', **font)

datapath = "."
savepath = "."

Numm = 949
RegionName, Long_min, Long_max, Lat_min, Lat_max = getRegion(Numm)


if True:
    #file 1
    FCGFED = pd.read_pickle(datapath + "/SIF/FLUXCOM/MonthMeans_FCGFED_AU949_newForSuppPaper.pkl")
    AllMMTrendyNBP = pd.read_pickle(datapath + "/TRENDY/dataframes/MonthFrames/AllMMTrendyNBPAU949_newforSuppPaper.pkl")   
    allSatModel =pd.read_pickle(datapath + "/TM5Inversion/dataframes_flux/MonthFrames/MonthMeansAllSatModels_newForSuppPaper_AU949.pkl")          
    allModel= pd.read_pickle(datapath + "/CT2019/dataframes/monthFrames/MonthMeansAllISModels_AU949_newForSuppPaper_.pkl")
    
    file1 = allSatModel[['Year','Month','TM5ACOSIS_nbp','TM5RTIS_nbp']]
    file1.insert(loc=2,column = 'TM5-4DVAR+RemoTeC/GOSAT',value = file1.TM5RTIS_nbp*12/44)
    file1.insert(loc=2,column = 'TM5-4DVAR+ACOS/GOSAT',value = file1.TM5ACOSIS_nbp*12/44)
    file1 = file1.drop(columns=['TM5ACOSIS_nbp','TM5RTIS_nbp'])
    file1 = pd.merge(file1, allModel[['Year','Month','TM5_nbp']], on=['Year','Month'])
    file1.insert(loc=2,column = 'TM5-4DVAR_in-situ',value = file1.TM5_nbp*12/44)
    file1 = file1.drop(columns=['TM5_nbp'])
    file1 = pd.merge(file1, AllMMTrendyNBP[['Year','Month','meanNBP','meanNBPgoodModel','meanNBPbadModel','meangppgoodModel','meangppbadModel','meanrhragoodModel','meanrhrabadModel']], on=['Year','Month'])
    file1.insert(loc=2,column = 'TRENDY_all',value = -file1.meanNBP*12/44)
    file1.insert(loc=2,column = 'TRENDY_selection',value = -file1.meanNBPgoodModel*12/44)
    file1.insert(loc=2,column = 'TRENDY_others',value = -file1.meanNBPbadModel*12/44)
    file1.insert(loc=2,column = 'TRENDY_selection_GPP',value = file1.meangppgoodModel*12/44)
    file1.insert(loc=2,column = 'TRENDY_others_GPP',value = file1.meangppbadModel*12/44)
    file1.insert(loc=2,column = 'TRENDY_selection_TER',value = file1.meanrhragoodModel*12/44)
    file1.insert(loc=2,column = 'TRENDY_others_TER',value = file1.meanrhrabadModel*12/44)
    file1 = file1.drop(columns=['meanNBP','meanNBPgoodModel','meanNBPbadModel','meangppgoodModel','meangppbadModel','meanrhragoodModel','meanrhrabadModel'])
    file1 = pd.merge(file1, FCGFED[['Year','Month','NEEtot']], on=['Year','Month'])
    file1.insert(loc=2,column = 'FLUXCOM_NEE',value = file1.NEEtot/10**12)
    file1 = file1.drop(columns=['NEEtot'])
    file1.to_csv(savepath + '/Results/Plots/Paper2021/final/Fluxes_TotalAustralia.csv')
    
    #file 2
    AllMMTrendyNBP950 = pd.read_pickle(datapath + "/TRENDY/dataframes/MonthFrames/AllMMTrendyNBPAU950_newforSuppPaper.pkl")   
    file2 = AllMMTrendyNBP950[['Year','Month','meanNBP','meanNBPgoodModel','meanNBPbadModel','meangppgoodModel','meangppbadModel','meanrhragoodModel','meanrhrabadModel','meanrhgoodModel','meanrabadModel','meanragoodModel','meanrhbadModel']]
    file2.insert(loc=2,column = 'TRENDY_all',value = -file2.meanNBP*12/44)
    file2.insert(loc=2,column = 'TRENDY_selection',value = -file2.meanNBPgoodModel*12/44)
    file2.insert(loc=2,column = 'TRENDY_others',value = -file2.meanNBPbadModel*12/44)
    file2.insert(loc=2,column = 'TRENDY_selection_GPP',value = file2.meangppgoodModel*12/44)
    file2.insert(loc=2,column = 'TRENDY_others_GPP',value = file2.meangppbadModel*12/44)
    file2.insert(loc=2,column = 'TRENDY_selection_TER',value = file2.meanrhragoodModel*12/44)
    file2.insert(loc=2,column = 'TRENDY_others_TER',value = file2.meanrhrabadModel*12/44)
    file2.insert(loc=2,column = 'TRENDY_selection_ra',value = file2.meanragoodModel*12/44)
    file2.insert(loc=2,column = 'TRENDY_others_ra',value = file2.meanrabadModel*12/44)
    file2.insert(loc=2,column = 'TRENDY_selection_rh',value = file2.meanrhgoodModel*12/44)
    file2.insert(loc=2,column = 'TRENDY_others_rh',value = file2.meanrhbadModel*12/44)
    file2 = file2.drop(columns=['meanNBP','meanNBPgoodModel','meanNBPbadModel','meangppgoodModel','meangppbadModel','meanrhragoodModel','meanrhrabadModel','meanrhgoodModel','meanrabadModel','meanragoodModel','meanrhbadModel'])
    file2.to_csv(savepath + '/Results/Plots/Paper2021/final/Fluxes_SemiAridAustralia.csv')
    
    #file 3
    AllMMTrendyNBP951 = pd.read_pickle(datapath + "/TRENDY/dataframes/MonthFrames/AllMMTrendyNBPAU951_newforSuppPaper.pkl")   
    file3 = AllMMTrendyNBP951[['Year','Month','meanNBP','meanNBPgoodModel','meanNBPbadModel','meangppgoodModel','meangppbadModel','meanrhragoodModel','meanrhrabadModel']]
    file3.insert(loc=2,column = 'TRENDY_selection',value = -file3.meanNBPgoodModel*12/44)
    file3.insert(loc=2,column = 'TRENDY_selection_GPP',value = file3.meangppgoodModel*12/44)
    file3.insert(loc=2,column = 'TRENDY_selection_TER',value = file3.meanrhragoodModel*12/44)
    file3 = file3.drop(columns=['meanNBP','meanNBPgoodModel','meanNBPbadModel','meangppgoodModel','meangppbadModel','meanrhragoodModel','meanrhrabadModel'])
    file3.to_csv(savepath + '/Results/Plots/Paper2021/final/Fluxes_NonSemiAridAustralia.csv')
    
    
    
    #file 4
    Month_means_TM5nocs =pd.read_pickle(datapath + "/TK5_4DVAR/dataframes/Unsampled/monthFrames/MonthMeans_TM5_newSuppPaper_IS_"+RegionName+ str(Numm)+".pkl")
    file4 = Month_means_TM5nocs[['Year','Month','CO2']]
    file4.insert(loc=2,column = 'TM5-4DVAR',value = file4.CO2)
    file4 = file4.drop(columns=['CO2'])
    file4.to_csv(savepath + '/Results/Plots/Paper2021/final/Concentrations_TotalAustralia.csv')
       
    
    
    
# calc IAV and annual fluxes
# dataframes created by Plot_CarbonFluxes.py
# make sure that thedatframes are created with 2009-01-01 as start date and not 2009-04-15!!!
if False:
    FCGFED = pd.read_pickle(datapath + "/SIF/FLUXCOM/MonthMeans_FCGFED_"+RegionName+str(Numm)+".pkl")
    AllMMTrendyNBP = pd.read_pickle(datapath + "/TRENDY/dataframes/MonthFrames/AllMMTrendyNBP"+RegionName+str(Numm)+".pkl")   
    allSatModel = pd.read_pickle(datapath + "/TM5Inversion/dataframes_flux/MonthFrames/MonthMeansAllSatModels_1x1_"+RegionName+str(Numm)+".pkl")          
    allModel = pd.read_pickle(datapath + "/CT2019/dataframes/monthFrames/MonthMeansAllISModels_"+RegionName+str(Numm)+".pkl")
    
    data = [FCGFED,AllMMTrendyNBP,AllMMTrendyNBP,AllMMTrendyNBP,allSatModel,allSatModel,allSatModel,allModel,allModel,allModel,allModel]
    var = ['nbp','meanNBPgoodModelSwitch','meanNBPbadModelSwitch','meanNBPSwitch','mean_nbp','TM5RTIS_nbp','TM5ACOSIS_nbp','mean_nbp','CAMSsur_nbp', 'CT_nbp', 'TM5_nbp']
    Yearvarl = ['_x','','','','','','','','','','']
    for i,model in enumerate(['FC+GFED','TRENDYgood','TRENDYbad','TRENDYall','allSat','RemoTeC','ACOS','allModel','CAMS','CT','TM5']):
        factor = 12/44
        dataset = data[i]
        varname = var[i]
        Yearvar = Yearvarl[i]
        YearsJJwo2009 = dataset[(dataset['JJYear'+Yearvar] >= 2010)&(dataset['JJYear'+Yearvar] < 2018)].groupby('JJYear'+Yearvar)[varname].sum().reset_index()#in TgCO2
        AnnualMeanJJwo2009 = YearsJJwo2009[varname].mean()*factor
        StDevJJwo2009 = YearsJJwo2009[varname].std(ddof=0)*factor
        YearsJJ = dataset[(dataset['JJYear'+Yearvar] >= 2009)&(dataset['JJYear'+Yearvar] < 2018)].groupby('JJYear'+Yearvar)[varname].sum().reset_index()#in TgCO2
        AnnualMeanJJ = YearsJJ[varname].mean()*factor
        StDevJJ = YearsJJ[varname].std(ddof=0)*factor
        Yearswo2009 = dataset[(dataset['Year'+Yearvar] > 2009)].groupby('Year'+Yearvar)[varname].sum().reset_index()#in TgCO2
        AnnualMeanwo2009 = Yearswo2009[varname].mean()*factor
        StDevwo2009 = Yearswo2009[varname].std(ddof=0)*factor
        Years = dataset.groupby('Year'+Yearvar)[varname].sum().reset_index()#in TgCO2
        AnnualMean = Years[varname].mean()*factor
        StDev = Years[varname].std(ddof=0)*factor
        
        YearsJJwo1011 = dataset[((dataset['JJYear'+Yearvar] > 2011)&(dataset['JJYear'+Yearvar] < 2018))|(dataset['JJYear'+Yearvar] == 2009)].groupby('JJYear'+Yearvar)[varname].sum().reset_index()#in TgCO2
        AnnualMeanJJwo1011 = YearsJJwo1011[varname].mean()*factor
        Yearswo200911 = dataset[(dataset['Year'+Yearvar] > 2011)].groupby('Year'+Yearvar)[varname].sum().reset_index()#in TgCO2
        AnnualMeanwo200911 = Yearswo200911[varname].mean()*factor
        
        
        dfIAVind = pd.DataFrame(data = {'Dataset':[model],
                                    'AnnualMean Jul-Jun without 2009':[AnnualMeanJJwo2009],
                                    'AnnualMean Jul-Jun':[AnnualMeanJJ],
                                    'AnnualMean without 2009':[AnnualMeanwo2009],
                                    'AnnualMean':[AnnualMean],
                                    'StDev Jul-Jun without 2009':[StDevJJwo2009],
                                    'StDev Jul-Jun':[StDevJJ],
                                    'StDev without 2009':[StDevwo2009],
                                    'StDev':[StDev],
                                    'AnnualMean Jul-Jun without 10/11':[AnnualMeanJJwo1011],
                                    'AnnualMean without 2009-2011':[AnnualMeanwo200911]})
        if i == 0:
            dfIAV = dfIAVind.copy(deep = True)
        else:
            dfIAV = dfIAV.append(dfIAVind,ignore_index=True)
    
        
        amplitude, StDevAmplitude, minimum, maximum = getMeanAmplitude(dataset,varname, 'JJYear'+Yearvar, 2009, 2017)

#Map Satellite measurements, S4
if False:
    gdfRemoTeC = pd.read_pickle(datapath + "/GOSAT_Markus/dataframes/GDF_"+RegionName+str(Numm)+".pkl")
    gdfRemoTeC = gdfRemoTeC[(gdfRemoTeC.meas_geom == '0')] #no change for AU949
    gdfRemoTeC = gdfRemoTeC[['Sec', 'Min', 'Hour', 'Day','Month', 'Year', 'Lat', 'Long','geometry']]
    gdfRemoTeC = gdfRemoTeC[(gdfRemoTeC.Year < 2019)]
    
    gdfACOS = pd.read_pickle(datapath + "/ACOS/DataFrames/GDF6_AU949.pkl")
    gdfACOS = gdfACOS[(gdfACOS.quality == 0)] #only measurements with good quality
    gdfACOS = gdfACOS[['Sec', 'Min', 'Hour', 'Day','Month', 'Year', 'Lat', 'Long','geometry']]
    gdfACOS = gdfACOS[(gdfACOS.Year < 2019)]
    
    #get 1x1 Grid with IDcolumn
    MaskDry = pd.read_pickle(datapath + "/ERA5/Australia/DataFrames/MaskDry0.02_4.pkl")
    MaskWet = pd.read_pickle(datapath + "/ERA5/Australia/DataFrames/MaskWet0.02_4.pkl")
    Grid = pd.concat([MaskDry,MaskWet])
    Grid.reset_index(inplace=True)
    Grid.drop(columns = ['aggr','monthlytp','index'],inplace=True)
    Grid.reset_index(inplace=True)
    Grid.rename(columns= {'index':'GridID'},inplace=True)
    
    #dissolve to get a 3x2 grid
    Grid2 = Grid.copy(deep = True)
    Grid2.insert(loc = 1, column = 'GridIDlat2',value =np.floor(Grid2.Lat/2))
    Grid2.insert(loc = 1, column = 'GridIDlong3',value =np.floor(Grid2.Long/3))
    Grid2 = Grid2.dissolve(by=['GridIDlat2','GridIDlong3']).reset_index()
    
    #get rectengular 1x1 grid
    Grid3 = pd.read_pickle(savepath + "/Grid1_1_-54_-1_10_179.pkl")
    Grid3 = Grid3.set_geometry('geomPoly')
    Grid3.insert(loc = 1, column = 'GridIDlat2',value =np.floor(Grid3.Lat/2))
    Grid3.insert(loc = 1, column = 'GridIDlong3',value =np.floor(Grid3.Long/3))
    Grid3 = Grid3.dissolve(by=['GridIDlat2','GridIDlong3']).reset_index()
    Grid4 = geopandas.overlay(Grid3,Grid2,how='intersection', keep_geom_type=True)
    Grid5 = geopandas.sjoin(left_df=Grid3, right_df=Grid2[['geometry']], how='left')
    Grid5 = Grid5.loc[Grid5['index_right'].notnull()]
    Grid5.drop(columns=['index', 'GridID','index_right'],inplace=True)
    Grid5 = Grid5.drop_duplicates(subset=['GridIDlat2','GridIDlong3'])
    Grid5 = Grid5.set_geometry('geomPoly')
    Grid5.insert(loc=1,column = 'GridID', value= Grid5.GridIDlat2*100+Grid5.GridIDlong3)
    
    gdfRemoTeC3x2 = geopandas.sjoin(gdfRemoTeC,Grid5[['GridID','geomPoly']],how='left')
    gdfACOS3x2 = geopandas.sjoin(gdfACOS,Grid5[['GridID','geomPoly']],how='left')
    
    gdfRemoTeC = geopandas.sjoin(gdfRemoTeC,Grid[['GridID','geometry']],how='left')
    gdfACOS = geopandas.sjoin(gdfACOS,Grid[['GridID','geometry']],how='left')
    
    
    TranscomPoly = pd.read_pickle(savepath + "/Transcom_Regions.pkl")
    AU = TranscomPoly[(TranscomPoly.transcom == 'AU')]
    
    #1x1
    gdfRemoTeC1x1wet = gdfRemoTeC[(gdfRemoTeC.Month >= 10)|(gdfRemoTeC.Month <= 3)]
    gdfRemoTeC1x1dry = gdfRemoTeC[(gdfRemoTeC.Month <= 9 )&(gdfRemoTeC.Month >= 4)]
    
    gdfACOS1x1wet = gdfACOS[(gdfACOS.Month >= 10)|(gdfACOS.Month <= 3)]
    gdfACOS1x1dry = gdfACOS[(gdfACOS.Month <= 9 )&(gdfACOS.Month >= 4)]
    
    #3x2
    gdfRemoTeC3x2wet = gdfRemoTeC3x2[(gdfRemoTeC3x2.Month >= 10)|(gdfRemoTeC3x2.Month <= 3)]
    gdfRemoTeC3x2dry = gdfRemoTeC3x2[(gdfRemoTeC3x2.Month <= 9 )&(gdfRemoTeC3x2.Month >= 4)]
    
    gdfACOS3x2wet = gdfACOS3x2[(gdfACOS3x2.Month >= 10)|(gdfACOS3x2.Month <= 3)]
    gdfACOS3x2dry = gdfACOS3x2[(gdfACOS3x2.Month <= 9 )&(gdfACOS3x2.Month >= 4)]
    
    
    
    res = '3x2'
    
    #map abundance in 1x1 grid
    if res == '1x1':
        gdfRemoTeCwet = gdfRemoTeC1x1wet.copy(deep=True)
        gdfRemoTeCdry = gdfRemoTeC1x1dry.copy(deep=True)
        gdfACOSwet = gdfACOS1x1wet.copy(deep=True)
        gdfACOSdry = gdfACOS1x1dry.copy(deep=True)
        
    else:
        gdfRemoTeCwet = gdfRemoTeC3x2wet.copy(deep=True)
        gdfRemoTeCdry = gdfRemoTeC3x2dry.copy(deep=True)
        gdfACOSwet = gdfACOS3x2wet.copy(deep=True)
        gdfACOSdry = gdfACOS3x2dry.copy(deep=True)
        Grid = Grid5.copy()
        
    gdfACOSwetNumber = gdfACOSwet.groupby('GridID')['Year'].count().reset_index()
    gdfACOSwetNumber.rename(columns= {'Year':'number'},inplace=True)
    gdfACOSwetNumber = pd.merge(Grid,gdfACOSwetNumber, on='GridID',how='outer')
    gdfACOSwetNumber['number'].fillna(0,inplace = True)
    
    gdfACOSdryNumber = gdfACOSdry.groupby('GridID')['Year'].count().reset_index()
    gdfACOSdryNumber.rename(columns= {'Year':'number'},inplace=True)
    gdfACOSdryNumber = pd.merge(Grid,gdfACOSdryNumber, on='GridID',how='outer')
    gdfACOSdryNumber['number'].fillna(0,inplace = True)
    
    gdfRemoTeCwetNumber = gdfRemoTeCwet.groupby('GridID')['Year'].count().reset_index()
    gdfRemoTeCwetNumber.rename(columns= {'Year':'number'},inplace=True)
    gdfRemoTeCwetNumber = pd.merge(Grid,gdfRemoTeCwetNumber, on='GridID',how='outer')
    gdfRemoTeCwetNumber['number'].fillna(0,inplace = True)
    
    gdfRemoTeCdryNumber = gdfRemoTeCdry.groupby('GridID')['Year'].count().reset_index()
    gdfRemoTeCdryNumber.rename(columns= {'Year':'number'},inplace=True)
    gdfRemoTeCdryNumber = pd.merge(Grid,gdfRemoTeCdryNumber, on='GridID',how='outer')
    gdfRemoTeCdryNumber['number'].fillna(0,inplace = True)
    #OrRd before
    
    #new
    factor =  12/44
    cm = 1/2.54  # centimeters in inches
    lw = 0.8
    cs = 1.3
    fig, ax = plt.subplots(2, 2, 
                            figsize = (18.0*cm,14*cm),
                            gridspec_kw={'width_ratios': [1, 1],'height_ratios':[1,1]})
    
    ax[0,0].text(0.06, 0.9, r'$\rm \bf{A}$', horizontalalignment='center', verticalalignment='center', transform=ax[0,0].transAxes,fontsize=10, weight='bold')
    Grid.plot(ax=ax[0,0], color = '#440154FF',linewidth = 0.5)
    gdfRemoTeCdryNumber.plot(ax=ax[0,0], column = 'number',cmap = 'viridis',legend = True,legend_kwds={'orientation': "horizontal",'shrink': 0.8},norm=mpl.colors.LogNorm(vmin=1, vmax=1000))
    Grid.boundary.plot(ax=ax[0,0], color = 'grey',linewidth = 0.5)
    AU.boundary.plot(ax=ax[0,0], color = 'red',linewidth = 1)
    ax[0,0].text(160,-13,'GOSAT/RemoTeC',fontsize=6)
    ax[0,0].text(160,-16,'April-September',fontsize=6)
    ax[0,0].tick_params(axis='both', which='major', labelsize=6)
    #ax[0,0].set_ylabel('Latitude',fontsize=6)
    figGP = ax[0,0].figure; cb_ax = figGP.axes[4]; cb_ax.tick_params(labelsize=6,length = 0.05)
    ax[0,0].text(0.5, -0.5, 'Number of measurements', horizontalalignment='center', verticalalignment='center', transform=ax[0,0].transAxes,fontsize=7)
    
    
    ax[0,1].text(0.06, 0.9, r'$\rm \bf{B}$', horizontalalignment='center', verticalalignment='center', transform=ax[0,1].transAxes,fontsize=10, weight='bold')
    Grid.plot(ax=ax[0,1], color = '#440154FF',linewidth = 0.5)
    gdfACOSdryNumber.plot(ax=ax[0,1], column = 'number',cmap = 'viridis',legend = True,legend_kwds={'orientation': "horizontal",'shrink': 0.8},norm=mpl.colors.LogNorm(vmin=1, vmax=1000))#gdfACOSdryNumber.number.min()
    Grid.boundary.plot(ax=ax[0,1], color = 'grey',linewidth = 0.5)
    AU.boundary.plot(ax=ax[0,1], color = 'red',linewidth = 1)
    ax[0,1].text(160,-13,'GOSAT/ACOS',fontsize=6)
    ax[0,1].text(160,-16,'April-September',fontsize=6)
    ax[0,1].tick_params(axis='both', which='major', labelsize=6)
    figGP = ax[0,1].figure; cb_ax = figGP.axes[5]; cb_ax.tick_params(labelsize=6,length = 0.05)
    ax[0,1].text(0.5, -0.5, 'Number of measurements', horizontalalignment='center', verticalalignment='center', transform=ax[0,1].transAxes,fontsize=7)
    
    ax[1,1].text(0.06, 0.9, r'$\rm \bf{D}$', horizontalalignment='center', verticalalignment='center', transform=ax[1,1].transAxes,fontsize=10, weight='bold')
    Grid.plot(ax=ax[1,1], color = '#440154FF',linewidth = 0.5)
    gdfACOSwetNumber.plot(ax=ax[1,1], column = 'number',cmap = 'viridis',legend = True,legend_kwds={'orientation': "horizontal",'shrink': 0.8},norm=mpl.colors.LogNorm(vmin=1, vmax=1000))
    Grid.boundary.plot(ax=ax[1,1], color = 'grey',linewidth = 0.5)
    AU.boundary.plot(ax=ax[1,1], color = 'red',linewidth = 1)
    ax[1,1].text(160,-13,'GOSAT/ACOS',fontsize=6)
    ax[1,1].text(160,-16,'October-March',fontsize=6)
    ax[1,1].tick_params(axis='both', which='major', labelsize=6)
    #ax[1,1].set_xlabel('Longitude',fontsize=6)
    figGP = ax[1,1].figure; cb_ax = figGP.axes[6]; cb_ax.tick_params(labelsize=6,length = 0.05)
    ax[1,1].text(0.5, -0.5, 'Number of measurements', horizontalalignment='center', verticalalignment='center', transform=ax[1,1].transAxes,fontsize=7)
    
    
    ax[1,0].text(0.06, 0.9, r'$\rm \bf{C}$', horizontalalignment='center', verticalalignment='center', transform=ax[1,0].transAxes,fontsize=10, weight='bold')
    Grid.plot(ax=ax[1,0], color = '#440154FF',linewidth = 0.5)
    gdfRemoTeCwetNumber.plot(ax=ax[1,0], column = 'number',cmap = 'viridis',legend = True,legend_kwds={'orientation': "horizontal",'shrink': 0.8},norm=mpl.colors.LogNorm(vmin=1, vmax=1000))
    Grid.boundary.plot(ax=ax[1,0], color = 'grey',linewidth = 0.5)
    AU.boundary.plot(ax=ax[1,0], color = 'red',linewidth = 1)
    ax[1,0].text(160,-13,'GOSAT/RemoTeC',fontsize=6)
    ax[1,0].text(160,-16,'October-March',fontsize=6)
    #ax[1,0].set_xlabel('Longitude',fontsize=6)
    #ax[1,0].set_ylabel('Latitude',fontsize=6)
    ax[1,0].tick_params(axis='both', which='major', labelsize=6)
    figGP = ax[1,0].figure; cb_ax = figGP.axes[7]; cb_ax.tick_params(labelsize=6,length = 0.05)
    ax[1,0].text(0.5, -0.5, 'Number of measurements', horizontalalignment='center', verticalalignment='center', transform=ax[1,0].transAxes,fontsize=7)#, weight='bold')
    
    plt.subplots_adjust(wspace=0,  
                    hspace=0.2) 
    
    plt.savefig(savepath + "/Results/Plots/Paper2021/final/GOSAT_Amount_viridis_log.png", dpi=400,bbox_inches='tight')
    plt.savefig(savepath + "/Results/Plots/Paper2021/final/GOSAT_Amount_viridis_log.pdf", dpi=400,bbox_inches='tight',format = 'pdf')
