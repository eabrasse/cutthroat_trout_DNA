#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot cutthroat trout DNA concentrations
"""

# importing modules, the beginning of all python code
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import pandas as pd
import string
import numpy as np

# some helpful plotting commands
plt.close('all')
tab10 = plt.get_cmap('tab10',10)

# load in data
home = '/Users/elizabethbrasseale/job materials/application to eDNA/'
data_fn = home+'cutthroat_trout_qPCR[UW2].csv'
df = pd.read_csv(data_fn,sep=',',engine='python')

# convert date from "321" to March 2021
df['date'] = pd.to_datetime(df['time'],format='%m%y')

# get lists of creeks, stations, and bioreps
creeks = df['creek'].unique()
creeks.sort()
stations = df['station'].unique()
stations = stations[::-1]
bioreps = df['biorep'].unique()
bioreps.sort()

# generate figure and axis handle
fig, axs = plt.subplots(len(bioreps),len(creeks),figsize=(12,6))

# loop through creeks
count=0
creek_count = 0
for creek in creeks:
    # create a dataframe of just data correpsonding to the current creek
    g = df[df['creek']==creek]
    
    # loop through bioreps
    biorep_count=0
    for biorep in bioreps:
        
        # create a dataframe of this biorep from dataframe for this creek
        gg = g[g['biorep']==biorep]
        
        # new axis for each biorep/creek combo
        ax = axs[biorep_count,creek_count]
        
        # keep time axis consistent for each creek
        if biorep_count>0:
            ax.sharex(axs[0,creek_count])
        
        # sort samples by date
        gg=gg.sort_values('date')
        
        # create two new dataframes representing upstream and downstream data
        gg_up = gg[gg['station']=='Up']
        gg_dn = gg[gg['station']=='Dn']
        
        # plot upstream and downstream data with different line/marker styles, color coded by biorep
        ax.plot(gg_up['date'],gg_up['DNA_concentation'],linestyle='solid',color=tab10(biorep_count),marker='o',mec=tab10(biorep_count),mfc=tab10(biorep_count))
        ax.plot(gg_dn['date'],gg_dn['DNA_concentation'],linestyle='dashed',color=tab10(biorep_count),marker='o',mec=tab10(biorep_count),mfc='none')
        
        # everything else is axis formatting & annotating
        # conditional statements are used to avoid redundant labeling & keep things pretty
        if creek_count==0:
            ax.set_ylabel('DNA concentration')
            ax.text(0.01,0.85,f'biorep {biorep}',fontweight='normal',fontsize=10,ha='left',va='top',transform=ax.transAxes,color=tab10(biorep_count))
        if biorep_count==0:
            ax.text(0.05,0.95,f'{creek}',fontweight='bold',fontsize=12,ha='left',va='top',transform=ax.transAxes)
        ax.grid()
            
        if biorep==bioreps[-1]:
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%-m/%y"))
            ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
            plt.setp( ax.xaxis.get_majorticklabels(), rotation=30, ha="right",rotation_mode='anchor')
            ax.set_xlabel('Date')
        else:
            plt.setp(ax.get_xticklabels(), visible=False)
            
        
        if count==0:
            ax.text(0.65,0.65,'Up',color=tab10(biorep_count),rotation=65,transform=ax.transAxes,fontweight='bold')
            ax.text(0.8,0.18,'Dn',color=tab10(biorep_count),rotation=0,transform=ax.transAxes,fontweight='bold')
        
        ax.set_ylim([-2,100])
        biorep_count+=1
        count+=1
    
    creek_count+=1

# show plot
plt.tight_layout()
plt.show(block=False)
plt.pause(0.1)
