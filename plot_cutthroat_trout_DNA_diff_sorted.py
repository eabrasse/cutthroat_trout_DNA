#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot cutthroat trout DNA concentrations change from upstream to downstream
as a function of upstream concentration
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
props = dict(boxstyle='round', fc='white',ec='None',alpha=1)
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
fig = plt.figure(figsize=(4,5.5))
ax = fig.gca()

# initialize lists for data
DNA_up_list = []
DNA_diff_list = []

# loop through creeks
creek_count = 0
for creek in creeks:
    # create a dataframe of just data correpsonding to the current creek
    g = df[df['creek']==creek]
    
    # skip 2Brn, it doesn't have two stations to compare
    if creek=='2Brn':
        continue
    
    # loop through bioreps
    biorep_count=0
    for biorep in bioreps:
        
        # create a dataframe of this biorep from dataframe for this creek
        gg = g[g['biorep']==biorep]
        
        # sort dates, since we'll use them to align upstream/downstream data
        gg=gg.sort_values('date')
        
        # create two new dataframes representing upstream and downstream data
        gg_up = gg[gg['station']=='Up']
        gg_dn = gg[gg['station']=='Dn']
        
        # get lists of dates - only want to include comparisons with data coincident at both stations
        date_list_up = [date for date in gg_up['date']]
        date_list_dn = [date for date in gg_dn['date']]
        
        # find samples in upstream data that don't have a downstream component, and vice versa
        mask_up = [date in date_list_dn for date in date_list_up]
        mask_dn = [date in date_list_up for date in date_list_dn]
        
        # convert dataframes to numpy array
        # dataframes are easier to plot, numpy arrays are easier for quantitative analysis
        DNA_up = np.array([dna for dna in gg_up['DNA_concentation']])
        DNA_dn = np.array([dna for dna in gg_dn['DNA_concentation']])
        
        # mask data
        DNA_up = DNA_up[mask_up]
        DNA_dn = DNA_dn[mask_dn]
        
        # store upstream data in a master array
        DNA_up_list= np.append(DNA_up_list,DNA_up,axis=0)
        
        # calculate normalized change in concentration from upstream to downstream
        # store difference in master array
        DNA_diff = (DNA_dn-DNA_up)/DNA_up
        DNA_diff_list = np.append(DNA_diff_list,DNA_diff,axis=0)
        
        # plot normalized change by initial concentration
        # do this here so that we can maintain color key for bio replicates
        ax.scatter(DNA_up,DNA_diff,color=tab10(biorep_count),facecolors='None',s=15,zorder=100)
        
        # generate a legend
        # I like doing this manually, I think it looks nicer
        if creek_count==0: # only do it once
            ax.text(0.9,0.55-0.05*biorep,f'biorep = {biorep}',fontweight='normal',fontsize=10,ha='right',va='top',transform=ax.transAxes,color=tab10(biorep_count),bbox=props)

        biorep_count+=1
    
    creek_count+=1

# BINNED MEANS
# calculate binned means
bin_edges = np.logspace(np.log(DNA_up_list.min()), np.log(DNA_up_list.max()), 10,base=np.e)
bins = np.exp(0.5*(bin_edges[:-1]+bin_edges[1:]))
digitized = np.digitize(DNA_up_list, bin_edges)
bin_means = [DNA_diff_list[digitized == i].mean() for i in range(1, len(bins))]

# plot binned means
ax.scatter(bin_edges[1:-1],bin_means,color='k',s=50,zorder=500,marker='s')

# add binned means to legend
ax.text(0.9,0.55-0.05*(biorep+1),'Binned means',fontweight='bold',fontsize=10,ha='right',va='top',transform=ax.transAxes,color='k',bbox=props)

# LINEAR FIT
# grab y-axis limits before plotting the line
ylim = ax.get_ylim()

# sort data in order by x-axis before fitting
# in this case, x-axis is DNA_up_list
DNA_up_inds = DNA_up_list.argsort()
DNA_up_list = DNA_up_list[DNA_up_inds]
DNA_diff_list = DNA_diff_list[DNA_up_inds]

# calculate linear fit, maintaining log scale of x-axis
DNA_up_log = np.log(DNA_up_list)
p = np.polyfit(DNA_up_log,DNA_diff_list,1)
y = p[1] + p[0]*DNA_up_log

# plot line
ax.semilogx(DNA_up_list,y,color='k',linestyle='dashed',zorder=150)

# maintain y-axis limits from before
# this is to keep the plot frame around the raw data, not the line
# which I consider an annotation, i.e. secondary info
ax.set_ylim(ylim)

# add linear fit to legend
ax.text(0.9,0.55-0.05*(biorep+2),'linear fit\n'+f'y = {p[1]:.2f} - {np.abs(p[0]):.2f}x',fontweight='normal',fontsize=10,ha='right',va='top',transform=ax.transAxes,color='k',bbox=props)

# AXIS FORMATTING
ax.set_ylabel(r'$(\mathrm{DNA}_{\mathrm{Up}}-\mathrm{DNA}_{\mathrm{Dn}})/\mathrm{DNA}_{\mathrm{Up}}$')
ax.grid()
ax.set_xlabel(r'$\mathrm{DNA}_{\mathrm{Up}}$')
ax.set_xscale('log')

# show plot
plt.tight_layout()
plt.show(block=False)
plt.pause(0.1)
