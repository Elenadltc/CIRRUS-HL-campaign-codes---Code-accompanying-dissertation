# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 12:37:37 2021

@ author: kirs_so
"""
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
import numpy as np
import os
import pandas as pd
import tkinter as tk

from pylab import cm
from matplotlib.colors import LogNorm
from Utility import (get_dNdD,
                     get_dNdD_mean,
                     get_dNdlogD_mean,
                     get_dNdlogD
                     )
from pandas.plotting import register_matplotlib_converters

# Register converters for using pandas in combination with matplotlib
register_matplotlib_converters()

current_working_directory = os.getcwd()

def check_start_end(timestart, timeend, data_df):
    
    if timestart == None:
        t_start = data_df.index[0]
    else:
        if data_df.index[0] > timestart:
            t_start = data_df.index[0]
        else:
            t_start = timestart
            
    if timeend == None:
        t_end = data_df.index[-1]
    else:
        if data_df.index[-1] < timeend:
            t_end = data_df.index[-1]
        else:
            t_end = timeend

    return t_start, t_end


def plot_Quicklook(
        data_df, bin_df, title, typ, log='off', campaign=None, members=None,
        timestart = None, timeend = None):
    """
    This function plots a timeseries of the particle size distribution, the
    total number concentration, and the LWC, MVD and EC derived from the PSD.
    
    Parameters
    ----------
    data_df: dataframe of floats
        containing the following columns:    
        dNdD: dataframe of floats
            number concentration per bin in m^-4
        N: iterable of floats
            total number concentration in m^-3
        MVD: iterable of floats
            Median volume diameter in m
        ED: iterable of floats
            Effective diameter in m
        LWC: iterable of floats
            Liquid Water Content in kg/cm³
    optionally also 
        IWC: iterable of floats
            Ice Water Content in kg/cm³
    bins: iterable of floats
        bin midpoints in mu
    timestart: datetime.datetime
        starting time of plot 
    timeend: datetime.datetime
        starting time of plot
    title: str
        title of plot
    category: str
        column name of total number concentration column
    
    Returns
    -------
    fig: figure
    """
    if campaign is None:
        campaign = ''

    edges = bin_df['bin_min'].to_numpy()
    mids = bin_df['bin_mid'].to_numpy()

    t_start, t_end = check_start_end(timestart, timeend, data_df) 

    data_df = data_df.loc[(data_df.index >= t_start) & (data_df.index <= t_end)]
            
    dNdlogD_m = get_dNdlogD_mean(data_df, bin_df *1e-6, t_start, t_end) * 1e-6
    dNdlogD = get_dNdlogD(data_df, bin_df*1e-6) * 1e-6
    dNdlogD[dNdlogD < 1e-15] = 1e-15
    dNdlogD = dNdlogD.loc[(data_df.index >= t_start) & (
            data_df.index <= t_end)].transpose()

    dNdD = get_dNdD(data_df) * 1e-12
    dNdD[dNdD < 1e-15] = 1e-15
    dNdD = dNdD[(data_df.index >= t_start) & (
            data_df.index <= t_end)]
    ###
    #The following creates a dataframe with bin edge to bin edge dNdlogD values
    bin_min = list(map(float,bin_df['bin_min']))
    bin_max = list(map(float,bin_df['bin_max']))
    bin_edges = [bin_min[0]]+bin_max
    double_edge=[bin_min[0]]
    double_edge.append(bin_max[-1])
    for w in range(1,len(bin_edges)-1):
        double_edge.insert(w,bin_edges[w])
        double_edge.insert(w+1,bin_edges[w])
    double_edge.sort()
    
    dNdlogD_m_float = list(map(float, dNdlogD_m))
    
    double_dlog = dNdlogD_m_float.copy()
    dlog_subset = dNdlogD_m_float[0:len(dNdlogD_m_float)].copy()
    for j in range(0,len(dlog_subset)):
        double_dlog.insert(j+j+1, dlog_subset[j])

    double_log_df  = pd.DataFrame({'Bin_edges': double_edge,
                                   'dNdlogD_m': double_dlog})
    
    dNdD = dNdD.transpose()
    dNdD_m = get_dNdD_mean(data_df, t_start, t_end) * 1e-12
    MVD = data_df['MVD'][(data_df.index >= t_start) & (
        data_df.index <= t_end)] * 1e6
    N = data_df['N'][(data_df.index >= t_start) & (
        data_df.index <= t_end)] * 1e-6
    ED = data_df['ED'][(data_df.index >= t_start) & (
        data_df.index <= t_end)] * 1e6
    LWC = data_df['LWC'][(data_df.index >= t_start) & (
        data_df.index <= t_end)] * 1e3
    if '(S-P)/(S+P)' in data_df.columns: 
        # This criterion only exists for the BCPD, for which the ratio
        # of the polarization directions shall be plotted instead of the
        # LWC
        pol_ratio = data_df['(S-P)/(S+P)'][(data_df.index >= t_start) & (
            data_df.index <= t_end)]

    if 'IWC' in data_df.columns:
        IWC = data_df['IWC'][(data_df.index >= t_start) & (
        data_df.index <= t_end)] *1e3 
        
    if typ == 'OAP':
        if log == 'off':
            t = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2]
        else:
            t = [1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3]
    elif typ == 'FS':
        if log =='off':
            t = [1e-3, 1e-2, 1e-1, 1, 1e1, 1e2]
        else:
            t = [1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3]
    else:
        print("""
              Please specify which typ your probe is. 'OAP' and 'FS' are 
              available.
              """)
        return
    
    fig = plt.figure(figsize=(28.5, 15.5))
    gs = gridspec.GridSpec(6, 4, width_ratios=[16, 1.1, 0.4, 3], 
                           height_ratios=[1, 1, 1, 1, 2, 0.3], wspace= 0.)
    axt = fig.add_subplot(gs[0, 0])
    axq = fig.add_subplot(gs[0, 1:])
    ax1 = fig.add_subplot(gs[1, 0])
    axa = fig.add_subplot(gs[1, 1:4], sharey=ax1)
    axa.yaxis.tick_right()
    ax2 = fig.add_subplot(gs[2, 0], sharex=ax1)
    axb = fig.add_subplot(gs[2, 1:4], sharey=ax2)
    axb.yaxis.tick_right()
    ax3 = fig.add_subplot(gs[3, 0], sharex=ax2)
    axc = fig.add_subplot(gs[3, 1:4], sharey=ax3)
    axc.yaxis.tick_right()
    ax0 = fig.add_subplot(gs[4, 0], sharex=ax3)
    axx = fig.add_subplot(gs[4, 1], sharey=ax0)
    axx.yaxis.tick_right()
    axy = fig.add_subplot(gs[4, 3])
    axy.yaxis.tick_right()
    axy.yaxis.set_label_position('right')
    axe = fig.add_subplot(gs[5, 0:])
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    
    author = ''
    if members is not None:
            author = ', '.join(name for name in members.keys())
    
    axt.text(0, 0.5, 'Quicklook ' + campaign + '      ' + title, fontsize = 40)
    axt.text(0, 0.2, 'preliminary data, only for quicklook use', fontsize = 10)
    axt.text(0, 0.0, author, fontsize = 12)
    axt.axis('off')
    
    img = mpimg.imread(os.path.dirname(os.path.abspath(__file__))\
                       + '/Icons/DLR.png')
    axq.imshow(img)
    axq.axis('off')

    ax1.plot(N, c='k')
    ax1.set_yscale('log')
    ax1.set_ylabel(r'N [$cm^{-3}$]', fontsize=12)
    ax1.set_title(title + ' ' + t_start.strftime('%d/%m/%Y') + ' ' +
    t_start.strftime('%H:%M:%S') + '-' + t_end.strftime('%H:%M:%S'))
    
    axa = plot_BOX(N, ['N'], fig=fig, ax=axa, sharey=True)
    
    mvd_line, = ax2.plot(MVD, c='k', label='MVD')
    ed_line, = ax2.plot(ED, c='b', label='ED')
    ax2.legend(handles=[mvd_line, ed_line], loc='upper right')
    ax2.set_ylabel(r'MVD, ED [$\mu$m]', fontsize=12)
    
    axb = plot_BOX([MVD[MVD!=0.], ED[ED!=0.]], ['MVD', 'ED'], fig=fig,
                    ax=axb, sharey=True)
    
    if '(S-P)/(S+P)' in data_df.columns: 
        # for the BCPD only
        ax3.plot(pol_ratio, c = 'k', label = '(S-P)/(S+P)')
        axc_name = ['(S-P)/(S+P)']
        axc_list = [pol_ratio]
        ax3.set_ylim(bottom = -1, top = 1)
        ax3.set_ylabel(r'(S-P)/(S+P)', fontsize = 12)
    else:     
        # Plot LWC and IWC
        LWC_line, = ax3.plot(LWC, c='k', label='LWC')
        axc_list = [LWC[LWC!=0.]]
        axc_name = ['LWC']
        if 'IWC' in data_df.columns:
            axc_list.append(IWC[IWC!=0.])
            axc_name.append('PIWC')
            IWC_line, = ax3.plot(IWC, '--', c='gray', label='PIWC')
            ax3.legend(handles=[LWC_line, IWC_line], loc='upper right')
            ax3.set_ylabel(r'LWC, PIWC [g/m$^3$]', fontsize=12)
        else:
            ax3.set_ylabel(r'LWC [g/m$^3$]', fontsize=12)
        ax3.set_ylim(bottom=0, top=5)
        ax3.set_yscale('log')
    
    axc = plot_BOX(axc_list, axc_name, fig=fig, ax=axc, sharey=True)
    

    for tick in ax0.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax0.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax2.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax3.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)

    if log == 'off':
        im = ax0.pcolormesh(data_df.index, edges, dNdD, norm=LogNorm(vmin=t[0], vmax=t[-1]), cmap=cm.jet,
                            rasterized=True, shading='flat')

        ax0.set_yscale('log')
        
        axx.semilogy(dNdD_m, edges, drawstyle='steps-pre', color='r',
                           label = 'PSD')
        axx.semilogy(dNdD_m, mids, 's', markersize=3, color='r')
        axx.set_xlabel(r'dN/dD [$\mathrm{cm}^{-3} \mathrm{\mu m}^{-1}$]')
        
        axy.semilogy(edges, dNdD_m, drawstyle='steps-post', color='r',
                           label = 'PSD')
        axy.semilogy(mids, dNdD_m, 's', markersize=3, color='r')
        axy.set_ylabel(r'dN/dD [$\mathrm{cm}^{-3} \mathrm{\mu m}^{-1}$]')
    elif log == 'on':

        im = ax0.pcolormesh(data_df.index, bin_edges, dNdlogD, norm=LogNorm(vmin=t[0], vmax=t[-1]), cmap=cm.jet)
        ax0.set_yscale('log')
        
        axx.semilogy(double_log_df['dNdlogD_m'], double_log_df['Bin_edges'], drawstyle='steps-mid', color='r',
                           label = 'PSD')
        axx.semilogy(dNdlogD_m_float, mids, 's', markersize=3, color='r')
        axx.set_xlabel(r'dN/dlogD [$\mathrm{cm}^{-3}$]')
        axy.semilogy(double_log_df['Bin_edges'], double_log_df['dNdlogD_m'], drawstyle='steps-mid', color='r',
                           label = 'PSD')
        axy.semilogy(mids, dNdlogD_m_float, 's', markersize=3, color='r')
        axy.set_ylabel(r'dN/dlogD [$\mathrm{cm}^{-3}$]')
    else:
        print("""
              Only use 'on' and 'off' for log option!
              """)
    
    ax0.set_ylabel(r'diameter [$\mu$m]', fontsize=14)
    ax0.xaxis_date()
    ax0.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    ax0.set_xlim([timestart, timeend])
    
    axx.set_xscale('log')
    axx.set_ylim((edges[0], edges[-1]))
    axx.legend(loc='lower left')
    
    axy.set_xlabel(r'diameter [$\mu$m]', fontsize=14)
    axy.set_xscale('log')
    axy.legend(loc='upper right')
    
    pos0 = ax0.get_position()
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.08, pos0.y0, 0.005, pos0.height])
    fig.colorbar(im, cax=cbar_ax, ticks=t)
    cbar_ax.yaxis.set_ticks_position('left')
    if log == 'off':
        cbar_ax.set_xlabel(r'[$\mathrm{cm}^{-3} \mathrm{\mu m}^{-1}$]')
    else:
        cbar_ax.set_xlabel(r'[$\mathrm{cm}^{-3}$]')
    
    email = ''
    if members is not None:
            email = ', '.join(name for name in members.values())
    
    axe.text(0, 0.0, email, fontsize = 8)
    axe.axis('off')

    return fig, t_start, t_end, mids, double_log_df['dNdlogD_m'], double_log_df['Bin_edges'] , dNdlogD_m_float

def plot_BOX(Series, xlabels, ylabels=None, fig=None, ax=None, sharey=None):
    """
    Input: 
    Series: pandas.Series or List of pandas.Series (obligatory)
    xlabels: list of strings (Name of each Graph) (obligatory)
    ylables: yaxes if not shared
    fig: the figure to which the plot will be connected
    ax: the subplot the plot will be connected
    sharey: True if plot will share y axes with other Subplot
    
    Funktion: 
    """
    if fig is None:
        fig = plt.figure()
    if ax is None:
        ax = fig.add_subplot()
        
    colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
              '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
              '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
              '#000075', '#808080', '#ffffff', '#000000']   
              
    bp = ax.boxplot(Series, showfliers=False, labels=xlabels, whis=[10,90],
               patch_artist=True, showmeans=True,
               meanprops=dict(marker='+', markeredgecolor='k',
                              markerfacecolor='tomato', markersize=13,
                              linewidth=20))
    
    for median in bp['medians']:
        median.set(color='snow', linewidth=2)
    for cap in bp['caps']:
        cap.set(color='black', linewidth=2)
    for whisker in bp['whiskers']:
        whisker.set(color='black', linewidth=2)
    for box,color in zip(bp['boxes'], colors):
        box.set(color='black', linewidth=2)
        box.set(facecolor=color)
        
    if sharey is None:
        if ylabels is not None:
            bp.set_ylabel(ylabels)
            bp.set_yscale('log')
            bp.set_ylim(min(Series)+0.1, max(Series))
    elif sharey is True:
        pass
    else:
        print("""
              Please only use sharey=True only if you want to share the y axes.
              Otherwise don't specify sharey!
              """)
    return bp
        
def plot_Aircraft(
        data_df, bin_df, title, typ, log='off', LWC_IWC='IWC', campaign=None, 
        members=None, timestart = None, timeend = None, mean=None):

    edges = bin_df['bin_min'].to_numpy()
    
    t_start, t_end = check_start_end(timestart, timeend, data_df)
    
    data_df = data_df.loc[(data_df.index >= t_start) & (data_df.index <= t_end)]

    dNdlogD = get_dNdlogD(data_df, bin_df*1e-6) * 1e-6
    dNdlogD[dNdlogD < 1e-15] = 1e-15
    dNdlogD = dNdlogD[(data_df.index >= t_start) & (
            data_df.index <= t_end)].transpose()

    dNdD = get_dNdD(data_df) * 1e-12
    dNdD[dNdD < 1e-15] = 1e-15
    dNdD = dNdD[(data_df.index >= t_start) & (
            data_df.index <= t_end)].transpose()

    N = data_df['N'][(data_df.index >= t_start) & (
        data_df.index <= t_end)] * 1e-6
    ED = data_df['ED'][(data_df.index >= t_start) & (
        data_df.index <= t_end)] * 1e6
    LWC = data_df['LWC'][(data_df.index >= t_start) & (
        data_df.index <= t_end)] * 1000
    if 'IWC' in data_df.columns:
        IWC = data_df['IWC'][(data_df.index >= t_start) & (
        data_df.index <= t_end)] *1e3 
    try:
        Alt = data_df['altitude'][(data_df.index >= t_start) & (
            data_df.index <= t_end)]
    except KeyError:
        Alt = None
    try:
        temp = data_df['static_temp'][(data_df.index >= t_start) & (
            data_df.index <= t_end)]
    except KeyError:
        temp = None

    if mean is not None:
        scale = str(mean) + 'S'
        N = N.resample(scale).mean()
        ED = ED.resample(scale).mean()
        LWC = LWC.resample(scale).mean()
        IWC = IWC.resample(scale).mean()
        
        
    if typ == 'OAP':
        if log == 'off':
            t = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2]
        else:
            t = [1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3]
    elif typ == 'FS':
        if log =='off':
            t = [1e-3, 1e-2, 1e-1, 1, 1e1, 1e2]
        else:
            t = [1e-1, 1, 1e1, 1e2, 1e3, 1e4]
    else:
        print("""
              Please specify which typ your probe is. 'OAP' and 'FS' are 
              available.
              """)
        return
    
    fig = plt.figure(figsize=(18.5, 15.5))
    gs = gridspec.GridSpec(7, 2, width_ratios=[3, 1], hspace=0.05, 
                           height_ratios=[1, 1, 1, 1, 1, 2, 0.3], wspace= 0.)
    axt = fig.add_subplot(gs[0, 0])
    axq = fig.add_subplot(gs[0, 1])
    ax6 = fig.add_subplot(gs[1, :])
    ax1 = fig.add_subplot(gs[2, :], sharex=ax6)
    ax2 = fig.add_subplot(gs[3, :], sharex=ax1)
    ax3 = fig.add_subplot(gs[4, :], sharex=ax2)
    ax0 = fig.add_subplot(gs[5, :], sharex=ax3)
    axe = fig.add_subplot(gs[6, :])
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax6.get_xticklabels(), visible=False)
    
    author = ''
    if members is not None:
            author = ', '.join(name for name in members.keys())
    
    axt.text(0, 0.6, 'Quicklook ' + campaign + '      ' + title, fontsize = 24)
    axt.text(0, 0.4, 'preliminary data, only for quicklook use', fontsize = 10)
    axt.text(0, 0.2, author, fontsize = 12)
    axt.axis('off')
    
    img = mpimg.imread(os.path.dirname(os.path.abspath(__file__))\
                       + '/Icons/DLR.png')
    axq.imshow(img)
    axq.axis('off')

    if Alt is not None:
        ax6.plot(Alt, c='k')
        ax6.set_ylabel(r'Altitude [m]', fontsize=14)
        ax6.set_title(title + ' ' + t_start.strftime('%d/%m/%Y') + ' ' +
        t_start.strftime('%H:%M:%S') + '-' + t_end.strftime('%H:%M:%S'))
        ax6.set_xlim([t_start, t_end])
    
    if temp is not None:
        ax7 = ax6.twinx()
        ax7.set_ylabel(r'Static Temperature [°C]', fontsize=14, color='orange')
        ax7.plot(temp, c='orange')
        ax7.tick_params(axis='y', labelcolor='orange')


    ax1.plot(N, c='r')
    ax1.set_ylabel(r'N [$cm^{-3}$]', fontsize=14)
    ax1.set_yscale('log')
    
    ed_line, = ax2.plot(ED, c='g', label='ED')
    ax2.set_yscale('log')
    ax2.set_ylim(bottom=3)
    ax2.set_ylabel(r'ED [$\mu$m]', fontsize=14)
    
    if LWC_IWC == 'LWC' or 'IWC' not in data_df.columns:
        LWC = ax3.plot(LWC, c='b', label='LWC')
        ax3.set_ylabel(r'LWC [g/m$^3$]', fontsize=14)
    elif LWC_IWC == 'IWC' and 'IWC' in data_df.columns:
        IWC = ax3.plot(IWC, c='b', label='IWC')
        ax3.set_ylabel(r'IWC [g/m$^3$]', fontsize=14)
    else:
        print('Please, specify if IWC or LWC should be plotted!')
        
    ax3.set_ylim(top=5)
    ax3.set_yscale('log')

    for tick in ax0.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax0.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax1.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax2.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax3.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)

    if log == 'off':
        im = ax0.pcolormesh(data_df.index, edges, dNdD, norm=LogNorm(vmin=t[0], vmax=t[-1]), cmap=cm.jet,
                            rasterized=True, shading='flat')
        ax0.set_yscale('log')
    elif log == 'on':
        im = ax0.pcolormesh(data_df.index, edges, dNdlogD, norm=LogNorm(vmin=t[0], vmax=t[-1]), cmap=cm.jet,
                            rasterized=True, shading='flat')
        ax0.set_yscale('log')
    else:
        print("""
              Only use 'on' and 'off' for log option!
              """)
    
    ax0.set_ylabel(r'diameter [$\mu$m]', fontsize=14)
    ax0.xaxis_date()
    ax0.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    ax0.set_xlim([t_start, t_end])
    
    pos0 = ax0.get_position()
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.07, pos0.y0, 0.002, pos0.height])
    fig.colorbar(im, cax=cbar_ax, ticks=t)
    cbar_ax.yaxis.set_ticks_position('left')
    if log == 'off':
        cbar_ax.set_xlabel(r'[$\mathrm{cm}^{-3} \mathrm{\mu m}^{-1}$]')
    else:
        cbar_ax.set_xlabel(r'[$\mathrm{cm}^{-3}$]')
    
    email = ''
    if members is not None:
            email = ', '.join(name for name in members.values())
    
    axe.text(0, 0.0, email, fontsize = 10)
    axe.axis('off')

    return fig, t_start, t_end

def plot_Housekeeping(hk_df,default, probe, campaign=None, members=None,
        timestart=None, timeend=None):
    
    if timestart is None:
        t_start = hk_df.index[0]
    else:
        if hk_df.index[0] > timestart:
            t_start = hk_df.index[0]
        else:
            t_start = timestart
            
    if timeend is None:
        t_end = hk_df.index[-1]
    else:
        if hk_df.index[-1] < timeend:
            t_end = hk_df.index[-1]
        else:
            t_end = timeend
    
    colors = ['b','g','r','c','m','y','k'];  
    
    # Used in the default case to bypass the Checkboxes
    class fake_obj:
        def get(self):
            return 1
    
    def var_plot():
        if default == False: 
            master.destroy()
        num_plots = 0; max_ind = len(var_obj)-1
        # Find out how many parameters shall be plotted
        for i in range(len(var_obj)): 
            num_plots = num_plots + var_obj[i].get(); 
            if num_plots > 3: 
                max_ind = i; 
                break
        if num_plots == 0: 
            return
        var_obj_short = var_obj[:max_ind+1]; 
        
        
        fig = plt.figure(figsize=(28.5, 15.5))
        gs = gridspec.GridSpec(num_plots, 1, wspace= 0.)
        
        # Create axis in loop
        ax_obj = []; 
        for i, var in enumerate(var_obj_short): 
            if var.get():
                if len(ax_obj) == 0: 
                    ax = fig.add_subplot(gs[num_plots-1,0])
                else: 
                    ax = fig.add_subplot(gs[num_plots-1-len(ax_obj),0],
                                            sharex=ax_obj[0])                
                ax.plot(hk_df.iloc[:,i], c=colors[len(ax_obj)])                
                ax.set_ylabel(hk_df.columns[i], fontsize=12)
                #ax.set_title(title + ' ' + t_start.strftime('%d/%m/%Y') + ' ' +
                #t_start.strftime('%H:%M:%S') + '-' + t_end.strftime('%H:%M:%S'))
                ax_obj.append(ax);
                
        ax_obj[0].xaxis_date()
        ax_obj[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
        ax_obj[0].set_xlim([timestart, timeend])
    
    
    # Create var names: 
    var_names = []
    var_obj = []
    
    
    if default == True:
        for i in range(0,len(hk_df.columns)):
            var_obj.append(fake_obj()) 
        var_plot()        
    
    else: 
        # Create checkbox object
        master = tk.Tk()
        master.geometry("500x450+300+300")
        
        # Create checkboxes with the names of the housekeeping parameters
        for i in range(0,len(hk_df.columns)): 
            var_names.append('var'+str(i)); 
            var = tk.IntVar()
            var_obj.append(var)
            tk.Checkbutton(master, text=hk_df.columns[i],\
                        variable=var).grid(row=int(i/2), column=i%2, sticky=tk.W)        
        
        tk.Button(master, text='OK',
                  command=var_plot).grid(row=int(len(hk_df.columns)/2)+2,
                                        sticky=tk.W,padx=50, pady=10)
        tk.mainloop()

def plot_PBP_Hist(df, column, bins, t_start, t_end, density, logx):
    values = df[(df.loc[:,'Datetime'] >= t_start) & (df.loc[:,'Datetime'] <= t_end)][column]
    plt.figure(figsize=(10,7))
    plt.hist(values, bins, density = density)
    plt.xlabel(column)
    if logx == True:
        plt.gca().set_xscale("log")

