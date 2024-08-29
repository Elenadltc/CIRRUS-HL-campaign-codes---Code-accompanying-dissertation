# -*- coding: utf-8 -*-
"""
Data Object

Simon Kirschler
"""
import os
import sys
import numpy as np
import pandas as pd

from Utility import (load_data, save_data, determine_type,
                     get_aircraft_data, df_resample_fromdf)
from Plots import (plot_Quicklook,
                   plot_Aircraft, 
                   plot_Housekeeping
                   )
#from Plots import plot_LWC, plot_Housekeeping
from housekeeping import (get_cdp_housekeeping, get_cas_housekeeping,
                          get_pip_housekeeping, get_cip_housekeeping)


class Data():
    
    # Initialize the data object
    def __init__(self): 
        # Directory (path) in which the object shall be saved
        self.path = ''
        # Name under which the data shall be saved
        self.name = ''
        # Start time 
        self.t_start = None
        # End time
        self.t_end = None
        # Bin number of lowest size bin that shall be used (starting from 1)
        self.min_bin = None
        # Bin number of highest size bin that shall be used (starting from 1)
        self.max_bin = None
        # Probe type (e.g CIP, 2DS, CDP ...)
        self.probe = ''
        # Type (e.g OAP, FS)
        self.type = ''
        # Directory where the object data come from
        self.source = ''
    ''' More attributes can be added. Make sure to also add them in the save 
    and load functions in Utility.py, so that the Data object can be properly
    stored '''    
    
    def get_mean_value(self,df_key,t_start=None,t_end=None): 
        ''' Computes the mean value of the data frame column 'df_key' between 
        t_start and t_end'''
        if t_start == None: 
            t_start = self.t_start
        if t_end == None: 
            t_end = self.t_end
        # Find indices of t_start and t_end in data_df
        start_index = pd.Index(self.data_df.index).get_loc(str(t_start),
                               method='nearest')
        end_index = pd.Index(self.data_df.index).get_loc(str(t_end),
                             method='nearest')
        # Handle the different data types that start_index and end_index may 
        # have. 
        if type(start_index) == slice: 
            start_index = start_index.stop
        elif type(start_index) == np.ndarray: 
            start_index = start_index[0]
        else: 
            start_index = start_index
        if type(end_index) == slice: 
            end_index = end_index.start
        elif type(end_index) == np.ndarray:
            end_index = end_index[0]
        else: 
            end_index = end_index
        val = np.mean(self.data_df[df_key][start_index:end_index])
        return val
    
    def import_aircraft_data(self, aircraft_files, freq=1):
        
        list_data = []
        for i, file in enumerate(aircraft_files):
            aircraft_part = get_aircraft_data(file)
            list_data.append(aircraft_part)
        aircraft_data = pd.concat(list_data)
        aircraft_res = df_resample_fromdf(self.data_df, aircraft_data, freq)
        if self.probe == 'CDP' or self.probe == 'CAS':
            self.data_df = df_resample_fromdf(aircraft_res, self.data_df, freq)
            self.data_df = self.data_df.join(aircraft_res)
            self.data_df = df_resample_fromdf(aircraft_res, self.data_df, freq)
            self.data_df = self.data_df.fillna(0)
        else:
            self.data_df = self.data_df.join(aircraft_res).fillna(0)
    
    def save(self):
        save_data(self)
        
    def plot_quicklook(self,t_start=None, t_end=None, campaign='', members=None, log='off'):

        if t_start is None:
            timestart = self.t_start
        else:
            timestart = t_start

        if t_end is None:
            timeend = self.t_end
        else:
            timeend = t_end
        
        self.fig, t_start, t_end, self.bin_mids, self.double_log_dNdlogD_m, self.double_bin_edges, self.dNdlogD_m_float  = plot_Quicklook(
                self.data_df, self.bin_df, self.probe, self.type, log=log,
                campaign=campaign, members=members,
                timestart=timestart, timeend=timeend)

    def plot_aircraft(self,t_start=None, t_end=None, campaign='', members=None, 
                      log='off', LWC_IWC='IWC', mean=None):

        if t_start is None:
            timestart = self.t_start
        else:
            timestart = t_start

        if t_end is None:
            timeend = self.t_end
        else:
            timeend = t_end
        
        self.fig, t_start, t_end = plot_Aircraft(
                self.data_df, self.bin_df, self.probe, self.type, log=log,
                LWC_IWC=LWC_IWC, campaign=campaign, members=members, mean=mean,
                timestart=timestart, timeend=timeend)

    def plot_housekeeping(self, hk_file=None, default=True):
        
        if self.probe == 'CDP':
            if hk_file is None:
                self.hk_df = get_cdp_housekeeping(self.source)
            else:
                self.hk_df = get_cdp_housekeeping(hk_file)
        
        elif self.probe == 'CAS':
            self.hk_df = get_cas_housekeeping(self.source)
        
        # Housekeeping in PIP are in csv file and not in the nc file saved in source
        elif self.probe == 'PIP':
            self.hk_df = get_pip_housekeeping(hk_file)
            
        elif self.probe == 'CIPG':
            self.hk_df = get_cip_housekeeping(hk_file)
            
        plot_Housekeeping(self.hk_df, default, self.probe,
                          timestart = self.t_start, timeend = self.t_end)
        
        
def load(path): 
    # Create a data object, then pass it to the function load data in 
    # Utility.py, where the data is added to the Data object
    d = Data()    
    d = load_data(d,path)
    return d 
    
        
        