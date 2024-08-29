# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 16:31:56 2021

Import Housekeeping data into Data Object

Input: .csv file
Output: hf_df 

@author: dela_el
"""

import numpy as np
import pandas as pd
from Utility import timestamp_dmt

def get_cdp_housekeeping(csv_file):
    
    # Line where the recorded data starts still has to be determined, until 
    # then it is set to nan
    data_start_line = np.nan
    
    with open(csv_file) as f: 
        for i, line in enumerate(f.readlines()): 
            if '****' in line: 
                data_start_line = i  # CSV file starts with index 1
                break
            else: 
                continue
    
    # Read all the data into data
    data = pd.read_csv(csv_file, sep=',', header=1, skiprows=data_start_line,
                       index_col='Time')
   
    # Get timestamped data from data 
    timestamp = timestamp_dmt(data)
    
    # Get Housekeeping data from the data frame
    
    hk = data.loc[:, 'Laser Current (mA)':'Control Board T (C)']
    
    hk_df = hk.set_index(timestamp)
    
    hk_df['Applied PAS (m/s)'] = data.loc[:, 'Applied PAS (m/s)']
    
    return hk_df

def get_cas_housekeeping(csv_file):
    
    # Line where the recorded data starts still has to be determined, until 
    # then it is set to nan
    data_start_line = np.nan
    
    with open(csv_file) as f: 
        for i, line in enumerate(f.readlines()): 
            if '****' in line: 
                data_start_line = i  # CSV file starts with index 1
                break
            else: 
                continue
    
    # Read all the data into data
    data = pd.read_csv(csv_file, sep=',', header=1, skiprows=data_start_line,
                       index_col='Time')
   
    # Get timestamped data from data 
    timestamp = timestamp_dmt(data)
    
    # Get Housekeeping data from the data frame
    
    hk = data.loc[:, 'Recovery Temp (C)':'Back Block Temp (C)']
    hk = hk.join(data.loc[:, 'Qual TEC Temp (C)':'Back High Gain Baseline (V)'])
    hk = hk.join(data.loc[:, 'Dpol Baseline (V)':'Laser Monitor'])
    
    # Convert to numpy array, then use the array to create a new DataFrame 
    hk_df = hk.set_index(timestamp)
    
    hk_df['Applied PAS (m/s)'] = data.loc[:, 'Applied PAS (m/s)']
    
    return hk_df

def get_pip_housekeeping(csv_file):
    
    # Line where the recorded data starts still has to be determined, until 
    # then it is set to nan
    data_start_line = np.nan
    
    with open(csv_file) as f: 
        for i, line in enumerate(f.readlines()): 
            if '****' in line: 
                data_start_line = i  # CSV file starts with index 1
                break
            else: 
                continue
    
    # Read all the data into data
    data = pd.read_csv(csv_file, sep=',', header=1, skiprows=data_start_line,
                       index_col='Time')
   
    # Get timestamped data from data 
    timestamp = timestamp_dmt(data)
    
    # Get Housekeeping data from the data frame
    
    hk = data.loc[:, 'Diode 1 Voltage':'Diode 32 Voltage']
    hk = hk.join(data.loc[:, 'Laser Temp (C)':'Laser Power'])
    hk = hk.join(data.loc[:, 'PIP T Ambient (C)':'PIP PAS (m/s)'])
    
    # Convert to numpy array, then use the array to create a new DataFrame 
    hk_df = hk.set_index(timestamp)
    
    hk_df['Applied PAS (m/s)'] = data.loc[:, 'Applied PAS (m/s)']
    hk_df['DSP Temp (C)'] = data.loc[:, 'DSP Temp (C)']
    
    return hk_df

def get_cip_housekeeping(csv_file):
    
    # Line where the recorded data starts still has to be determined, until 
    # then it is set to nan
    data_start_line = np.nan
    
    with open(csv_file) as f: 
        for i, line in enumerate(f.readlines()): 
            if '****' in line: 
                data_start_line = i  # CSV file starts with index 1
                break
            else: 
                continue
    
    # Read all the data into data
    data = pd.read_csv(csv_file, sep=',', header=1, skiprows=data_start_line,
                       index_col='Time')
   
    # Get timestamped data from data 
    timestamp = timestamp_dmt(data)
    
    # Get Housekeeping data from the data frame
    
    hk = data.loc[:, 'Diode 1 Voltage':'Diode 64 Voltage']
    hk = hk.join(data.loc[:, 'Board Temp (C)':'Laser Current (mA)'])
    hk = hk.join(data.loc[:, 'CIP Grayscale T Ambient (C)':'CIP Grayscale PAS (m/s)'])
    
    # Convert to numpy array, then use the array to create a new DataFrame 
    hk_df = hk.set_index(timestamp)
    
    hk_df['Applied PAS (m/s)'] = data.loc[:, 'Applied PAS (m/s)']
    hk_df['Laser Temp (C)'] = data.loc[:, 'Laser Temp (C)']
    hk_df['Recovery Temp (C)'] = data.loc[:, 'Recovery Temp (C)']
    
    return hk_df