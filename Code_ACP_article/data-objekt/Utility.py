# -*- coding: utf-8 -*-
"""
Utilities for Data Object

Simon Kirschler
"""
import os
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import netCDF4 as nc


'''NOTE: ALL INPUTS AND OUTPUTS SHOULD BE IN SI UNITS'''


def average_data_df(data_df,averaging_time=5):
    new_data_df = data_df.resample(str(averaging_time)+'S',
                                   convention='end').mean()
    return new_data_df
            
        
    

def bin_df_from_edges(bin_edges):
    ''' Input: Bin edges including lowermost and uppermost bin edge. That means
    for n bins n+1 bin edges must be provided. '''

    bin_df = pd.DataFrame({'bin_min': bin_edges[:-1],
                        'bin_max': bin_edges[1:]},
                        index=np.arange(1,len(bin_edges)))
    
    bin_df['bin_mid'] = (bin_df['bin_min'] + bin_df['bin_max'])/ 2
    bin_df['bin_width'] = bin_df['bin_max'] - bin_df['bin_min']

    return bin_df

def bin_to_1um_res(dn,dNdD,bin_df,logarithmic=True): 
    '''
    Inputs: 
    Pandas data frames bin_df, dn and dNdD, plus the option to 
    select between logarithmic and linear interpolation    
    
    
    Outputs: 
        
    Pandas data frame bin_df at 1 um resolution and Pandas data frames
    dn_1um (the counts at 1 um resolution) and dNdD_1um, number concentration
    per m^4 at 1 um resolution
    
    '''
    # Find upper and lower limit of the 1um bin df
    bottom_edge = np.min(bin_df['bin_min']); 
    top_edge = np.max(bin_df['bin_max']); 
    
    # Compute the bin edges
    edges_bin = np.arange(bottom_edge, top_edge+1)
    # From the bin edges find the new bin_df
    bin_df_1um = bin_df_from_edges(edges_bin)  
    
    # Create empty array for dndD, the counts divided by the bin width
    dndD = np.empty(np.shape(dn))
                                            
    # Create empty arrays for the new dNdD_1um and dn_1um arrays 
    dNdD_1um = np.zeros([len(dNdD.index),len(edges_bin)-1])
    dn_1um = np.zeros([len(dn.index),len(edges_bin)-1]) 

    # Loop along the time axis                                   
    for t in range(0,len(dNdD.index)): 
        # The counts have to be transformed to counts divided by bin width
        # before they can be interpolated to 1um resolution
        dndD[t,:] = dn.to_numpy()[t,:]/bin_df['bin_width'].to_numpy()
        # The data can be interpolated linearly or logarithmically
        if logarithmic == True: 
            dNdD_1um[t,:] = log_interp(bin_df_1um.loc[:,'bin_mid'],
                    bin_df.loc[:,'bin_mid'],dNdD.iloc[t,:])
            dn_1um[t,:] = log_interp(bin_df_1um['bin_mid'].to_numpy(),
                  bin_df['bin_mid'].to_numpy(),dndD[t,:])
        else:
            dNdD_1um[t,:] = np.interp(bin_df_1um.loc[:,'bin_mid'].to_numpy(),
                    bin_df.loc[:,'bin_mid'].to_numpy(),dNdD.iloc[t,:])
            dn_1um[t,:] = np.interp(bin_df_1um['bin_mid'].to_numpy(),
                  bin_df['bin_mid'].to_numpy(),dndD[t,:])
    
    # Remove nan values
    dNdD_1um[np.isnan(dNdD_1um)] = 0; 
    dn_1um[np.isnan(dn_1um)] = 0; 
    
    # Get the column names
    dNdD_1um_names,dn_1um_names = column_names(bin_df_1um)
    
    # Create dNdD data frame at 1 um resolution
    dNdD_1um_df = pd.DataFrame(data = dNdD_1um, index = dNdD.index, 
                               columns = dNdD_1um_names)
    # Create dN data frame at 1 um resolution
    dn_1um_df = pd.DataFrame(data = dn_1um, index = dn.index,
                             columns = dn_1um_names)
    
    return bin_df_1um, dNdD_1um_df, dn_1um_df

def calc_ED(data_df, bin_df):
    """ Calculate the effective diameter in mu for OAP data
    Parameters
    ----------
    data_df: Must at least contain dNdD
    
    bin_df: bin boundaries, widths and midpoint, all in [m]
    
    Returns
    -------
    ED : float or iterable of floats
         effective diameter in mu
    
    Notes
    -----
    
    References
    ----------
    http://www.dropletmeasurement.com/PADS_Help/ED_(Effective_Diameter)_in_um.htm
    
    """
    dN = get_dN(data_df, bin_df)
    bin_mid = bin_df['bin_mid']
    
    ED = np.empty(np.shape(dN)[0])

    for i in range(len(dN)):
        if np.nansum(dN[i]) != 0.:
            ED[i] = (2 * np.nansum(dN[i] * (bin_mid / 2)**3) / \
                      np.nansum(dN[i] * (bin_mid / 2)**2))
        else:
            ED[i] = 0.

    return ED # in [m]

def calc_LWC(dN, bin_mid):
    """ Calculate the liquid water content
    Parameters
    ----------
    dN : float or iterable of floats 
         Aerosol particle concentration in each bin in 1/m³
    bin_mid: bin centers in m !
        Can be taken from bin_df['bin_mid'].
        DO NOT FORGET TO CONVERT TO METERS!
          
    Returns
    -------     
    LWC : float or iterable of floats
          total LWC in kg/m³
          
    LWC_i : float or iterable of floats
          LWC in each bin in kg/m³
    Notes
    -----
    Density of water is 1g/cm3
    LWC = sum(dN_i * pi/6*m_i³ * 10^-12)
    
    References
    ----------
    http://www.dropletmeasurement.com/PADS_Help/LWC_(Liquid_Water_Content)_in_g_cm%5E3.htm'
    
    """
    rho_water = 1000 # Density of water [kg/m³]
    LWC_i = np.empty(np.shape(dN))
    LWC = np.empty(np.shape(dN)[0])
    for i in range(np.shape(dN)[0]):           # time axes
        LWC_i[i] = dN[i] * rho_water * np.pi / 6 * bin_mid**3
        LWC[i] = np.nansum(LWC_i[i])

    return LWC, LWC_i

def calc_MVD(data_df, bin_df):
    """ Calculate the median volume diameter (MVD) in mu. The MVD is the 
    diameter where 50% of the LWC are in particles smaller than this diameter.
    
    Parameters
    ----------
    data_df : Data Frame that contains counts per bin (dn) and counts per 
              cubic meter ber micrometer (dNdD)
    bins : float or iterable of floats 
           Size bin edges in mu
    
    Returns
    -------
    MVD : float or iterable of floats
          Median volume diameter in mu
    
    Notes
    -----
   
    
    """
    
    dN = get_dN(data_df, bin_df)
    bin_min = bin_df['bin_min'].to_numpy()
    bin_max = bin_df['bin_max'].to_numpy()
    bin_mid = bin_df['bin_mid'].to_numpy()
    [LWC, LWC_i] = calc_LWC(dN, bin_mid)
    MVD = np.empty(np.shape(dN)[0])
    
    # Loop over the time axes
    for i in range(np.shape(dN)[0]):
        # Make sure that LWC[i] is not zero, otherwise NANs will be present in 
        # the MVD           
        if LWC[i] != 0:
            # Compute the cumulative LWC
            cum_LWC_i = [np.nansum(LWC_i[i,:j]) / LWC[i]\
                         for j in range(1,len(LWC_i[i])+1)]
        else:
            cum_LWC_i = [0] * (len(LWC_i[i]))
        # Obtain MVD by interpolating the bin_max at the position cum_i = .5
        MVD[i] = np.interp(.5, [0] + cum_LWC_i, [bin_min[0]] + list(bin_max),
                           left=0, right=0)
            
    return MVD


def column_names(bin_df):

    length = len(bin_df)

    dn_names = []
    dNdD_names = []
    for i in range(1,length + 1):
        if i < 10:
            dn_names.append('dn_00' + str(i))
            dNdD_names.append('dNdD_00' + str(i))
        elif 10 <= i < 100:
            dn_names.append('dn_0' + str(i))
            dNdD_names.append('dNdD_0' + str(i))
        elif i >= 100:
            dn_names.append('dn_' + str(i))
            dNdD_names.append('dNdD_' + str(i))

    return dNdD_names, dn_names
    
def convert_time(year, days, seconds):
    'Returns the timestamp for the data frame'
    timestamp = pd.Timestamp(int(year),1,1) + pd.Timedelta(days=int(days-1), 
                            seconds=int(seconds),
                            microseconds=int(np.mod(seconds,1)*1e6))

    return timestamp

def pads_time_to_datetime(date,pads_seconds):
    pads_time = timedelta(seconds = pads_seconds)
    pads_time = (datetime.min + pads_time).time()
    return datetime.combine(date,pads_time)

def determine_type(probe):
        ''' Determines whether the probe is an OAP or FS probe'''
        if 'CIP' in probe: 
            typ = 'OAP'
        elif 'CDP' in probe: 
            typ = 'FS'        
        elif 'BCPD' in probe: 
            typ = 'BS'  # backscatter
        elif 'CAS' in probe:
            typ = 'FS'
        elif 'CAPS' in probe:
            typ = 'FS'
        elif 'PIP' in probe:
            typ = 'OAP'
        elif '2DS' in probe:
            typ = 'OAP'
        return typ
    
def df_resample_fromdf(df_source, df, freq):
    """
    Resample a dataframe (df) to be adapted to another dataframe (df_source)
    with a certain sampling frequency

    Parameters
    ----------
    df_source : pd.DataFrame
        Dataframe source
    df : pd.DataFrame
        Dataframe to be resampled
    freq : int
        sampling frequency

    Returns
    -------
    df_resampled : pd.DataFrame
        df after resampling process

    """
    
    df = df.resample(str(freq) + 'S').mean()
    start = df_source.index[0]
    end = df_source.index[-1]
    df_resampled = df[start:end]
    
    df_resampled = df.resample(str(freq) + 'S', 
                           origin = start).mean()
    
    return df_resampled
    
def get_aircraft_data(nc_file):
    """
    Get altitude, static temperature and static pressure pandas DataFrame
    from aircraft data in netCDF format

    Parameters
    ----------
    nc_file : string
        path to aircraft data nc file

    Returns
    -------
    aircraft_df : pandas DataFrame
        Dataframe with altitude, static temperature and static pressure

    """
    
    f = nc.Dataset(nc_file)
    time = f.variables['TIME'][:]
    
    aircraft_df = pd.DataFrame({'altitude':f.variables['H'][:],
                                'static_temp':(f.variables['TS'][:] -273.15),
                                'static_p': f.variables['PS'][:]}, index=time)
    aircraft_df.index = pd.to_datetime(aircraft_df.index, unit='s')
    
    return aircraft_df

def get_dN(data_df, bin_df):
    """
    Input:  OAP Daten Dictionary
    
    Funktion:   gibt die Konzentration per Bin multipliziert mit der jeweiligen
                Binbreite aus. Um von 1/m^4 auf 1/m^3 zu kommen
            
    Output: dN [1/m^3]
    Format: numpy.array
    """
    bin_width = bin_df['bin_width'].to_numpy()
    dNdD = get_dNdD(data_df).to_numpy()
    dN = np.empty(np.shape(dNdD))
    for i in range(len(dN[0])):
        for p in range(len(dN)):
            dN[p, i] = dNdD[p, i] * bin_width[i]
    return dN


def get_dNdD(df):
    """
    Input:  OAP dataframe
    Format  pandas.Dataframe
    
    Funktion:   gibt die Konzentration per Bin als Dataframe aus
    
    Output: dNdD dataframe
    Format: pandas.Dataframe
    """
    dNdD_names = []
    for i in df.columns:
        if i[:2] == 'dN':
            dNdD_names.append(i)
    dNdD = df.loc[:, dNdD_names]
    return dNdD

def get_dlogD(bin_df):
    """
    Input: bin_df
    
    Output: list of logarithmic bin widths
    """
    bin_min = bin_df['bin_min'].to_numpy()
    bin_max = bin_df['bin_max'].to_numpy()
    
    dlogD = []
    for x in range(len(bin_min)):
        dlogD.append(np.log10(bin_max[x] / bin_min[x]))
    
    return dlogD


def get_dlnD(bin_df):
    """
    Input: bin_df

    Output: list of logarithmic bin widths
    """
    bin_min = bin_df['bin_min'].to_numpy()
    bin_max = bin_df['bin_max'].to_numpy()

    dlogD = []
    for x in range(len(bin_min)):
        dlogD.append(np.log(bin_max[x] / bin_min[x]))

    return dlogD

def get_dNdlogD(data_df, bin_df):
    """
    
    """
    dlogD = get_dlogD(bin_df)
    dN = get_dN(data_df, bin_df)
    
    pdN = get_pandas_dN(data_df, dN)
    dNdlogD = pd.DataFrame().reindex_like(pdN)
    
    for i in range(len(pdN.columns)):
        dNdlogD[dNdlogD.columns[i]] = pdN[pdN.columns[i]] / dlogD[i]
        
    return dNdlogD

def get_dNdlogD_mean(data_df, bin_df, t_start, t_end):
    """
    Input:
    data_df: pandas.DataFrame 
    t_start: datetime.datetime
    t_end: datetime.datetime
    
    Funktion: Berechnet die gemittelte Anzahlkonzentration [#/m³]
              
    Output:
    dN_mean: numpy.array mit Anzahl Bin Einträge
    """
    dNdlogD = get_dNdlogD(data_df, bin_df)
    
    dNdlogD_mean = np.zeros(np.shape(dNdlogD)[1])
                
    for x in range(len(dNdlogD.columns)):
        dNdlogD_mean[x] = np.mean(dNdlogD[dNdlogD.columns[x]][t_start:t_end])
                        
    return dNdlogD_mean

def get_dNdD_mean(data_df, t_start, t_end):
    """
    Input:
    data_df: pandas.DataFrame 
    t_start: datetime.datetime
    t_end: datetime.datetime
    
    Funktion: Berechnet die gemittelte Anzahlkonzentration [#/m³] pro m über 
              Zeitraum.
              
    Output:
    dNdD_mean: numpy.array mit Anzahl Bin Einträge
    """
    dNdD = get_dNdD(data_df)
    
    dNdD_mean = np.zeros(np.shape(dNdD)[1])
    
    for x in range(len(dNdD.columns)):
        dNdD_mean[x] = np.mean(dNdD[dNdD.columns[x]][t_start:t_end])
        
    return dNdD_mean

def get_pandas_dN(data_df, dN):
    """
    
    """
    dNdD = get_dNdD(data_df)

    pdN = pd.DataFrame().reindex_like(dNdD)
    for i in range(len(pdN.columns)):
        pdN[pdN.columns[i]] = dN[:, i]  
    
    return pdN

def get_variable(var_name, csv_file):
    """
    Search for a certain variable in a given csv file

    Parameters
    ----------
    var_name : string
        name of the desired variable
    csv_file : string
        path to csv file

    Returns
    -------
    variable : pandas DataFrame
        data from csv file corresponding to the searched variable
        End Seconds from csv file

    """
    
    with open(csv_file) as f:
        for i, line in enumerate(f.readlines()):
            if '****' in line:
                data_start_line = i
    data = pd.read_csv(csv_file, sep=',',header=1,skiprows=data_start_line,
                       index_col='Time')
    variable = data[var_name]
    variable = pd.DataFrame({var_name:variable, 'End Seconds': data['End Seconds']}, index=data.index)
    
    return variable

def timestamp_dmt(df_raw, typ=None, seconds=None):
    '''
    Input: Data frame with timestamp in DMT format (Year, Day of Year, End 
    Seconds)
    
    Output: Timestamp in datetime format
    '''
    if seconds == None:
        seconds = 'End Seconds'

    if typ=='CAPS':
        return df_raw.apply(lambda x: convert_time(x['Year'], x['Day Of Year'], 
                                                   x[seconds]), axis=1,
                                                    result_type='expand')    
    else:
        return df_raw.apply(lambda x: convert_time(x['Year'], x['Day of Year'], 
                                                   x[seconds]), axis=1,
                                                    result_type='expand')

def load_data(data_object, path): 
    ''' 
    Loads data from the data files contained in 'path' into the Data object
    data_object.
    
    Output: Data object with attributes, bin_df and data_df
    '''
    
    # Get the names of all the files in the directory (third argument 
    # of os.walk)
    filenames = [f for f in os.listdir(path) 
                if os.path.isfile(os.path.join(path, f))]
    # Get the names of the info, bin and data files
    for filename in filenames: 
        if 'INFO' in filename:
            infofile = filename
        elif 'BIN_DF' in filename:
            binfile = filename
        elif 'DATA_DF' in filename:
            datafile = filename
    
    # Get the data from the info file
    f = open(path + '/' + infofile,'r')
    for line in f.readlines(): 
        if 'PATH' in line: 
            data_object.path = line.split('=')[1].strip()
        elif 'NAME' in line: 
            data_object.name = line.split('=')[1].strip()
        elif 'PROBE' in line: 
            data_object.probe = line.split('=')[1].strip()
        elif 'SOURCE' in line: 
            data_object.source = line.split('=')[1].strip()
        elif 'TYPE' in line: 
            data_object.type = line.split('=')[1].strip()
        elif 'STARTTIME' in line: 
            try: 
                data_object.t_start = datetime.strptime(line.split('=')[1].strip(),
                                                        '%Y-%m-%d %H:%M:%S.%f')
            except ValueError: 
                data_object.t_start = datetime.strptime(line.split('=')[1].strip(),
                                                        '%Y-%m-%d %H:%M:%S')
        elif 'ENDTIME' in line: 
            try: 
                data_object.t_end = datetime.strptime(line.split('=')[1].strip(),
                                                      '%Y-%m-%d %H:%M:%S.%f')
            except ValueError: 
                data_object.t_end = datetime.strptime(line.split('=')[1].strip(),
                                                      '%Y-%m-%d %H:%M:%S')
        elif 'MINBIN' in line: 
            data_object.min_bin = float(line.split('=')[1].strip())
        elif 'MAXBIN' in line: 
            data_object.max_bin = float(line.split('=')[1].strip())
        else: 
            continue
    
    f.close()
    
    # Read in the bin_df
    data_object.bin_df = pd.read_csv(path + '/' + binfile, index_col=0)
    
    # Read in the data_df
    data_object.data_df = pd.read_csv(path + '/' + datafile, index_col=0, 
                                      parse_dates=True)
    
    return data_object


def log_interp(zz, xx, yy):
    ''' Logaritmic interpolation at the locations zz, of the data (xx,yy)''' 
    logz = np.log10(zz)
    logx = np.log10(xx)
    logy = np.log10(yy)
    return np.power(10.0, np.interp(logz, logx, logy))


def save_data(data):
    '''
    The user has the possibility to set a path and a name for the data 
    object. If no path is set, a new folder containing files related to the 
    data object will be created in the current working directory. If no name is
    set the timestamp of the beginning of the measurement will be used as 
    the name. A probe identifier (e.g. 'CDP') will always be added to the 
    beginning of the directory and filename. 
    '''
    if data.path == '':
        path = ''
    else: 
        path = data.path + '/'
        
    if data.name == '':
        name = data.probe + '_' \
                + data.data_df.index[0].strftime('%Y%m%d%H%M%S')                
    else: 
        name = data.probe + '_' + data.name
        
    # Create the directory where the data will be saved
    if not os.path.exists(path+name): 
        os.mkdir(path+name)
    
    ''' 
    Three files will be created: 
        1. A file containing information about the data object such as start 
            time, end time, min bin and max bin. 
        2. A file containing the bin_df
        3. A file containing the data_df
    '''
    
    # Write Info to file
    f = open(path+name+'/'+name+'_'+'INFO'+'.csv','w')
    f.write('PATH='+ data.path +'\n')
    f.write('NAME='+ data.name +'\n')
    f.write('STARTTIME=' + data.t_start.strftime('%Y-%m-%d %H:%M:%S.%f'+'\n'))
    f.write('ENDTIME=' + data.t_end.strftime('%Y-%m-%d %H:%M:%S.%f'+'\n'))
    f.write('MINBIN='+str(data.min_bin)+'\n')
    f.write('MAXBIN='+str(data.max_bin)+'\n')
    f.write('PROBE='+str(data.probe+'\n'))
    f.write('TYPE='+str(data.type+'\n'))
    f.write('SOURCE='+str(data.source+'\n'))
    f.write('N='+str(data.get_mean_value('N'))+'\n')
    f.write('LWC='+str(data.get_mean_value('LWC'))+'\n')
    f.write('MVD='+str(data.get_mean_value('MVD'))+'\n')
    f.write('ED='+str(data.get_mean_value('ED'))+'\n')
    f.close()
    print('MVD='+str(data.get_mean_value('MVD')))
    
    # Write bin_df to file
    data.bin_df.to_csv(path_or_buf = path + name + '/' + name + '_'\
                       + 'BIN_DF' + '.csv', index_label = 'Index')
            
    # Write data_df to file
    data.data_df.to_csv(path_or_buf = path + name + '/' + name + '_'\
                        + 'DATA_DF' + '.csv', index_label='Index')
                 
    
    # Save the quicklook
    fig = data.plot_quicklook(log='on')
    plt.savefig(path + name + '/' + name + '_' + 'quicklook' + '.png')
    plt.close(fig)
    
    # Save housekeeping
    if data.probe in ['CDP', 'CAS', 'PIP', 'CIP']:
        fig = data.plot_housekeeping(False)
        plt.savefig(path + name + '/' + name + '_' + 'housekeeping' + '.png')
        plt.close(fig)
    
def save_psd(psd): 
    try:
        psd.path
    except: 
        psd.path = psd.Data[list(psd.Data.keys())[0]].path
        
    try:
        psd.name 
        if psd.name == '':
            psd.name = 'PSD' + '_' + psd.Data[list(psd.Data.keys())[0]].name
    except:
        psd.name = 'PSD' + '_' + psd.Data[list(psd.Data.keys())[0]].name
        
    path = psd.path
    name = psd.name
    
    
    # Create the directory where the data will be saved
    if not os.path.exists(path+'/' + name): 
        os.mkdir(path+ '/' + name)
    
    ''' 
    Three files will be created: 
        1. A file containing information about the data object such as start 
            time, end time, min bin and max bin. 
        2. A file containing the bin_df
        3. A file containing the data_df
    '''
    
    # Write Info to file
    f = open(path+'/' + name + '/'+name+'_'+'INFO'+'.csv','w')
    f.write('PATH='+ psd.path +'\n')
    f.write('NAME='+ psd.name +'\n')
    f.write('STARTTIME=' + psd.t_start.strftime('%Y-%m-%d %H:%M:%S.%f'+'\n'))
    f.write('ENDTIME=' + psd.t_end.strftime('%Y-%m-%d %H:%M:%S.%f'+'\n'))
    f.write('N='+str(psd.N)+'\n')
    f.write('LWC='+str(np.mean(psd.data_df['LWC']))+'\n')
    f.write('MVD='+str(np.mean(psd.data_df['MVD']))+'\n')
    f.write('ED='+str(np.mean(psd.data_df['ED']))+'\n')
    for probe in psd.Data.keys(): 
        f.write(probe + ' MINBIN=' + str(psd.Data[probe].min_bin) + '\n')
        f.write(probe + ' MAXBIN=' + str(psd.Data[probe].max_bin) + '\n')   
    f.close()
    
    # Write bin_df to file
    psd.bin_df.to_csv(path_or_buf = path + '/' +  name + '/' + name + '_'\
                       + 'BIN_DF' + '.csv',index_label='Index')
            
    # Write data_df to file
    psd.data_df.to_csv(path_or_buf = path + '/' +  name + '/' + name + '_'\
                        + 'DATA_DF' + '.csv',index_label='Index')
                 
    
    # Save the quicklook
    fig = psd.plot_combined()
    plt.savefig(path + '/' +  name + '/' + name + '_' \
                + 'PSD_COMBINED' + '.png')
    plt.close(fig)
    fig_cum = psd.plot_cumulative_distribution()
    plt.savefig(path + '/' +  name + '/' + name + '_' \
                + 'MASS_DISTRIBUTION' + '.png')
    plt.close(fig_cum)
    fig_port = psd.plot_pdf_distribution()
    plt.savefig(path + '/' + name + '/' + name + '_' + 'PDF' + '.png')
    plt.close(fig_port)
    
    fig_sep,ts,te = psd.plot()
    plt.savefig(path + '/' + name + '/' + name + '_' + 'PSD_SEPERATE' + '.png')
    plt.close(fig_sep)
    
