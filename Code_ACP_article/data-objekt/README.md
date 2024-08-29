**Read first:** only the necessary parts of the data-objekt module for running the codes of the repository 
accompanying the dissertation "Microphysical properties and interplay of natural cirrus, contrail cirrus 
and aerosol at different latitudes" by Elena De La Torre Castro has been uploaded in the respository. 
Functions for loading raw data from each instruments have been supress here, otherwise, the description
below still applies.

Authors: Simon Kirschler (ORCID: 0000-0003-4232-8277), Johannes Lucke (ORCID: 0000-0001-6724-864X), 
Raphael Märkl (ORCID: 0000-0002-0451-2084), and Elena De La Torre Castro (ORCID: 0000-0001-9848-0900)

The data-objekt is a program that facilitates the creation of quicklooks and data evaluation. It sets a standard data format that can be used for all optical and scattering instruments, thus making data from different instruments comparable. 
The data-objekt implements a class called "Data", the purpose of an object of the class Data is to store measurement data. 

The first steps when using the data-objekt are usually: 

1. Run the file Data_object.py. 
2. Create a new Data object. Example: dobj = Data() . Now an empty Data object with name "dobj" exists. 
3. Load your data into the data object. For each instrument a seperate function to load the data exists within the class Data, i.e. for loading CAS data you would use the following line of code: dobj.get_cas_data(your_data_file) . 
All data processed with SODA can be loaded using: dobj.get_soda_data(your_nc_file) . 
4. You can look at the data from your probe by calling: dobj.plot_quicklook() . 

**Attributes of a Data object**

A Data object contains the following attributes: 
- data_df, a pandas data frame that contains the number concentration per bin (dNdD), the particle counts per bin (dN), the MVD, ED, LWC and total number concentration, for each timestamp. 
- bin_df, a pandas data frame that contains the bin boundaries, bin midpoints and bin widths.
- probe, the name of the instrument.
- type, information whether the instrument is an optical array probe (OAP), a forward scattering probe (FS), or a backscatter probe (BS).
- t_start, the start time of the measurement. By default it is the first timestamp of the data frame (data_df). You can set an alternative start by calling: d_obj.t_start = your_start_time . Note that the start time has to be a datetime object. You can create a datetime object with: your_start_time = datetime(2020,12,4,15,55,0). 
- t_end, the end time of the measurement. By default it is the last timestamp of the data frame (data_df). You can set an alternative end by calling: d_obj.t_end = your_end_time .
- min_bin, the smallest bin of the instrument data that shall be used. By default it is the smallest existing bin. You can set a different lowest bin by calling: d_obj.min_bin = your_min_bin_index. Note that the number you need to provide for your_min_bin_index is the index that your desired smallest bin has in the bin_df. Setting of the min_bin does not affect calculation made within the class Data (e.g. The MVD calculation within Data still takes into account all bins. The min_bin attribute is mostly useful when combining two Data objects into a PSD.)
- max_bin, the largest bin of the instrument data that shall be used. By default it is the largest existing bin. You can set a different largest bin by calling: d_obj.max_bin = your_max_bin_index. Note that the number you need to provide for your_max_bin_index is the index that your desired largest bin has in the bin_df. Setting of the max_bin does not affect calculation made within the class Data (e.g. The MVD calculation within Data still takes into account all bins. The max_bin attribute is mostly useful when combining two Data objects into a PSD.)
- path, the file path under which the data object shall be saved. 
- name, the name under which the data object shall be saved. 

**Functions of a Data object**

The Data object implements the following functions: 
- plot_quicklook(), plots dNdD, MVD, ED, LWC and number concentration as well as a PSD against the time stamps. 
- get_mean_value(key), calculates the mean value within the interval t_start, t_end of the data frame column with the name 'key'
- save(), saves the data object under the set path and name. 
- load(path), loads a data object saved under 'path'. Usage: dobj = load(path_to_your_data_object)

