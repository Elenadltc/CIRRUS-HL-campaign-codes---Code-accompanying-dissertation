# Model name or Project Name

Dissertation: "Microphysical properties and interplay of natural cirrus, contrail cirrus and aerosol at different latitudes" by Elena De La Torre Castro

## Description

This repository provides codes relevant for the data analysis and key plots generation of the dissertation it accompanies: "Microphysical properties and interplay of natural cirrus, contrail cirrus and aerosol at different latitudes" by Elena De La Torre Castro
The data it handles consists of cirrus ice particle in situ measurements performed during the CIRRUS-HL campaign in June and July 2021 with the german research aircraft HALO. It includes the microphysical properties of these cirrus measurements.
It also combines this data with the backward trajectories from the air masses, in order to analyse the cloud properties according to their cloud history and formation location. Further data from aerosol in situ measurements,
simulations of aerosol and ice particle properties, nitrogen oxides (NOy) measurements and reference measurements from the aircraft systems.
The cloud particle data was obtained with a combination of cloud probes: Cloud Droplet Probe (CDP), Cloud Imaging Probe (CIPG), and Precipitation Imaging Probe (PIP). In this repository, the method to obtain a combined particle size distribution is included.
For that purpose, example files of one flight are provided. For creating final results plots, combined files of measurements and model data are provided. The required data to use this codes is given in the Data folder of this dataset.
The authors of the data, variables and basic description of the data are included in the corresponding files within the Data folder.
The repository basically containes Python scripts and jupyter notebooks. The particular explanation of each code is provided in this README in the section Structure.

## Authors or Maintainers

- Elena De La Torre Castro (ORCID: 0000-0001-9848-0900), Elena.delaTorreCastro@dlr.de / E.delaTorreCastro@tudelft.nl, Institute of Atmospheric Physics, German Aerospace Center (DLR) (also Faculty of Aerospace Engineering, Delft University of Technology)

Authors of data-objekt module: 
- Simon Kirschler (ORCID: 0000-0003-4232-8277), Simon.Kirschler@dlr.de, Institute of Atmospheric Physics, German Aerospace Center (DLR)
- Johannes Lucke (ORCID: 0000-0001-6724-864X), Johannes.Lucke@dlr.de, Institute of Atmospheric Physics, German Aerospace Center (DLR) (also Faculty of Aerospace Engineering, Delft University of Technology)
- Raphael Märkl (ORCID: 0000-0002-0451-2084), Raphael.Maerkl@dlr.de, Institute of Atmospheric Physics, German Aerospace Center (DLR)
- Elena De La Torre Castro (ORCID: 0000-0001-9848-0900), Elena.delaTorreCastro@dlr.de / E.delaTorreCastro@tudelft.nl, Institute of Atmospheric Physics, German Aerospace Center (DLR) (also Faculty of Aerospace Engineering, Delft University of Technology)

## Table of contents

- Code_ACP_article
	- Combine_probes
		- Combi.py
		- Combi_cluster.py
		- CDP_CIPPIPCombi.ipynb
		- README
	- data-objekt
		- Icons
		- Data_object.py
		- housekeeping.py
		- Plots.py
		- README
		- Utility.py
	- CIRRUS_HL_Cloud_Statistics.iypnb
	- Trajectory_Plots.ipynb
	
- Code_Extra_thesis
	- Aerosol_profiles.ipynb
	- Trajectory_INP_EMAC.ipynb
	
- export_environment.yml
- LICENSE
- README

## Requirements  

The codes have been created and run using Python version 3.7.9. All the packages required in the environment used, have been exported to the file export_environment.yml.

To install environment use the following command:

    ```
        conda env create -f export_environment.yml
    ```


## Structure

Be aware that the paths indicated in the Notebooks and codes need to be adapted to the local path to the Code and Data folders of the user!

- Code_ACP_article: This folder contains the codes used to create the results of the ACP article (De La Torre Castro et al. 2023) included as Chapter 5 of the dissertation.
Additional code for results included in Chapter 6, which has not been published, are also included here, due to the topic.

	- Combine_probes: this folder provides the necessary codes to compute the combination of the particle size distribution of the three cloud probes CDP, CIPG and PIP. 
	Refer to the README of this folder for concrete details on the code. The methods are detailed in Chapter 4 of the dissertation.
		- Combi.py
		- Combi_cluster.py
		- CDP_CIPPIPCombi.ipynb
		- README
	
	- data-objekt: this folder provides a module required by the codes of Combine_probes and CIRRUS_HL_Cloud_Statistics.ipynb. Refer to the README of the file for further details.
		- Icons
		- Data_object.py
		- housekeeping.py
		- Plots.py
		- README
		- Utility.py
		
	- CIRRUS_HL_Cloud_Statistics.iypnb

	This jupyter notebook is used to create results from Chapter 5 and Chapter 6 of the disseration and the article De La Torre Castro et al. (2023) of the evaluation of the cirrus microphysical properties with respect to latitude. 
	Data from /Data/Data_ACP_article/Flights_meteo/ and /Data/Data_ACP_article/Extinction/ and tropopause data from the following repository are required:
	Hoffmann, Lars; Spang, Reinhold, 2021, "Reanalysis Tropopause Data Repository", https://doi.org/10.26165/JUELICH-DATA/UBNGI2, Jülich DATA, V1
	All files from /Data/Data_ACP_article/Flights_meteo/ are compiled in one data set, several filters are applied (e.g. measurements below -38°C), convection flights of July, 8, and July, 13 are substracted as they are not considered in the study. The data is resampled to 2-sec.
	Data from Extinction are also loaded and combined in the data set. The data of the tropopause heights form the mention repository is also incorporated in the data set and a difference of altitude of the measurement with respect to tropopause is calculated.
	The frequency distribution of the altitude relative to tropopause (TP) is plotted. Median and mean values of the cirrus properties for high (latitude >= 60°N) and mid-latitudes (latitude < 60°N) and 25th and 75th percentiles.
	Normalized frequency distributions of the ice crystal number concentration (N), effective diameter (ED), ice water content (IWC), and extinction coefficient are plotted across latitude bins of 1°N. The frequency of events are normalized by the total frequency pro latitude bin.
	Medians per latitude bin are calculated and represented with triangles. A linear fit of the medians and the 25th and 75th percentiles for the effective diameter are calculated and plotted.

	- Trajectory_Plots.ipynb

	This jupyter notebook is used to create results from Chapter 5 and 6 of the dissertation and the article De La Torre Castro et al. (2023) and it evaluates the cirrus measurements with the help of the information from the backward trajectories,
	it defines a new classification of cirrus according to their latitude of formation, and investigates potential contrail cirrus measurements with the help of the NOy measurements.
	Data from /Data/Data_ACP_article/Trajectory_data/, /Data/Data_ACP_article/Extinction/, and /Data/Data_ACP_article/Example_files_combine_probes/Final_combined/F04/CDPCIPGPIP_Combi/ are required.
	The data from /Data/Trajectory_data/ and /Data/Extinction/ is loaded. The updraft from the trajectories are translated to IS units and an analysis of the updrafts is performed, to calculate the updraft at formation and mean and standard
	deviation of the updrafts along the trajectories. Statistical plots follow on the occurrence of in situ or liquid origin cirrus per day and the events classified as HNLED, which means, high number concentration and low effective diameter,
	which are detected with the criterion N > 0.1 cm-3 and ED < 40 um. A further analysis of the relative humidity with respect to ice is provided, including the contribution of HNLED events. The updrafts of all in situ or liquid origin cirrus
	measurements at high or mid-latitudes (>= or < 60°N) along the trajectories (from formation to measurement point) are stored in lists and plotted together with the updrafts from the measurements.
	In a further plot, the updrafts from the trajectories together with the relative humidity, NO to NOy ratio and cloud age in order to provide statistics on the HNLED events.
	An example of two trajectories of a liquid origin cirrus and in situ origin cirrus are plotted, in order to exemplify the method of finding the first point at IWC=0 and containing or not LWC in that time. A map and a matrix are plotted
	to identify the classification of the cirrus into M-M, M-H and H-H (for more details refer to De La Torre Castro et al., 2023). Further statistics of N, ED and RHice are given for the cirrus measurements according to the latitude groups.
	The extinction profiles of in situ and liquid origin cirrus for high and mid-latitudes are also plotted and finally, the particle size distributions of the same groups with two methods: mean particle size distribution and median
	over an avering of 180 seconds.
	
- Code_Extra_thesis: this folder contains the codes used to create additional results included in the dissertation, in particular, the section 5.6, and chapter 7.
	
	- Aerosol_profiles.ipynb

	This jupyter notebook is used to create results from Chapter 7 of the analysis of background aerosol concentrations. It loads the in situ measurements from Data/Data_Extra_thesis/Aerosol/CPC_data and Data/Data_Extra_thesis/Aerosol/allopc.csv and the aerosol model
	data from Data/Data_Extra_thesis/Aerosol/EMAC/S4D. Since the model output is provided along the flight tracks for all pressure levels defined in the model, the relevant data at the flight track level needs to be extracted, in order to increase
	the statistics and make it more comparable to the changes of altitude of the higher resolution of the in situ measurements, one level above and below are also considered. The data is shown in standard temperature and pressure conditions (STP),
	but in the model they are given only at ambient conditions, so this is also calculated. The observations and model data are binned in pressure bins and plotted as function of pressure. The temperature profiles are also drawn, 
	and the difference between model and observations is calculated and plotted. Three further plots of comparison between the different size modes of the observations is generated (for D>12 nm, for non-volatile particles nvPM for D>14nm, and for D > 250nm)

	- Trajectory_INP_EMAC.ipynb

	This jupyter notebook is used to create results from section 5.5 of the dissertation, which analyses the aerosol influence (ice nucleating particles, INPs) along the backward trajectories of the air masses sampled in situ.
	Data from /Data/Data_ACP_article/Trajectory_data/ and /Data/Data_Extra_thesis/EMAC/ are required. The data from /Data/Data_ACP_article/Trajectory_data/ is first loaded. The file alltraject_v2 contains all the trajectories along the flight tracks, however, only the ones corresponding
	to the cirrus measurements are needed, and they are identified by the index "n" in both the cirrus data set and the trajectory data set. Then, the files from /Data/Data_Extra_thesis/EMAC/Climatology_NoAmmSu/ are loaded. Three maps of the concentrations
	of the three types of INPs are generated, with three example trajectories of the cirrus definition grouping in M-M, M-H and H-H (for more details refer to De La Torre Castro et al., 2023): dust (DU), black carbon from aviation (BCtag) and black carbon from other sources (BC).
	The variables of the model files are stored in arrays, but an adjustment to the longitude definition of the backward trajectories is done. Here, the same process would be applicable to the Hourly_NoAmmSu data, but with 4 dimensions
	(for the time) instead of 3. Since the Hourly_NoAmmSu requires more time for the computation and this was usually done in a cluster, here the example for the climatology data is rather given, however the computed data of the hourly 
	data set is also given in the file /EMAC/combi_emac_traj_1h.parquet, which can be further used for plot generation. For the latitude, longitude and pressure of the data points in the backward trajectories, the variables of the model 
	are computed accordingly and stored in a pandas.DataFrame (this takes time). An adjustment of the units and of the blank data is done. This result should be the same as provided in /EMAC/combi_emac_traj_clim.parquet. From here onwards,
	the .parquet files already provided for the climatology and the hourly output can be used for the generation of the plots, the steps before lead to those results. In the next step, several metrics are calculated, for the INP concentrations,
	the median value along all values in the trajectory is computed, all values of the newly formed ice crystals (Nice) along the trajectories are stored in dictionaries, and the parameters from the cirrus measurements are also compiled together
	with the INP and Nice information from the model. A frequency distribution of the lifetime of the three cirrus types is plotted, and further statistics are plotted, which are shown in the section 5.5 of the dissertation.


- export_environment.yml : file to recreate the employed python environment with the required packages
- LICENSE.md: file of license for this repository
- README.md


## License

The contents of this repository are licensed under a **MIT** license (see LICENSE file).

Copyright notice:
	
	Technische Universiteit Delft hereby disclaims all copyright interest in the program
	“[name_program]” (provide one line description of the content or function) written by the
	Author(s).
	Henri Werij, Faculty of Aerospace Engineering, Technische Universiteit Delft.
	© 2024, Elena De La Torre Castro, CIRRUS-HL project




## References

De La Torre Castro, E., Jurkat-Witschas, T., Afchine, A., Grewe, V., Hahn, V., Kirschler, S., Krämer, M., Lucke, J., Spelten, N., Wernli, H., Zöger, M., and Voigt, C.: Differences in 
microphysical properties of cirrus at high and mid-latitudes, Atmos. Chem. Phys., 23, 13167–13189, https://doi.org/10.5194/acp-23-13167-2023, 2023.

Hoffmann, Lars; Spang, Reinhold, 2021, "Reanalysis Tropopause Data Repository", https://doi.org/10.26165/JUELICH-DATA/UBNGI2, Jülich DATA, V1

Wernli, H., Boettcher, M., Joos, H., Miltenberger, A. K., and Spichtinger, P.: A trajectory-based classification of ERA-Interim 
ice clouds in the region of the North Atlantic storm track, Geophys. Res. Lett., 43, 6657–6664, https://doi.org/10.1002/2016GL068922, 2016.


## Citation

DOI: 10.4121/34ee4130-d189-42ec-aa64-9a310ea2912b		


