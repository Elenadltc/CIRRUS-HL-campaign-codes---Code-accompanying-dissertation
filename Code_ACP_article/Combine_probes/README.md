The codes in "Combine_probes" facilitate the combination of the particle size distributions of the 
Cloud Droplet Probe (CDP), the Cloud Imaging Probe in grayscale (CIPG) and the Precipitation Imaging 
Probe (PIP). Example files are provided in the folder Data/Data_ACP_article/Example_files_combine_probes. 
Please, be aware that the paths need to be adjusted!!

## The steps to obtain the combined size distribution are the following:

1. The data within the folders "CDP", "CIPG", and "PIP" under the folder Data/Data_ACP_article/Example_files_combine_probes
need to be downloaded.
2. The combination of the size distribution between the CIPG and PIP is obtained using the codes "Combi.py"
and "Combi_cluster.py" of this folder. Make sure that the paths in "Combi_cluster.py" are adjusted to your 
paths. The module "data-objekt" of the "Code" folder is needed, make sure to adapt the import to your path.
3. This code was designed to be run in a cluster, but it can also be run locally. The command is:
python <path_to_Code_repository>/Code/Code_ACP_article/Combine_probes/Combi_cluster.py "04"
"04" indicates the number of flight that will be processed, the given example is the 04.
4. The combination of the size distribution between the CDP and the already combined CIPGPIP is obtained
with the jupyter notebook "CDP_CIPPIPCombi.ipynb" of this folder. Make sure that the correct paths are provided.
To run the code the data of the folders "CDP", "CIPG", and "PIP" under the folder Example_files_combine_probes
are needed, as well as "CIPGPIP", which would be the result of the previous step. Here, the module "data-objekt"
is needed, make sure to adapt the path in the imports.

## Details

- CIPG and PIP combination with Combi_cluster.py and Combi.py

The method and motivation for combining of both instruments is deeply explained in the section 
4.3 Combined size distibution of the thesis that this repository accompanies. 
A mean concentration between the concentration of the bins of the CIPG and PIP is calculated between the sizes
of 247.5 and 637.5 um. Before and after these thresholds, the CIPG bins or PIP bins are taken, respectively.
For the average calculation, a polynomial function is calculated to approximate the concentration in the
corresponding PIP bins and create a new division of 15 um per bin, in order to adjust the binning to the CIPG.
In this way, the averages between CIPG and PIP' (new) bins are calculated.

Combi_cluster.py loads the data and define the thresholds for the combination and calls the class Combi of
Combi.py to with the combination definition. The attributes of this class are:
- pip_obj: it is a data object. It is loaded in Combi_cluster.py with the function "load" from data-objekt
module and requires the path to the instrument data 
(Example_files_combine_probes/PIP/F04/PIP_20210628075104)
- cip_obj: Same but for CIPG. (Example_files_combine_probes/CIPG/F04/CIPG_20210628075338)
- cip_in: first CIPG bin of the average region (16 for bin 247.5-262.5 um)
- cip_fin: last CIPG bin for the average region (42 for 637.5-652.5 um)
- pip_in: first PIP bin of the average region (2 for 250-350 um)
- pip_fin: first PIP bin of the average region (6 for bin 650-750 um)
- bin_edges: calculated from the above specified attributes. It is an array of bin edges for the CIPG 
for the average
- edges_pip: same for the PIP

Calculations for the average are included in the functions of Combi.py. The calculations continue in 
Combi_cluster.py where the bins of the new computed averaged size distribution (between 247.5 and 637.5) 
are merged in bins of 15, 30, 45, 60, 75 and 90 um of width. Before 247.5 the bins of the CIPG are directly
added. After 637.5, PIP bins are added. With the new particle size distribution, the microphysical properties
number concentration (N), effective diameter (ED), and ice water content (IWC), which are included in Utility.py
from the data-objekt module. Finally, a data object (CIPGPIP) is created with the new particle 
size distribution. The result is what is contained in the folder Example_files_combine_probes/CIPGPIP/

- CDP combination with CIPGPIP for final distribution with notebook CDP_CIPPIPCombi.ipynb

Specific details and motivation for this method are found in the section 4.3 Combined size distibution of
the thesis. In this notebook, the data objects CDP, CIPG, PIP, and CIPGPIP are loaded. CIPG and PIP are
actually only needed for visualization of the combinations. In this notebook, the CDP bins are combined
in wider bins. The final size distribution is determined by taking the size distribution of the CDP up to
37.5 um and after that, the CIPGPIP size distribution. The microphysical properties number concentration (N), 
effective diameter (ED), and ice water content (IWC), which are included in Utility.py are calculated, and
a data object (Combi) is created with the finalized distribution. The result is provided already in 
Example_files_combine_probes/Final_combined/.
The particle size distributions for different time periods can be visualized for the different data objects.
A quicklook overview of the histograms of the microphysical properties during the flight for the final
combination is provided in the notebook.
