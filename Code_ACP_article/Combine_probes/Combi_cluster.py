import os
import sys
import pandas as pd
import numpy as np
from Combi import Combi
sys.path.insert(0, ''+'/Code/Code_ACP_article/data-objekt/')
from Utility import (get_dlogD, get_dNdD, get_dNdlogD_mean,
                     bin_df_from_edges, column_names, get_dNdlogD,
                     get_dN, calc_ED, get_dlnD)
from Data_object import Data, load
from datetime import datetime, timedelta


def calc_IWC(data_df, bin_df):
    dN = get_dN(data_df, bin_df)
    IWC_i = np.empty(np.shape(dN))
    IWC = np.empty(np.shape(dN)[0])
    for i in range(np.shape(dN)[0]):  # time axes
        IWC_i[i] = dN[i] * (0.00528 * (bin_df.bin_mid * 100) ** 2.1) * 1e-3
        IWC[i] = np.nansum(IWC_i[i])

    return IWC, IWC_i

fl_nr = sys.argv[1]

print('Start flight:', fl_nr)

# The path for the CIPG and PIP data objects need to be adapted accordingly
pathcip = ''+'/Data/Data_ACP_article/Example_files_combine_probes/CIPG/F' + fl_nr +'/'
pathpip = ''+'/Data/Data_ACP_article/Example_files_combine_probes/PIP/F' + fl_nr +'/'
list_files = os.listdir(pathcip)
CIPG = load(pathcip + '/' + list_files[0])
list_files = os.listdir(pathpip)
PIP = load(pathpip + '/' + list_files[0])

test = Combi(pip_obj=PIP, cip_obj=CIPG, cip_in=16, cip_fin=42, pip_in=2, pip_fin=6)
edgespip_new = test.get_edges_pip_init()
bin_pip_new = bin_df_from_edges(edgespip_new)

t_interval = timedelta(seconds=10)
t_in = min(CIPG.t_start, PIP.t_start)
t_final = max(CIPG.t_end, PIP.t_end)
n = int((t_final-t_in)/t_interval)

list_join = [1, 8, 16, 24]
list_drop = [2, 9, 17, 25]

sequence = [15] + [30]* 2 + [45]*2 + [60]*1 + [75]*1 + [90]*1

dNdlogD = np.zeros((len(CIPG.data_df.loc[t_in:t_final]), len(sequence) + len(CIPG.bin_df[:16]) + len(PIP.bin_df[6:])), dtype=float)
dNdD = np.zeros((len(CIPG.data_df.loc[t_in:t_final]), len(sequence) + len(CIPG.bin_df[:16]) + len(PIP.bin_df[6:])), dtype=float)
dNdlogD_int = np.zeros((len(CIPG.data_df), len(bin_pip_new)), dtype=float)
for i in range(1, n + 1):
    t_start = t_in + (i - 1) * t_interval
    t_end = t_in + i * t_interval - timedelta(seconds=1)
    x = CIPG.bin_df.bin_mid[15:44]
    y = get_dNdlogD_mean(CIPG.data_df, CIPG.bin_df * 1e-6, t_start, t_end)[15:44] * 1e-6
    x_lines, y_lines, coefs = test.get_coefficients(xi=x, yi=y, bindf_pip_new=bin_pip_new)

    print(t_start, t_end)
    if PIP.data_df.loc[t_start:t_end].to_numpy().size != 10:
        cip_in = 16
        cip_fin = 42
        dlogDcip = get_dlogD(CIPG.bin_df)[cip_in:cip_fin]
        dNcippip = get_dNdlogD(CIPG.data_df, CIPG.bin_df*1e-6).loc[t_start:t_end] *1e-6
        dNcippip = dNcippip.to_numpy()[:, cip_in:cip_fin] * dlogDcip
        bin_pipjoint = CIPG.bin_df.loc[cip_in+1:cip_fin]
        bin_pipjoint.index = np.arange(1, len(bin_pipjoint)+1)

    else:
        dNnew, dNdlogDnew = test.get_new_dNdD(bin_pip_new, coefs, t_start, t_end)
        dNdlogD_int[(i - 1) * 10: i * 10, :] = dNdlogDnew
        dNjoint, dNdlogDjoint, bin_pipjoint = test.get_dNjoint(bin_pip_new, dNnew, dNdlogDnew, list_join, list_drop)
        dNcippip, dNdlogDcippip = test.get_mean_dNs(bin_pipjoint, dNjoint, t_start, t_end)

    bin_df_final, dNdD_final, dNdlogD_final = test.get_final_dfs(sequence, bin_pipjoint, dNcippip, t_start, t_end)
    dNdlogD[(i - 1) * 10: i * 10, :] = dNdlogD_final * 1e6
    dNdD[(i - 1) * 10: i * 10, :] = dNdD_final * 1e6

def column_names(bin_df, name1, name2):

    length = len(bin_df)

    dn_names = []
    dNdD_names = []
    for i in range(1,length + 1):
        if i < 10:
            dn_names.append(name1+'_00' + str(i))
            dNdD_names.append(name2+'_00' + str(i))
        elif 10 <= i < 100:
            dn_names.append(name1+'_0' + str(i))
            dNdD_names.append(name2+'_0' + str(i))
        elif i >= 100:
            dn_names.append(name1+'_' + str(i))
            dNdD_names.append(name2+'_' + str(i))

    return dn_names, dNdD_names

dN = np.zeros(np.shape(dNdD))
for i in range(len(dNdD[0])):
    dN[:, i] = dNdD[:, i] * np.round(bin_df_final.bin_width[i+1], 1)*1e-6
N = dN.sum(axis=1)

dNdlogD_names, dNdD_names = column_names(bin_df_final, 'dNdlogD', 'dNdD')
dNdD_df = pd.DataFrame(dNdD, columns=dNdD_names, index = CIPG.data_df.loc[t_in:t_final].index)

ED = calc_ED(dNdD_df, bin_df_final*1e-6)
IWC, IWC_i = calc_IWC(dNdD_df, bin_df_final*1e-6)

data_df = dNdD_df
data_df['N'] = N
data_df['ED'] = ED
data_df['Applied PAS (m/s)'] = CIPG.data_df['Applied PAS (m/s)'].loc[t_in:t_final]
data_df['IWC'] = IWC
data_df['Cloud_flag'] = data_df['IWC'].apply(lambda x: 1 if x > 0 else 0)
data_df['MVD'] = 0
data_df['LWC'] = 0

CIPGPIP = Data()
CIPGPIP.data_df = data_df
CIPGPIP.bin_df = bin_df_final
CIPGPIP.t_start = data_df.index[0]
CIPGPIP.t_end = data_df.index[-1]
CIPGPIP.type = 'OAP'
CIPGPIP.probe = 'CIPGPIP'

os.mkdir(''+'/Data/Data_ACP_article/Example_files_combine_probes/Combi/F' + fl_nr)
CIPGPIP.max_bin = CIPGPIP.bin_df.index[-1]
CIPGPIP.min_bin = CIPGPIP.bin_df.index[0]
CIPGPIP.path = ''+'/Data/Data_ACP_article/Example_files_combine_probes/Combi/F'+fl_nr
CIPGPIP.save()