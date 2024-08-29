import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
import matplotlib.image as mpimg
sys.path.insert(0, ''+'/Code/Code_ACP_article/data-objekt')
from Utility import get_variable, get_dlogD, get_dNdD, get_dNdlogD_mean, bin_df_from_edges, column_names, get_dNdlogD
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.stats import gmean

pipres = 100*1e-6

class Combi:

    def __init__(self, pip_obj, cip_obj, cip_in, cip_fin, pip_in, pip_fin):
        # cip_in, cip_fin etc are one number less of the index of the bin_df
        self.pip_obj = pip_obj
        self.cip_obj = cip_obj
        self.cip_in = cip_in
        self.cip_fin = cip_fin
        self.pip_in = pip_in
        self.pip_fin = pip_fin
        self.bin_edges = np.concatenate((self.cip_obj.bin_df.bin_min[cip_in:cip_fin + 1],
                                         np.array([self.cip_obj.bin_df.bin_max[cip_fin + 1]])))
        self.edges_pip = np.arange(pipres / 2 * (2 * pip_in + 1), pipres * (pip_fin + 1), pipres) * 1e6

    def get_edges_pip_init(self):
        bins1 = np.concatenate((np.array([self.bin_edges[0]]), np.array([self.edges_pip[0]])))
        bins2 = np.concatenate((np.arange(self.bin_edges[1], self.bin_edges[7], 15), np.array([self.edges_pip[1]])))
        bins3 = np.concatenate((np.arange(self.bin_edges[7], self.bin_edges[14], 15), np.array([self.edges_pip[2]])))
        bins4 = np.concatenate((np.arange(self.bin_edges[14], self.bin_edges[21], 15), np.array([self.edges_pip[3]])))
        bins5 = np.concatenate((np.arange(self.bin_edges[21], self.bin_edges[27], 15), np.array([self.edges_pip[4]])))

        edgespip_new = np.concatenate((bins1, bins2, bins3, bins4, bins5))

        return edgespip_new

    def get_coefficients(self, xi, yi, bindf_pip_new):

        def function(x, a, b, c, d, e, f):
            return a * x ** 5 + b * x ** 4 + c * x ** 3 + d * x ** 2 + e * x + f

        last = 1
        coefs = np.zeros(len(bindf_pip_new))

        popt, _ = curve_fit(function, xi, yi)
        # summarize the parameter values
        a, b, c, d, e, f = popt
        x_lines = np.arange(min(xi), max(xi), 0.5)
        y_lines = function(x_lines, a, b, c, d, e, f)

        for i in range(len(self.edges_pip) - 1):

            int_tot_i, err = quad(function, self.edges_pip[i], self.edges_pip[i + 1], args=(a, b, c, d, e, f))

            for j in range(last, len(bindf_pip_new) + 1):
                lim1 = np.round(bindf_pip_new.bin_min[j], 1)
                lim2 = np.round(bindf_pip_new.bin_max[j], 1)
                # print(i, j, lim1, lim2)
                if lim1 in [0, np.round(self.edges_pip[i + 1], 1)]:
                    if (coefs < 0).any():
                        coefs[last-1: j-1] = 1 / len(coefs[last-1: j-1])
                    last = j
                    break

                else:
                    ans, err = quad(function, lim1, lim2, args=(a, b, c, d, e, f))
                    coef = ans / int_tot_i
                    #print(j, coef)
                    coefs[j-1] = coef
                if j == len(bindf_pip_new) and (coefs < 0).any():
                    coefs[last - 1: j] = 1 / len(coefs[last - 1: j])




        return x_lines, y_lines, coefs

    def get_new_dNdD(self, bindf_pip_new, coefs, t_start, t_end):

        dlogDpip = get_dlogD(self.pip_obj.bin_df)[self.pip_in:self.pip_fin]
        dNdlogD_pip = get_dNdlogD(self.pip_obj.data_df, self.pip_obj.bin_df * 1e-6).loc[t_start:t_end]
        dNdlogD_pip = dNdlogD_pip[dNdlogD_pip.columns[self.pip_in:self.pip_fin]].to_numpy() * 1e-6
        dN_i_t = dNdlogD_pip * dlogDpip
        dlogDpipnew = get_dlogD(bindf_pip_new)

        last = 1
        dNdlogDnew = np.zeros((len(dN_i_t), len(bindf_pip_new)), dtype=float)
        dNnew = np.zeros((len(dN_i_t), len(bindf_pip_new)), dtype=float)
        for i in range(len(dN_i_t[0])):

            for j in range(last, len(bindf_pip_new) + 1):
                if np.round(bindf_pip_new.bin_min[j], 1) in [0, np.round(self.edges_pip[i + 1], 1)]:
                    last = j
                    break
                else:

                    dNmod = coefs[j-1] * dN_i_t[:, i]
                    dNdlogD_i = dNmod / dlogDpipnew[j - 1]
                    dNdlogDnew[:, j - 1] = dNdlogD_i
                    dNnew[:, j - 1] = dNmod
        self.dNnewlast = dNnew[:, -1]

        return dNnew, dNdlogDnew

    def get_dNjoint(self, bindf_pip_new, dNnew, dNdlogDnew, list_join, list_drop):

        bin_pipjoint = pd.DataFrame(columns=bindf_pip_new.columns,
                                    index=np.arange(1, len(bindf_pip_new) - len(list_join), 1))
        pointer = 1
        dNjoint = np.zeros((len(dNnew), len(bin_pipjoint)), dtype=float)
        dNdlogDjoint = np.zeros((len(dNnew), len(bin_pipjoint)), dtype=float)
        for i in range(1, len(bindf_pip_new)):
            if bindf_pip_new.index[i - 1] in list_join:
                bin_min = bindf_pip_new.bin_min[i]
                bin_max = bindf_pip_new.bin_max[i + 1]
                bin_width = bin_max - bin_min
                bin_mid = (bin_max + bin_min) / 2
                bin_pipjoint.loc[pointer] = [bin_min, bin_max, bin_mid, bin_width]
                dNi = dNnew[:, i - 1] + dNnew[:, i]
                dlogD = np.log10(bin_max / bin_min)
                dNjoint[:, pointer - 1] = dNi
                dNdlogDjoint[:, pointer - 1] = dNi/ dlogD
                pointer += 1
            elif bindf_pip_new.index[i - 1] in list_drop:
                continue
            else:
                bin_pipjoint.loc[pointer] = bindf_pip_new.loc[i]
                dNjoint[:, pointer - 1] = dNnew[:, i - 1]
                dNdlogDjoint[:, pointer - 1] = dNdlogDnew[:, i - 1]
                pointer += 1

        return dNjoint, dNdlogDjoint, bin_pipjoint

    def get_mean_dNs(self, bin_pipjoint, dNjoint, t_start, t_end):

        dlogDcip = get_dlogD(self.cip_obj.bin_df)[self.cip_in:self.cip_fin]
        dNcip = get_dNdlogD(self.cip_obj.data_df, self.cip_obj.bin_df*1e-6).loc[t_start:t_end] *1e-6
        dNcip = dNcip.to_numpy()[:, self.cip_in:self.cip_fin] * dlogDcip

        dNcippip = np.zeros((len(dNcip), len(bin_pipjoint)), dtype=float)
        dNdlogDcippip = np.zeros((len(dNcip), len(bin_pipjoint)), dtype=float)
        for i in range(len(bin_pipjoint)):
            dNcippip[:, i] = (dNcip[:, i] + dNjoint[:, i]) / 2
            dNdlogDcippip[:, i] = dNcippip[:, i] / dlogDcip[i]

        return dNcippip, dNdlogDcippip

    def get_final_dfs(self, sequence, bin_pipjoint, dNcippip, t_start, t_end):

        last = 1
        bin_joint = pd.DataFrame(columns=bin_pipjoint.columns, index=np.arange(1, len(sequence) + 1))
        dN_ol = np.zeros((len(dNcippip), len(bin_joint)), dtype=float)
        dNdlogD_ol = np.zeros((len(dNcippip), len(bin_joint)), dtype=float)
        for i, item in enumerate(sequence):
            n = item // 15
            bin_in = bin_pipjoint.loc[last]
            bin_fin = bin_pipjoint.loc[last + n - 1]
            bin_joint.bin_min.loc[i + 1] = bin_in.bin_min
            bin_joint.bin_max.loc[i + 1] = bin_fin.bin_max
            bin_joint.bin_width.loc[i + 1] = item
            bin_joint.bin_mid.loc[i + 1] = (bin_in.bin_min + bin_fin.bin_max) / 2

            dN_ol[:, i] = dNcippip[:, (last - 1):(last + n - 1)].sum(axis=1)
            dNdlogD_ol[:, i] = dN_ol[:, i] / np.log10(bin_fin.bin_max / bin_in.bin_min)

            last += n

        nbins = len(self.cip_obj.bin_df[:self.cip_in]) + len(bin_joint) + len(self.pip_obj.bin_df[self.pip_fin:])
        bin_df_final = pd.DataFrame(columns=bin_joint.columns, index=np.arange(1, nbins + 1))
        dN_final = np.zeros((len(dNcippip), nbins), dtype=float)
        bin_df_final.loc[:self.cip_in] = self.cip_obj.bin_df.loc[:self.cip_in]
        bin_df_final.loc[self.cip_in + 1:(self.cip_in + (len(bin_joint)))] = bin_joint.values

        minbin = np.round(bin_joint.bin_max.loc[len(bin_joint)], 1)
        maxbin = self.pip_obj.bin_df.bin_max.loc[self.pip_fin+1]
        bin_df_final.bin_min.loc[(len(bin_joint) + self.cip_in + 1)] = minbin
        bin_df_final.bin_max.loc[(len(bin_joint) + self.cip_in + 1)] = maxbin
        bin_df_final.bin_width.loc[(len(bin_joint) + self.cip_in + 1)] = maxbin - minbin
        bin_df_final.bin_mid.loc[(len(bin_joint) + self.cip_in + 1)] = (maxbin + minbin) / 2

        bin_df_final.loc[(len(bin_joint) + self.cip_in + 2):] = self.pip_obj.bin_df[self.pip_fin+1:].values

        dlogDcip = get_dlogD(self.cip_obj.bin_df)[:self.cip_in]
        dNcip = get_dNdlogD(self.cip_obj.data_df, self.cip_obj.bin_df * 1e-6).loc[t_start:t_end].to_numpy()[:, :self.cip_in] * 1e-6 * dlogDcip

        dN_final[:, :self.cip_in] = dNcip[:, :self.cip_in]
        dN_final[:, self.cip_in:self.cip_in + len(bin_joint)] = dN_ol

        if self.pip_obj.data_df.loc[t_start:t_end].to_numpy().size == 10:

            dlogDpip = get_dlogD(self.pip_obj.bin_df)[self.pip_fin:]
            dNpip = get_dNdlogD(self.pip_obj.data_df, self.pip_obj.bin_df * 1e-6).loc[t_start:t_end].to_numpy()[:, self.pip_fin:] * 1e-6 * dlogDpip

            dN_final[:, self.cip_in + len(bin_joint)] = dNpip[:, 0] + self.dNnewlast
            dN_final[:, self.cip_in + len(bin_joint) + 1:] = dNpip[:, 1:]

        dlogD_final = get_dlogD(bin_df_final)
        dNdlogD_final = dN_final / dlogD_final

        dNdD_final = dN_final / (bin_df_final.bin_width.to_numpy() * 1e-6)

        return bin_df_final, dNdD_final, dNdlogD_final


