#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 15:58:16 2021
modified June 2, 2021
@author: davidnaranjo
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle


def load_clockdrift(snr_trh, dist_trh,
                    dir_ClockDrifts="./precomputed_clockdrift_objects"):
    """
    This function loads the clock drift object.

    Parameters
    ----------
    snr_trh : float
        The SNR threshold.
    dist_trh : float
        The distance threshold.
    dir_ClockDrifts : str
        The directory where the clock drift objects are stored.

    Returns
    -------
    cd : ClockDrift
        The clock drift object.

    """
    iteration = "snr_trh_" + str(snr_trh) + "__dist_trh_" + str(dist_trh)
    file_name = "ClockDrift_" + iteration + ".obj"
    path_2_clockdrift = os.path.join(dir_ClockDrifts, file_name)
    with open(path_2_clockdrift, 'rb') as f:
        cd = pickle.load(f) 
    return cd


# function to count the non nan correlation values
def count_non_nan_corr(cd):
    """
    This function counts the non nan correlation values.

    Parameters
    ----------
    cd : ClockDrift
        The clock drift object.

    Returns
    -------
    count : int
        The number of non nan correlation values.

    """
    count = 0
    for corr in cd.correlations:
        if not np.isnan(corr.t_app[-1]):
            count += 1
    return count


def remove_nan_corr(cd):
    """
    This function removes the correlation in cd.correlations that have a
    nan value in t_app[-1].
    This means, that the correlation did not meet the SNR and distance 
    thresholds, or that the program could not compute the 
    $t_{i, j, k}^{\mathrm{(+, app)}} + t_{i, j, k}^{\mathrm{(-, app)}}$.

    Parameters
    ----------
    cd : ClockDrift
        The clock drift object.

    Returns
    -------
    cd : ClockDrift
        The clock drift object.

    """
    cd.correlations = [corr for corr in cd.correlations if not np.isnan(corr.t_app[-1])]
    return cd

def get_last_a_value(self):
    """
    This function creates a dictionary where the keys are the station names
    and the values are the a[-1] values. Station names can be taken from
    cd.stations[i].code. Only stations that .needs_correction == True
    should be included in the dictionary.

    Parameters
    ----------
    cd : ClockDrift
        The clock drift object.

    Returns
    -------
    last_a_value : dict
        The dictionary with the a[-1] values as a list.

    """
    last_a_value = {}
    for station in self.stations:
        if station.needs_correction:
            if station.included_in_last_inversion:
                last_a_value[station.code] = [station.a[-1]]
            else:
                last_a_value[station.code] = [np.nan]
    return last_a_value

def get_last_b_value(self):
    """
    This function creates a dictionary where the keys are the station names
    and the values are the b[-1] values. Station names can be taken from
    cd.stations[i].code. Only stations that .needs_correction == True
    should be included in the dictionary.

    Parameters
    ----------
    cd : ClockDrift
        The clock drift object.

    Returns
    -------
    last_b_value : dict
        The dictionary with the b[-1] values as a list.

    """
    last_b_value = {}
    for station in self.stations:
        if station.needs_correction:
            if station.included_in_last_inversion:
                last_b_value[station.code] = [station.b[-1]]
            else:
                last_b_value[station.code] = [np.nan]
    return last_b_value

def combine_dicts(dict1, dict2):
    """
    This function combines two dictionaries with the same keys, where
    the values are lists, and returns a dictionary with the same keys
    and the values are the concatenation of the values of the two
    dictionaries.

    Parameters
    ----------
    dict1 : dict
        The first dictionary.
    dict2 : dict
        The second dictionary.

    Returns
    -------
    combined_dict : dict
        The combined dictionary.

    """
    combined_dict = {}
    for key in dict1.keys():
        combined_dict[key] = dict1[key] + dict2[key]
    return combined_dict


def plot_hist(data, title, xlabel, ylabel, savefig=False, filename=None):
    """
    This function plots a histogram from a list with the 5% and 95% quantiles 
    and the mean, and standard deviation.

    Parameters
    ----------
    data : list
        The list with the data.
    title : str
        The title of the plot.
    xlabel : str
        The label of the x-axis.
    ylabel : str
        The label of the y-axis.
    savefig : bool, optional
        If True, the figure is saved as a png file. The default is False.
    filename : str, optional
        The filename of the saved figure. The default is None.

    Returns
    -------
    None.

    """
    # Calculate the mean and standard deviation.
    mean = np.mean(data)
    std = np.std(data)
    
    # Calculate the 5% and 95% quantiles.
    q5 = np.quantile(data, 0.05)
    q95 = np.quantile(data, 0.95)

    # Calculate the 1% and 99% quantiles for plotting.
    q1 = np.quantile(data, 0.01)
    q99 = np.quantile(data, 0.99)

    # Plot the histogram.
    plt.figure(figsize=(8, 6))
    plt.hist(data, bins=100, density=True)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.axvline(q5, color="r", label="5% quantile")
    plt.axvline(q95, color="r", label="95% quantile")
    plt.axvline(mean, color="k", label="mean")
    plt.axvline(mean + std, color="k", linestyle="--", label="std")
    plt.axvline(mean - std, color="k", linestyle="--")
    plt.xlim(q1, q99)
    plt.legend()
    if savefig:
        plt.savefig(filename, dpi=300)
    plt.show()

def plot_jackknife_results(station, a_vals, b_vals, snr_trh, dist_trh, 
                           savefig=False, filename=None):
    """
    This function plots the jackknife results for a given station.

    Parameters
    ----------
    station : str
        The station name.
    a_vals : dict
        The dictionary with the a values for each station.
    b_vals : dict
        The dictionary with the b values for each station.
    snr_trh : float
        The signal-to-noise ratio threshold.
    dist_trh : float
        The distance threshold.
    savefig : bool, optional
        If True, the figure is saved as a png file. The default is False.
    filename : str, optional
        The filename of the saved figure. The default is None.

    Returns
    -------
    None.

    """
    # Plot the a values.
    plot_hist(a_vals[station], f"Jacknife a values for {station} (SNR > {snr_trh}, "
              f"dist > {dist_trh})", "a", "Frequency", savefig=savefig, 
              filename=filename)
    
    # Plot the b values.
    plot_hist(b_vals[station], f"Jacknife b values for {station} (SNR > {snr_trh}, "
              f"dist > {dist_trh})", "b", "Frequency", savefig=savefig, 
              filename=filename)
