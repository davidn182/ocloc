#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 13:34:39 2022

@author: davidnaranjo
"""
from ocloc import *

def _start_recording_time_of_correlation(station1, station2, correlation,
                                         reference_time):
    """
    It will get the start recording time as the earliest deployment from the
    OBS involved in the cross-correlation. If we do not know the time
    when the instrument started recording we will use the reference time
    given.

    Parameters
    ----------
    station1 : ocloc.Station()
        DESCRIPTION.
    station2 : ocloc.Station()
        DESCRIPTION.
    correlation : ocloc.Correlation()
        DESCRIPTION.

    Returns
    -------
    start_recording_time : obspy.UTCDateTime

    """
    start_recording_time = obspy.UTCDateTime(reference_time)
    if sta1.needs_correction and sta2.needs_correction:
        if hasattr(sta1, 'starttime') and hasattr(sta1, 'starttime'):
            start_recording_time = min[sta1.starttime, sta2.starttime]
        elif hasattr(sta1, 'starttime'):
            start_recording_time = obspy.UTCDateTime(sta1.starttime)
        elif hasattr(sta2, 'starttime'):
            start_recording_time = obspy.UTCDateTime(sta2.starttime)
    return start_recording_time

def check_snr_thr(correlation):
    """
    Function to check whether the minimum signal to noise ratio is met for
    a single correlation object. If the correlation has the minimum SNR
    the function returns True, otherwise False.
    

    Parameters
    ----------
    correlation : TYPE
        DESCRIPTION.
    snr_trh : TYPE
        DESCRIPTION.

    Returns
    -------
    bool
        DESCRIPTION.

    """
    snr_trh = correlation.processing_parameters.snr_trh
    if correlation.snr_a < snr_trh or correlation.snr_c < snr_trh:
        return False
    return True


def check_station_separation_trh(correlation):
    """
    Function to check whether the minimum station separation is met for
    a single correlation object. If the correlation has the minimum station
    separation the function returns True, otherwise False.
    
    Parameters
    ----------
    correlation : ocloc.Correlatio object
        DESCRIPTION.

    Returns
    -------
    bool
        DESCRIPTION.
    """
    
    dist_trh = correlation.processing_parameters.dist_trh
    cpl_dist = correlation.cpl_dist
    ref_vel = correlation.processing_parameters.ref_vel
    freqmax = correlation.processing_parameters.freqmax
    min_wl = ref_vel / freqmax  # Minimum wavelength separation.

    if cpl_dist / min_wl < dist_trh:
        return False
    return True

def check_xcorr_have_same_parameters(correlations):
    check_input_correlation_list(correlations)
    freqmins = [c.processing_parameters.freqmin for c in correlations]
    freqmaxs = [c.processing_parameters.freqmax for c in correlations]
    freqmin = list(set(freqmins))
    freqmax = list(set(freqmaxs))
    if len(freqmin) != 1 or len(freqmax) != 1:
        raise Exception(
            "The processing parameters are different for each" + "correlation"
        )
    freqmin = freqmin[0]
    freqmax = freqmax[0]
    if len(correlations) < 2:
        msg = "There should be at least two correlations to use this method"
        raise Exception(msg)
    
    sta1 = list(
        set([correlation.station1_code for correlation in correlations])
    )
    sta2 = list(
        set([correlation.station2_code for correlation in correlations])
    )
    if len(sta1) != 1 or len(sta2) != 1:
        msg = "The first and second station in the correlations are not the "
        msg += "same for all the correlations."
        raise Exception(msg)
    return freqmin, freqmax, sta1, sta2


def _calculate_apriori_shift_rate(correlations):
    check_xcorr_have_same_parameters(correlations)
    earliest_time = obspy.UTCDateTime("3000-12-12")
    latest_time = obspy.UTCDateTime("1000-12-12")
    for c in correlations:
        # This checks that the correlation meets the SNR threshold.
        if not check_snr_thr(correlation):
            continue
    
        if c.average_date < earliest_time:
            earliest_time = c.average_date
            earliest_correlation = c
        if c.average_date > latest_time:
            latest_time = c.average_date
            latest_correlation = c
    if earliest_time == obspy.UTCDateTime("3000-12-12"):
        return np.nan
    if latest_time == obspy.UTCDateTime("1000-12-12"):
        return np.nan
    if earliest_correlation is latest_correlation:
        return np.nan

    earliest_tr = read_correlation_file(earliest_correlation.file_path)
    earliest_tr = earliest_tr.filter(
        "bandpass", 
        freqmin=earliest_correlation.processing_parameters.freqmin,
        freqmax=earliest_correlation.processing_parameters.freqmax,
        corners=4,
        zerophase=True)
    
    latest_tr = read_correlation_file(latest_correlation.file_path)
    latest_tr = latest_tr.filter(
        "bandpass", 
        freqmin=latest_correlation.processing_parameters.freqmin,
        freqmax=latest_correlation.processing_parameters.freqmax,
        corners=4,
        zerophase=True)

    cc = correlate(earliest_tr.data, latest_tr.data, 1000)
    shift, value = xcorr_max(cc, abs_max=False)
    time_shift = shift / earliest_tr.stats.sampling_rate

    delta_t = latest_time - earliest_time
    shift_rate = time_shift / delta_t
    return shift_rate

def calculate_apriori_dt_ins(cd, correlations, plot=False, **kwargs):
    """
    Calculates de apriori estimate of given several correlation files of the
    same station pair, given that the correlation was perform in the same
    order for all the files (meaning station 1 is the same and station 2 is
    the same)

    Parameters
    ----------
    cd: ClockDrift()
        DESCRIPTION.
    correlations: list
      list of Correlations object. You can use the following function
      to retrieve all the correlations for a given station pair:
      correlations = ClockDrift.get_correlations_of_stationpair(
          station1_code,
          station2_code)
    if plot is set to tru provide a min_t and t_max to trim the correlation
    in the times you want to check

    Returns
    -------
    None.

    """
    # Check that all correlations have the same processing params.
    freqmin, freqmax, _, _ = check_xcorr_have_same_parameters(correlations)
    
    # Add the snr value as an attribute to each correlation in the list.
    _calculate_SNR(correlations)
    
    shift_rate = _calculate_apriori_shift_rate(correlations)
    
    for correlation in correlations:
        if np.isnan(shift_rate):
            correlation.apriori_dt1 = np.nan
            correlation.apriori_dt2 = np.nan
            continue
    
        if not check_snr_thr(correlation):
            correlation.apriori_dt1 = np.nan
            correlation.apriori_dt2 = np.nan
            continue
    
        if not check_station_separation_trh(correlation):
            correlation.apriori_dt1 = np.nan
            correlation.apriori_dt2 = np.nan
            continue
    
        sta1 = cd.get_station(correlation.station1_code)
        sta2 = cd.get_station(correlation.station2_code)
    
        start_recording_time = (
            _start_recording_time_of_correlation(sta1, sta2, correlation,
                                                 cd.reference_time))
        t = correlation.average_date
        dt = (t - start_recording_time) * shift_rate
        if sta1.needs_correction:
            if sta2.needs_correction:
                correlation.apriori_dt1 = -dt / 2
                correlation.apriori_dt2 = -dt / 2
            else:
                correlation.apriori_dt1 = -dt
                correlation.apriori_dt2 = 0
        elif sta2.needs_correction:
            correlation.apriori_dt1 = 0
            correlation.apriori_dt2 = -dt
        else:
            msg = "At least one of the stations should need correction"
            msg += correlation.station1_code + correlation.station2_code
            raise ValueError(msg)

def calculate_aprioridt_4_allcorrelations(self):
    """
    Function that calculates the apriori dt for all correlation files.
    Given the all the stations contained in the clock drift object, the
    function calculates all the possible station-pair combinations
    and then calculates the apriori estimate for each of the correlations.

    When using different processing parameters, we check if the apriori
    estimate is the same, if it is not then one of the two is wrong
    so we will use the value of 0 as an apriori estimate.
    Returns
    -------
    None.

    """
    print("Calculating the apriori estimates for each stationpair")
    for i, station1 in enumerate(self.stations):
        for station2 in self.stations[i + 1:]:
            sta1 = station1
            sta2 = station2
            correlations = self.get_correlations_of_stationpair(
                sta1.code, sta2.code
            )
            if len(correlations) == 0:
                continue
            for params in self.processing_parameters:
                correlations_params = correlations_with_parameters(
                    correlations, params
                )
                # If there are no corelations for station pair... skip
                if len(correlations_params) == 0:
                    continue

                # If there is only one corelation assume the first
                # estimate as nan time shift.
                if len(correlations_params) == 1:
                    correlations_params[0].apriori_dt1 = np.nan
                    correlations_params[0].apriori_dt2 = np.nan
                    continue
                # Else calculate the apriori estimate.
                calculate_apriori_dt_ins(self, correlations_params)

#%%

# Parameters for locating the files where the correlation files and station 
# information is contained.
path2data_dir = "/Users/localadmin/Dropbox/GitHub/data"
# path2data_dir = "/Users/localadmin/Dropbox/GitHub/ocloc/tutorials/correlations_O20"
# station_file = "/Users/localadmin/Dropbox/GitHub/ocloc/tutorials/station_info"
station_file = "/Users/localadmin/Dropbox/GitHub/ocloc/tutorials/metadata/station_file.txt"

reference_time = '2014-08-21T00:00:00.000000Z'

params = ProcessingParameters(
    freqmin=0.2,   # Low freq. for the bandpass filter
    freqmax=0.4,   # High freq. for the bandpass filter 
    ref_vel=2500,  # m/s
    dist_trh=2.5,  # Minimum station separation in terms of wavelength
    snr_trh=35,    # Signal-to-noise ratio threshold
    noise_st=300,  # start of the noise window.
    dt_err=0.004,  # Sampling interval needs to be multiple of this value.
    resp_details=False)

cd = ClockDrift(station_file, path2data_dir, 
                reference_time='2014-08-21T00:00:00.000000Z',
                list_of_processing_parameters=[params])

calculate_aprioridt_4_allcorrelations(cd)
cd.calculate_dt_ins()
cd.calculate_tapp_4_allcorrelations()
cd.filter_stations()
cd.build_matrices()
cd.solve_eq(method="weighted_lstsq")
#%%
for i in range(2):
    cd.calculate_tapp_4_allcorrelations()
    cd.filter_stations()
    cd.build_matrices()
    cd.solve_eq(method="weighted_lstsq")
    cd.calculate_dt_ins()
    # cd.remove_outiers(max_error=1.5)
    # cd.calculate_tapp_4_allcorrelations()
#%%

cd.plot_obs_and_pred_shifts_with_land_stations("O08")

# for sta in cd.stations:
#     if sta.code=="O08":
#         continue
#     # cd.plot_before_n_after_t_app("O08", sta.code)
#     cd.plot_before_n_after_t_app("O08", sta.code)
#%%
cd.plot_correlation_beforeNafter_correction("O08", "KEF")
cd.plot_correlation_beforeNafter_correction("O21", "KEF")

cd.plot_matrix_diff_observed_predicted()