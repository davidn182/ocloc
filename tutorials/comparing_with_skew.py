# %%
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
import time
import os
import obspy

import sys
module_path = os.path.abspath(os.path.join('../src/ocloc'))
if module_path not in sys.path:
    sys.path.append(module_path)
    from ocloc import ProcessingParameters, ClockDrift, suppress_stdout
    from ocloc import ClockDrift, read_correlation_file, trim_correlation_trace, correlations_with_parameters
# Importing the main code.


def add_start_end_times_to_stations(self, filename):
    # Read the csv file
    df = pd.read_csv(filename)

    # Iterate over each row of the DataFrame
    for idx, row in df.iterrows():
        # Find the station with the same sensor code
        try:
            station = self.get_station(row['sensor_code'])
        except:
            print(f"Station {row['sensor_code']} not found")
            continue

        if station is not None:
            # List of all time-related fields
            time_fields = [
                'gps_sync_time',
                'start_time_gps',
                'st_date',
                'rt_date',
                'deployment_time',
                'end_time_recovery',
                'stop_recording',
                'skew_s_corrected_leap_second'
            ]

            # Iterate over all time-related fields
            for time_field in time_fields:
                if row[time_field] is not np.nan:
                    try:
                        setattr(station, time_field,
                                obspy.UTCDateTime(row[time_field]))
                    except ValueError:
                        print(
                            f"Problem reading Station {row['sensor_code']}")

            # Add skew value if it's not NaN
            if row['skew_s_corrected_leap_second'] is not np.nan:
                station.skew_value_s = row['skew_s_corrected_leap_second']


def calculate_drift_based_on_skew(cd, station_code):
    """
    Calculate the skew drift per day and skew at the reference time for a given station.

    Given the GPS sync time and end recovery time of a station, and its skew value at 
    the end recovery time, this function solves a system of equations to determine the
    skew drift per day and the skew value at a reference time.

    Mathematically, it solves the following system of equations:

    \[
    \begin{align*}
    a \cdot \text{skew_drift}_s\_per\_day + b &= \text{skew_on_deployment_time} \\
    c \cdot \text{skew_drift}_s\_per\_day + d &= \text{skew_on_end_recovery_time}
    \end{align*}
    \]

    Where:
    \( a = \text{time_since_sync_to_ref_in_days} \)
    \( b = 1 \) (coefficient for skew_at_reference_time)
    \( c = \text{time_since_end_recovery_to_ref_in_days} \)
    \( d = 1 \) (coefficient for skew_at_reference_time)

    The system can be represented in matrix form \(Ax = B\) where:
    \( A = \left[\begin{array}{cc} a & b \\ c & d \end{array}\right] \)
    \( x = \left[\begin{array}{c} \text{skew_drift}_s\_per\_day \\ \text{skew}_s\_at\_reference\_time \end{array}\right] \)
    \( B = \left[\begin{array}{c} \text{skew_on_deployment_time} \\ \text{skew_on_end_recovery_time} \end{array}\right] \)

    The solution is given by \(x = A^{-1}B\)
    """
    station = cd.get_station(station_code)
    time_since_sync_to_ref_in_days = (
        station.gps_sync_time - cd.reference_time) / 86400.0
    time_since_end_recovery_to_ref_in_days = (
        station.end_time_recovery - cd.reference_time) / 86400.0

    skew_on_deployment_time = 0
    skew_on_end_recovery_time = station.skew_value_s

    # Represent the coefficients as matrices
    coefficient_matrix = np.array([[time_since_sync_to_ref_in_days, 1],
                                   [time_since_end_recovery_to_ref_in_days, 1]])
    skew_vector = np.array(
        [skew_on_deployment_time, skew_on_end_recovery_time])
    # Solve the system of equations

    skew_drift_s_per_day, skew_s_at_reference_time = np.linalg.solve(
        coefficient_matrix, skew_vector)
    station.skew_drift_s_per_day = skew_drift_s_per_day
    station.skew_s_at_reference_time = skew_s_at_reference_time
    station.skew_s_at_gps_sync_time = 0
    station.skew_s_at_end_recovery_time = station.skew_value_s


def calculate_drift_based_on_skew_all_stations(cd):
    for station in cd.stations:
        if station.needs_correction:
            calculate_drift_based_on_skew(cd, station.code)
# %%


with open("bootstrap/clockdrift_rebuttal_boots.obj", 'rb') as f:
    cd = pickle.load(f)

filename = 'metadata/iceland_metadata_extended.csv'
add_start_end_times_to_stations(cd, filename)
calculate_drift_based_on_skew_all_stations(cd)

# %%


def calculate_shift_on_date(cd, date, station_code):
    station = cd.get_station(station_code)

    dt = (date - cd.reference_time) / 86400.0
    skew_drift_s_per_day = station.skew_drift_s_per_day
    skew_s_at_ref = station.skew_s_at_reference_time
    shift_s_based_on_skew = skew_drift_s_per_day * (dt) + skew_s_at_ref

    ocloc_drift_s_per_day = station.a[-1]
    ocloc_s_at_ref = station.b[-1]
    shift_s_based_on_ocloc = ocloc_drift_s_per_day * (dt) + ocloc_s_at_ref
    return shift_s_based_on_skew, shift_s_based_on_ocloc


# %%
plt.figure(figsize=(10, 5))
for i, station in enumerate(cd.stations):
    if station.needs_correction:
        plt.scatter(station.code, station.skew_drift_s_per_day *
                    365, c='C'+str(i))
        plt.scatter(station.code, station.a[-1]
                    * 365, c='C'+str(i), marker='x')

plt.scatter([], [], c='k', label="Skew drift")
plt.scatter([], [], c='k', marker='x', label="OCLOC drift")
plt.title("Clock drift [s/year]")
# Add more y-ticks
# Add more y-ticks
y_min, y_max = plt.ylim()  # get the current y-axis limits
new_ticks = np.linspace(y_min, y_max, num=20)  # create 20 evenly spaced ticks
plt.yticks(new_ticks)
plt.grid()
plt.legend()
plt.show()
# %%
plt.figure(figsize=(10, 5))
for i, station in enumerate(cd.stations):
    if station.needs_correction:
        plt.scatter(
            station.code, station.skew_s_at_reference_time, c='C'+str(i))
        plt.scatter(station.code, station.b[-1], c='C'+str(i), marker='x')
plt.scatter([], [], c='k', label="Skew at reference time")
plt.scatter([], [], c='k', marker='x', label="OCLOC at reference time")
plt.title("Skew at reference time")
plt.ylabel("dt [s]")
plt.show()

# %%
for i, station in enumerate(cd.stations):
    if station.needs_correction:
        shift_s_based_on_skew, shift_s_based_on_ocloc = calculate_shift_on_date(
            cd, station.gps_sync_time, station.code)
        plt.scatter(
            station.code, shift_s_based_on_skew, c='C'+str(i),
            label="Based on skew")
        plt.scatter(station.code, shift_s_based_on_ocloc, c='C'+str(i), marker='x',
                    label="Based on OCLOC")
plt.title("Skew at gps sync time")
plt.show()

# %%
plt.figure(figsize=(10, 5))
for i, station in enumerate(cd.stations):
    if station.needs_correction:
        shift_s_based_on_skew, shift_s_based_on_ocloc = calculate_shift_on_date(
            cd, station.end_time_recovery, station.code)
        assert np.isclose(shift_s_based_on_skew, station.skew_value_s)
        plt.scatter(
            station.code, shift_s_based_on_skew, c='C'+str(i))
        plt.scatter(station.code, shift_s_based_on_ocloc,
                    c='C'+str(i), marker='x')
plt.scatter([], [], c='k', label="Skew")
plt.scatter([], [], c='k', marker='x', label="OCLOC")

plt.title("Skew at recovery time")
y_min, y_max = plt.ylim()  # get the current y-axis limits
# new_ticks = np.linspace(y_min, y_max, num=20)  # create 20 evenly spaced ticks
# plt.yticks(new_ticks)
plt.grid()
plt.legend()
plt.show()
# %%


def plot_correlation_beforeNafter_correction(
    self, station1_code, station2_code, iteration=-1, min_t=-50, max_t=30
):
    """
    Parameters
    ----------
    station1_code: TYPE
        DESCRIPTION.
    station2_code: TYPE
        DESCRIPTION.
    iteration: TYPE, optional
        DESCRIPTION. The default is -1.
    min_t: TYPE, optional
        DESCRIPTION. The default is -50.
    max_t: TYPE, optional
        DESCRIPTION. The default is 30.

    Returns
    -------
    None.

    """
    correlation_list = self.get_correlations_of_stationpair(
        station1_code, station2_code
    )
    station1 = self.get_station(correlation_list[0].station1_code)
    station2 = self.get_station(correlation_list[0].station2_code)
    if not station1.needs_correction:
        if not station2.needs_correction:
            raise Exception("Neither station needs correction.")

        if len(station2.a) == 0:
            raise Exception(
                station2.code + " station has no correction yet."
            )
    elif len(station1.a) == 0:
        raise Exception(station1.code + " station has no correction yet.")
    ref_t = self.reference_time

    # As some traces were processed differently we will separate the plots
    # into the groups of traces that share the same processing parameters.
    for params in self.processing_parameters:
        correlations = correlations_with_parameters(
            correlation_list, params
        )

        f, (ax1, ax2, ax3) = plt.subplots(
            3, 1, sharey=True, sharex=True, figsize=(10, 12), dpi=300)
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes

        # Create inset axes for zoomed plots on the upper right corner
        axins1 = inset_axes(ax1, width="30%", height="60%", loc=1)
        axins2 = inset_axes(ax2, width="30%", height="60%", loc=2)
        axins3 = inset_axes(ax3, width="30%", height="60%", loc=3)

        for correlation in correlations:
            average_date = correlation.average_date
            shift_ocloc = 0
            shift_skew = 0
            if station1.needs_correction:
                corr_skew_sta1, corr_ocloc_sta1 = calculate_shift_on_date(
                    self, average_date, station1.code)
                shift_ocloc -= corr_ocloc_sta1
                shift_skew -= corr_skew_sta1

            if station2.needs_correction:
                corr_skew_sta2, corr_ocloc_sta2 = calculate_shift_on_date(
                    self, average_date, station2.code)
                shift_ocloc += corr_ocloc_sta2
                shift_skew += corr_skew_sta2

            freqmin = correlation.processing_parameters.freqmin
            freqmax = correlation.processing_parameters.freqmax
            tr = read_correlation_file(correlation.file_path)
            t1, data = trim_correlation_trace(
                tr, min_t, max_t, freqmin, freqmax
            )
            ax1.plot(t1, data, label=str(
                tr.stats.average_date)[:10], alpha=0.7)
            ax2.plot(t1 + shift_ocloc, data, alpha=0.7)
            ax3.plot(t1 + shift_skew, data, alpha=0.7)

            # Plotting the same data on the insets
            axins1.plot(t1, data, alpha=0.7)
            axins2.plot(t1 + shift_ocloc, data, alpha=0.7)
            axins3.plot(t1 + shift_skew, data, alpha=0.7)

        # Setting the x-limits for the insets to zoom between -20 and -10
        # seconds
        axins1.set_xlim(-20, -10)
        axins2.set_xlim(-20, -10)
        axins3.set_xlim(-20, -10)

        # Remove tick labels for the insets
        axins1.set_xticklabels([])
        axins1.set_yticklabels([])
        axins2.set_xticklabels([])
        axins2.set_yticklabels([])
        axins3.set_xticklabels([])
        axins3.set_yticklabels([])

        # Indicate the zoom effect on the main plots
        ax1.indicate_inset_zoom(axins1, edgecolor="black")
        ax2.indicate_inset_zoom(axins2, edgecolor="black")
        ax3.indicate_inset_zoom(axins3, edgecolor="black")

        f.suptitle(
            "Correction applied before and after inversion no: "
            + str(iteration)
        )
        ax1.set_title("Before correction " + tr.stats.station_pair)
        ax2.set_title("After correction with ocloc " + tr.stats.station_pair)
        ax3.set_title("After correction with skew" + tr.stats.station_pair)
        ax2.set_xlabel("Time [s]")
        ax2.set_ylabel("Amplitudes")
        ax1.set_ylabel("Amplitudes")
        ax1.set_xlim(min_t, max_t)
        # ax1.legend(loc=2)
        # ax2.legend(loc=2)
        plt.tight_layout()
        plt.show()


plot_correlation_beforeNafter_correction(
    cd, "RET", "O01", iteration=-1, min_t=-30, max_t=50)

# %%
