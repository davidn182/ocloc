
# %%
from obspy import Stream
import pandas as pd
from obspy.signal.cross_correlation import correlate_template
import seaborn as sns
from obspy.signal.cross_correlation import correlate, xcorr_max
import numpy as np
# from seismos import obspy_ext
import obspy
import os
import matplotlib.pyplot as plt
import importlib

# Function to read event and stream data


def read_event_data(event_path, stream_path):
    event = obspy.read_events(event_path)[0]
    stream = obspy.read(stream_path)
    return event, stream

# Function to preprocess a stream


def preprocess_stream(st, event_starttime, trim_start, trim_end,
                      freqmin=0.02, freqmax=0.08):
    st = st.detrend('demean').detrend('linear')
    st = st.resample(50)
    st = st.taper(max_percentage=0.05, type='cosine')
    st = st.filter('bandpass', freqmin=freqmin, freqmax=freqmax)
    st = st.trim(starttime=event_starttime + trim_start,
                 endtime=event_starttime + trim_end)
    st = st.taper(max_percentage=0.1, type='cosine', side='both')
    return st

# Function to remove specific stations from a stream


def remove_stations(st, stations_to_remove):
    st_clean = obspy.Stream(
        [tr for tr in st if tr.stats.station not in stations_to_remove])
    return st_clean

# Function to compute maximum cross-correlations and time shifts


def compute_pairwise_correlation(trace_i, trace_j, max_time_to_shift):
    """Compute the cross-correlation and time shift for a pair of traces.

    Parameters:
        trace_i, trace_j : ObsPy Trace objects
            The two traces to compare.
        max_time_to_shift : float
            Maximum time to shift in seconds.

    Returns:
        max_corr : float
            Maximum correlation value.
        time_shift : float
            Time shift in seconds.
    """
    assert trace_i.stats.sampling_rate == trace_j.stats.sampling_rate

    max_shift_in_samples = int(max_time_to_shift * trace_i.stats.sampling_rate)
    cc = correlate(trace_i.data, trace_j.data, max_shift_in_samples)
    shift, max_corr = xcorr_max(cc, abs_max=False)

    # Convert shift to seconds
    time_shift = shift / trace_i.stats.sampling_rate

    return max_corr, time_shift


def compute_correlations(st, max_time_to_shift=100):
    """Compute the maximum cross-correlations and time shifts for all pairs of traces in a stream.

    Parameters:
        st : ObsPy Stream object
            The stream containing all traces.
        max_time_to_shift : float
            Maximum time to shift in seconds.

    Returns:
        max_correlations : DataFrame
            DataFrame containing maximum correlation values.
        time_shifts : DataFrame
            DataFrame containing time shifts in seconds.
    """
    stations = [tr.stats.station for tr in st]
    max_correlations = pd.DataFrame(np.nan, index=stations, columns=stations)
    time_shifts = pd.DataFrame(np.nan, index=stations, columns=stations)

    for i, trace_i in enumerate(st):
        for j, trace_j in enumerate(st):
            if i == j:  # Skip duplicate and self-comparisons
                continue

            max_corr, time_shift = compute_pairwise_correlation(
                trace_i, trace_j, max_time_to_shift)

            max_correlations.loc[trace_i.stats.station,
                                 trace_j.stats.station] = max_corr
            time_shifts.loc[trace_i.stats.station,
                            trace_j.stats.station] = time_shift

    return max_correlations, time_shifts

# Function to plot heatmap


def plot_heatmap(dataframe, cmap='viridis', center=0, title=None):
    plt.figure(figsize=(len(dataframe) / 1.2, len(dataframe) / 1.2))
    sns.heatmap(dataframe, annot=True, cmap=cmap, center=center)
    if title:
        plt.title(title)
    plt.show()


def plot_correlation(
        trace_i, trace_j, time_shift, max_corr_value):
    """
    Plot the correlation between two traces before and after time correction.

    Parameters:
    -----------
    - trace_i, trace_j (obspy.Trace): Traces to compare.
    - time_shift (float): Time shift in seconds.
    - max_corr_value (float): Maximum correlation value.

    Returns:
    -------
    - None

    Example:
    --------
    station1 = "SARN"
    station2 = "BS19"
    plot_correlation(
        st.select(station=station1)[0],
        st.select(station=station2)[0],
        time_shifts.loc[station1, station2],
        max_correlations.loc[station1, station2])
    """
    # Make a copy of the traces and normalize them
    trace_i = trace_i.copy().normalize()
    trace_j = trace_j.copy().normalize()
    # Plot before and after time correction
    fig, axs = plt.subplots(2, 1, figsize=(10, 12))

    # Title with station pair
    fig.suptitle(
        f'Station Pair: {trace_i.stats.station} - {trace_j.stats.station}')

    # Plot before correction
    axs[0].plot(trace_i.times(), trace_i.data, 'k',
                label=f"{trace_i.stats.station} (Original)")
    axs[0].plot(trace_j.times(), trace_j.data, 'r',
                label=f"{trace_j.stats.station} (Original)")
    axs[0].set_title('Before Time Correction')
    axs[0].legend()

    # Plot after correction
    axs[1].plot(trace_i.times(), trace_i.data, 'k',
                label=f"{trace_i.stats.station} (Original)")
    axs[1].plot(trace_j.times() + time_shift, trace_j.data,
                'r', label=f"{trace_j.stats.station} (Shifted)")
    title = 'After Time Correction \n '
    title += 'Time Shift: {:.2f} s \n max corr value: {:.2f}'.format(
        time_shift, max_corr_value)
    axs[1].set_title(title)
    axs[1].legend()

    plt.tight_layout()
    plt.show()


def find_significant_shifts(df, threshold, greater=True):
    """
    Find station pairs where the time shift exceeds or is below a threshold.

    Parameters:
    - df (pd.DataFrame): DataFrame with time shifts. Rows and columns should
                         be station names.
    - threshold (float): Value to compare for significant time shifts.
    - greater (bool): If True, find time shifts > threshold.
                      If False, find time shifts < threshold.

    Returns:
    - list of tuple: List of tuples with station pairs (row, column) having
                     time shifts either exceeding or below the threshold.

    Example:
    --------
    station_pairs = find_significant_shifts(
    time_shifts, threshold=-10, greater=False)

    # And make a nice plot:
    for station1, station2 in station_pairs:
        plot_correlation(
            st.select(station=station1)[0],
            st.select(station=station2)[0],
            time_shifts.loc[station1, station2],
            max_correlations.loc[station1, station2])
    """
    if greater:
        significant_shifts = [
            (row, col) for row in df.index for col in df.columns
            if df.loc[row, col] > threshold
        ]
    else:
        significant_shifts = [
            (row, col) for row in df.index for col in df.columns
            if df.loc[row, col] < threshold
        ]

    return significant_shifts


def read_sac_files(base_dir, station_code, start_time, end_time):
    """
    load SAC files from a given directory and time range

    Parameters:
    -----------
    base_dir: str
        path to the directory containing SAC files
    station_code: str
        station code or "*" for all stations
    start_time: str
        start time in UTCDateTime format
    end_time: str
        end time in UTCDateTime format

    Returns:
    --------
    st: obspy.Stream
        Stream containing SAC files
    """
    start_time = obspy.UTCDateTime(start_time)
    end_time = obspy.UTCDateTime(end_time)
    start_epoch = int(start_time.timestamp * 1000)
    end_epoch = int(end_time.timestamp * 1000)
    julian_day = start_time.julday
    year = start_time.year
    st = obspy.Stream()

    if station_code == "*":
        station_dirs = [os.path.join(base_dir, d) for d in os.listdir(
            base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    else:
        station_dirs = [os.path.join(base_dir, station_code)]

    for station_dir in station_dirs:
        dir_path = os.path.join(station_dir, str(year), str(julian_day))
        if not os.path.exists(dir_path):
            print(f"Directory {dir_path} not found")
            continue

        for filename in os.listdir(dir_path):
            if filename.endswith(".SAC"):
                file_epoch = int(filename.split('.')[0])
                if file_epoch <= end_epoch:
                    if file_epoch + 3600000 >= start_epoch:
                        sac_file_path = os.path.join(dir_path, filename)
                        st += obspy.read(sac_file_path)

    if len(st) > 0:
        st.merge()
        st.trim(start_time, end_time)
        return st
    else:
        print("No SAC file found for the given criteria")
        return None


# Function to filter Stream by common station codes


def filter_stream_by_common_ids(st1, st2):
    common_ids = set(tr.stats.station for tr in st1).intersection(
        set(tr.stats.station for tr in st2))

    return Stream([tr for tr in st1 if tr.stats.station in common_ids]), \
        Stream([tr for tr in st2 if tr.stats.station in common_ids])


# %%
# Your main code
path = "/Users/localadmin/Dropbox/GitHub/teleseismic_test_ocloc/data/"
event1_path = os.path.join(path, "gfz2015iatp", "event_obspyck.xml")
event2_path = os.path.join(path, "gfz2015jfwy", "event_obspyck.xml")
st_event1_path = os.path.join(path, "gfz2015iatp", "stream.mseed")
st_event2_path = os.path.join(path, "gfz2015jfwy", "stream.mseed")
event1, st_event1 = read_event_data(event1_path, st_event1_path)
event2, st_event2 = read_event_data(event2_path, st_event2_path)
stations_to_remove = ['SVSR', 'BS20', 'BS22',
                      'BS17', 'BS11', 'BS08', 'BS04', 'BS03', 'SHOP']

# Analysis for Event 1
title = "Event 1: {}".format(event1.origins[0].time)
st1 = preprocess_stream(st_event1, event1.origins[0].time, 700, 760)
st1 = remove_stations(st1, stations_to_remove)
max_corr_event1, time_shifts_event1 = compute_correlations(st1)
plot_heatmap(max_corr_event1, cmap='RdBu_r',
             center=0, title=title + '\n Max Correlation value')
plot_heatmap(time_shifts_event1, cmap='RdBu_r', center=None,
             title=title + '\n Time Shifts (s)')

# Analysis for Event 2
title = "Event 2: {}".format(event2.origins[0].time)
st2 = preprocess_stream(st_event2, event2.origins[0].time, 700, 760)
st2 = remove_stations(st2, stations_to_remove)
max_corr_event2, time_shifts_event2 = compute_correlations(st2)
plot_heatmap(max_corr_event2, cmap='RdBu_r',
             center=0, title=title + '\n Max Correlation value')
plot_heatmap(time_shifts_event2, cmap='RdBu_r', center=None,
             title=title + '\n Time Shifts (s)')
# Analysis of differences between the two events
assert list(max_corr_event1.index) == list(max_corr_event2.index)
assert list(time_shifts_event1.index) == list(time_shifts_event2.index)
assert list(max_corr_event1.columns) == list(max_corr_event2.columns)
assert list(time_shifts_event1.columns) == list(time_shifts_event2.columns)

max_corr_diff = max_corr_event1 - max_corr_event2
time_shifts_diff = time_shifts_event1 - time_shifts_event2

plot_heatmap(time_shifts_diff, cmap='RdBu_r', center=0,
             title='Difference between time shift(s) matrices')

# %% ------------------- 2nd part -------------------
# Same but now with the data stored in disk
base_dir = "/Volumes/kwintsheul/reykjanes/DECIMATED_IMAGE"
st_event1_master = read_sac_files(
    base_dir, "*", event1.origins[0].time, event1.origins[0].time+3600)
st_event2_master = read_sac_files(
    base_dir, "*", event2.origins[0].time, event2.origins[0].time+3600)
# %%
# Apply the function to your Streams
st_event1, st_event2 = filter_stream_by_common_ids(
    st_event1_master, st_event2_master)

# Now you can assert that the lengths are the same
assert len(st_event1) == len(st_event2)
# %%
# stations_to_remove = ['O03', 'O04', 'O11']
st_event1 = remove_stations(st_event1, stations_to_remove)
st_event2 = remove_stations(st_event2, stations_to_remove)

# Analysis for Event 1
title = "Event 1: {}".format(event1.origins[0].time)
st1 = preprocess_stream(st_event1, event1.origins[0].time, 700, 760)
st1 = remove_stations(st1, stations_to_remove)
max_corr_event1, time_shifts_event1 = compute_correlations(st1)
plot_heatmap(max_corr_event1, cmap='RdBu_r',
             center=0, title=title + '\n Max Correlation value')
plot_heatmap(time_shifts_event1, cmap='RdBu_r', center=None,
             title=title + '\n Time Shifts (s)')
# Analysis for Event 2
title = "Event 2: {}".format(event2.origins[0].time)
st2 = preprocess_stream(st_event2, event2.origins[0].time, 700, 760)
st2 = remove_stations(st2, stations_to_remove)
max_corr_event2, time_shifts_event2 = compute_correlations(st2)
plot_heatmap(max_corr_event2, cmap='RdBu_r',
             center=0, title=title + '\n Max Correlation value')
plot_heatmap(time_shifts_event2, cmap='RdBu_r', center=None,
             title=title + '\n Time Shifts (s)')

# Analysis of differences between the two events
assert list(max_corr_event1.index) == list(max_corr_event2.index)
assert list(time_shifts_event1.index) == list(time_shifts_event2.index)
assert list(max_corr_event1.columns) == list(max_corr_event2.columns)
assert list(time_shifts_event1.columns) == list(time_shifts_event2.columns)

max_corr_diff = max_corr_event1 - max_corr_event2
time_shifts_diff = time_shifts_event1 - time_shifts_event2

plot_heatmap(time_shifts_diff, cmap='RdBu_r', center=0,
             title='Difference between time shift(s) matrices')

# %%
columns_to_exclude = [
    col for col in time_shifts_diff.columns if col.startswith('O')]
time_shift_diff_excluding_OBS = time_shifts_diff.drop(
    columns=columns_to_exclude, axis=1)

plot_heatmap(time_shift_diff_excluding_OBS, cmap='RdBu_r', center=0,
             title='Difference between time shift(s) matrices')

dt = (event2.origins[0].time - event1.origins[0].time) / 86400.0
df = pd.DataFrame()
df['mean_excluding_O'] = time_shifts_diff.drop(
    columns=columns_to_exclude, axis=1).mean(axis=1)
df['mean_excluding_O'] = df['mean_excluding_O']  # * 365 / dt
# Display the first few rows of the modified dataframe
df

# %%
station1 = "O02"
station2 = "O06"

plot_correlation(
    st1.select(station=station1)[0],
    st1.select(station=station2)[0],
    time_shifts_event1.loc[station1, station2],
    max_corr_event1.loc[station1, station2])

plot_correlation(
    st2.select(station=station1)[0],
    st2.select(station=station2)[0],
    time_shifts_event2.loc[station1, station2],
    max_corr_event2.loc[station1, station2])
# %%
