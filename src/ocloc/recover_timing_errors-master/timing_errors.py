#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 15:58:16 2021
modified June 2, 2021
@author: davidnaranjo
"""
#%%
passes=2
order=4
transition_bandwidth=0.0
attenuation=0.0
with open("params.txt", "r") as f:
    lines = f.readlines()
    xcorr_path = lines[1].strip()
    station1 = lines[2].strip()
    station2 = lines[3].strip()
    dir_name = lines[4].strip()
    nn = int(lines[5].strip())
    lf = float(lines[6].strip())
    hf = float(lines[7].strip())
    cpl_dist = float(lines[8].strip())
    ref_vel = float(lines[9].strip())
    dist_trh = float(lines[10].strip())
    snr_trh = float(lines[11].strip())
    noise_st = float(lines[12].strip())
    apr_dt_st1 = float(lines[13].strip())
    apr_dt_st2 = float(lines[14].strip())
    dt_err = float(lines[15].strip())
    resp_details = lines[16].strip()

wr_resp_details = True if resp_details == 'True' else False
nt = int(nn/2)

# Time domain tapers
import numpy as np

def init_leftcos_taper(n, tapered_window):
    """
    Cosine tapering function
    
    Parameters
    ----------
    n : int
        Length of the tapering window.
    tapered_window : numpy.ndarray
        Tapering window.
    
    Returns
    -------
    None.
    """

    pi2 = np.arctan(1) * 2
    pi = np.pi
    for i in range(n):
        tapered_window[i] = (np.cos(pi - (i-1)/(n-1)*pi) + 1)/2

def init_rightcos_taper(n, tapered_window):
    """
    Cosine tapering function
    
    Parameters
    ----------
    n : int
        Length of the tapering window.
    tapered_window : numpy.ndarray
        Tapering window.
        
    Returns
    -------
    None.
    """
    pi2 = np.arctan(1) * 2
    pi = np.pi
    for i in range(n):
        tapered_window[i] = (np.cos((i-1)/(n-1)*pi) + 1)/2


xcorr = np.empty((nn))
xcorr_filt = np.empty((nn))

# Time domain tapers
nofts = int(0.1*nn)
left_taper = np.empty((nofts))
right_taper = np.empty((nofts))

init_leftcos_taper(nofts, left_taper)
init_rightcos_taper(nofts, right_taper)


from obspy import read
try:
    st = read(xcorr_path)
    xcorr = st[0].data
    nlen = len(st[0].data)
    BEGIN = st[0].stats.starttime
    DELTA = st[0].stats.delta
    nerr = 0
except Exception as e:
    nerr = e.args[0]
    print(f'Error {nerr} reading in file: {xcorr_path}')
    raise
dt = DELTA

import os
dir_8 = 'temp/'
if not os.path.exists(dir_8):
    os.makedirs(dir_8)

v1 = 2 * ref_vel
v2 = ref_vel / 2
cfs = (hf - lf) / 2
nsmpls_per_hp = int(1/cfs/2/dt)
sig = int(1/cfs/3/dt)
sig_ext = 2 * sig
sig_int = int(sig * (dt/dt_err))
sig_int_ext = 2 * sig_int

from scipy.signal import butter, filtfilt

xcorr_filt = xcorr
xcorr_filt[:nofts] *= left_taper
xcorr_filt[-nofts:] *= right_taper

nyquist = 0.5 * 1 / dt
low = lf / nyquist
high = hf / nyquist
b, a = butter(order, [low, high], btype='band')
xcorr_filt = filtfilt(b, a, xcorr_filt, axis=-1)

wl_trh_char = "{:.2f}".format(dist_trh)

min_wl = ref_vel/hf
if cpl_dist/min_wl < dist_trh:
    raise Exception("Station couple does not exceed minimum separation.")

snr_trh_char = "{:.1f}".format(snr_trh)

apr_dt_cpl = 2*apr_dt_st1 - 2*apr_dt_st2


def SNR_with_shift(cc, n3, dist, apr_dt, d1, noise_st, v1, v2):
    """
    Calculate SNR of the cross-correlation function.
    
    Parameters
    ----------
    cc : numpy.ndarray
        Cross-correlation function.
    n3 : int
        Length of the cross-correlation function.
    dist : float
        Distance between the two stations.
    apr_dt : float
        Approximate arrival time difference between the two stations.
    d1 : float
        Sampling interval of the cross-correlation function.
    noise_st : float
        Noise start time.
    v1 : float
        Fast velocity.
    v2 : float
        Slow velocity.
    c_snr : float
        Causal SNR.
    ac_snr : float
        Acausal SNR.
        
    Returns
    -------
    None.
    """
    a1r=dist/v1+apr_dt
    a2r=dist/v2+apr_dt               # Ending time of causal signal window 
    b1r=noise_st+apr_dt              # Starting time of causal noise window
    b2r=2.0*noise_st+apr_dt               # Ending time of causal noise window 
    
    c_a1i=n3+int(a1r/d1)        # Starting index of causal signal window 
    c_a2i=n3+int(a2r/d1)        # Ending index of causal signal window
    c_b1i=n3+int(b1r/d1)        # Starting index of causal noise window 
    c_b2i=n3+int(b2r/d1)        # Ending index of acausal noise window

    a1r=dist/v1-apr_dt                # Starting time of acausal signal window 
    a2r=dist/v2-apr_dt                # Ending time of acausal signal window 
    b1r=noise_st-apr_dt               # Starting time of acausal noise window
    b2r=2.0*noise_st-apr_dt                # Ending time of acausal noise window 

    ac_a1i=n3-int(a1r/d1)        # Ending index of acausal signal window 
    ac_a2i=n3-int(a2r/d1)        # Starting index of acausal signal window 
    ac_b1i=n3-int(b1r/d1)        # Ending index of acausal noise window
    ac_b2i=n3-int(b2r/d1)        # Starting index of acausal noise window

    c_amp = np.max(np.abs(cc[c_a1i:c_a2i]))
    ac_amp = np.max(np.abs(cc[ac_a2i:ac_a1i]))
    c_sum = np.sum(cc[c_b1i:c_b2i])
    ac_sum = np.sum(cc[ac_b2i:ac_b1i])
    tsum = c_sum + ac_sum
    n = (c_b2i-c_b1i+1) + (ac_b1i-ac_b2i+1)
    mean = tsum/n
    var = 0.
    for i in range(c_b1i,c_b2i+1):
        var = var + ((cc[i]-mean)**2)
    for i in range(ac_b2i,ac_b1i+1):
        var = var + ((cc[i]-mean)**2)
    nvar = var/n
    std = np.sqrt(nvar)
    c_snr = c_amp/std
    ac_snr = ac_amp/std
    return c_snr, ac_snr

c_snr, ac_snr = SNR_with_shift(xcorr_filt, int(nt+1), cpl_dist,
                               apr_dt_cpl/2., dt, noise_st, v1, v2)

# Cycle if SNR's are not high enough at this frequncy
if c_snr < snr_trh or ac_snr < snr_trh:
    raise Exception("Station couple does not exceed SNR")

# Compute envelope
apr_dt_shift=int(apr_dt_cpl/2./dt)
c2_env=int(nt+1)+int(3*cpl_dist/v2/dt)+apr_dt_shift #causal part
ac2_env=int(nt+1)-int(3*cpl_dist/v2/dt)+apr_dt_shift # acausal part

# Determine number of peaks and troughs
no_peaks=0
no_troughs=0

for k in range(ac2_env, c2_env):
    if xcorr_filt[k]-xcorr_filt[k-1] > 0. and xcorr_filt[k+1]-xcorr_filt[k] < 0.:
        no_peaks = no_peaks+1
    if xcorr_filt[k]-xcorr_filt[k-1] < 0. and xcorr_filt[k+1]-xcorr_filt[k] > 0.:
        no_troughs = no_troughs+1

# Store peaks and troughs and their times
peaks = np.zeros(no_peaks)
troughs = np.zeros(no_troughs)
peak_times = np.zeros(no_peaks)
trough_times = np.zeros(no_troughs)
no_peaks = 0
no_troughs = 0

for k in range(ac2_env, c2_env):
    if (xcorr_filt[k] - xcorr_filt[k-1] > 0. and
            xcorr_filt[k+1] - xcorr_filt[k] < 0.):
        no_peaks = no_peaks+1
        peaks[0:no_peaks-1] = xcorr_filt[k]
        peak_times[no_peaks-1]=BEGIN+(k-1)*dt
    if (xcorr_filt[k] - xcorr_filt[k-1] < 0. and
            xcorr_filt[k+1] - xcorr_filt[k] > 0.):
        no_troughs = no_troughs+1
        troughs[0:no_troughs-1] = xcorr_filt[k]
        trough_times[no_troughs-1] = BEGIN+(k-1)*dt

# Import necessary libraries
from scipy.interpolate import CubicSpline

# Construct spline for peaks
peak_spline = CubicSpline(peak_times, peaks.T)

# Construct spline for troughs
trough_spline = CubicSpline(trough_times, troughs.T)

# Create arrays for peak and trough envelopes
peak_env = np.zeros(c2_env-ac2_env+1)
trough_env = np.zeros(c2_env-ac2_env+1)


# Iterate over range between ac2_env and c2_env
for k in range(ac2_env, c2_env + 1):
    t = BEGIN + (k-1) * dt
    
    # Check if t falls within the range of peak_times
    if peak_times[0] < t < peak_times[-1]:
        peak_env[k-ac2_env] = peak_spline(t)
        
    # Check if t falls within the range of trough_times
    if trough_times[0] < t < trough_times[-1]:
        trough_env[k-ac2_env] = trough_spline(t)

# Deallocate peak_times, trough_times, peaks, and troughs
#del peak_times, trough_times, peaks, troughs



# %%
from scipy.signal import hilbert
analytic_signal = hilbert(peaks)
amplitude_envelope = np.abs(analytic_signal)

fig, (ax0, ax1) = plt.subplots(nrows=2)
ax0.plot(peak_times, peaks, label='signal')
ax0.plot(peak_times, amplitude_envelope, label='envelope')
ax0.set_xlabel("time in seconds")
ax0.legend()
# %%
