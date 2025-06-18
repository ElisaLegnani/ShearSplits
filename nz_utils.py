import numpy as np
import copy
from scipy.signal import savgol_filter
from scipy import interpolate as interp

np.random.seed(1738)

# Photo-z calibration uncertainty from Y3 Myles paper 
dz_calib = np.array([0.015, 0.011, 0.008, 0.015])

def mean_z(pz):
    
    if np.sum(pz) == 0:
        return 0.0

    indices = np.arange(len(pz))
    return np.sum(pz * indices) / np.sum(pz) / 100


def var_z(pz):
    
    sum_pz = np.sum(pz)
    if sum_pz == 0:
        return 0.0

    indices = np.arange(len(pz))
    mean = np.sum(pz * indices) / sum_pz
    variance = np.sum(pz * (indices - mean) ** 2) / sum_pz
    return variance / (100**2)


def median_z(pz):
    if np.sum(pz) == 0:
        return 0.0

    pz_norm = pz / np.sum(pz)  # Normalize
    cdf = np.cumsum(pz_norm)   # Compute the cumulative distribution

    median_index = np.searchsorted(cdf, 0.5)  # Find index where CDF crosses 0.5
    return median_index / 100



def remove_spikes_savgol(z, pz, window_length=25, polyorder=3):
    
    smoothed = savgol_filter(pz, window_length, polyorder)
    interpolator = interp.interp1d(z, smoothed, bounds_error=False, fill_value=0.0)
    return interpolator(z)


def pz_normalization(dx, pz):
    
    norm = np.sum(dx * pz)
    return pz / norm if norm != 0 else pz


def pile_up(pz, z_bin_width=0.01, cutoff_z=3.0):

    z_cutoff_idx = int(cutoff_z / z_bin_width)
    piled = copy.copy(pz)

    # Add the high-z tail to the last valid bin
    piled[z_cutoff_idx - 1] += np.sum(pz[z_cutoff_idx:])
    piled[z_cutoff_idx:] = 0.0

    # Normalize
    return piled[:z_cutoff_idx] / np.sum(piled[:z_cutoff_idx] * z_bin_width)
    