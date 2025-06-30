import numpy as np
import os
import numpy as np
import ctypes
import argparse

def covariance_simple(output_dir):

    lib = ctypes.CDLL('curved_sky_covariances_tools.so')

    config_file_str = f'{output_dir}/CosmoLike_config.txt' 
    
    with open('CosmoLike_Cells_list_demo.txt', 'r') as f:
        lines = f.readlines()
    Cell_files_str = [line.replace('output_dir', f'{output_dir}/shear_planck') for line in lines]
    
    N_lens = 0
    N_source = 2
    N_field = N_lens + N_source # 5 bins of galaxy clustering and 5 bins for cosmic shear.
    N_2point = N_field * (N_field + 1) / 2  # second term accounts for xi_minus
    N_bin = 20  # log-spaced bins
    N_data = (N_source * (N_source + 1) + N_source * N_lens + N_lens) * N_bin
    N_data_squared = N_data * N_data
    theta_min = 2.5  # in arcmin
    theta_max = 250.0  # in arcmin
    f_sky = 4143.17 * (np.pi / 180.0) ** 2 / (4.0 * np.pi)

    return_covariance = lib.return_covariance
    # return_covariance        (                           double* covariance,                           double* Gaussian_part,                                  double* Cov_SS,                                  double* Cov_SN,                                  double* Cov_NN, const char *config_file_cstr, const char *Cell_files_cstr, double theta_min, double theta_max,    int N_bin,        double f_sky, int small_bin,    int LSZ_correction_on,       int Npair_modus)
    return_covariance.argtypes = [
        ctypes.POINTER(ctypes.c_double * N_data_squared),
        ctypes.POINTER(ctypes.c_double * N_data_squared),
        ctypes.POINTER(ctypes.c_double * N_data_squared),
        ctypes.POINTER(ctypes.c_double * N_data_squared),
        ctypes.POINTER(ctypes.c_double * N_data_squared),
        ctypes.c_char_p,
        ctypes.c_char_p,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
    ]
    return_covariance.restype = None

    covariance_array = (ctypes.c_double * N_data_squared)()
    Gaussian_array = (ctypes.c_double * N_data_squared)()
    Cov_SS_array = (ctypes.c_double * N_data_squared)()
    Cov_SN_array = (ctypes.c_double * N_data_squared)()
    Cov_NN_array = (ctypes.c_double * N_data_squared)()

    Npair_modus = 0
    LSZ_correction_on = 1
    Legendre_modus = 0  # 0: curved sky & finite bin ; 1: curved sky & small bin ; 2: flat sky & finite bin

    return_covariance(
        ctypes.byref(covariance_array),
        ctypes.byref(Gaussian_array),
        ctypes.byref(Cov_SS_array),
        ctypes.byref(Cov_SN_array),
        ctypes.byref(Cov_NN_array),
        ctypes.c_char_p(config_file_str.encode("utf-8")),
        ctypes.c_char_p(Cell_files_str.encode("utf-8")),
        theta_min,
        theta_max,
        N_bin,
        f_sky,
        Legendre_modus,
        LSZ_correction_on,
        Npair_modus,
    )
    Gaussian_array = np.array(Gaussian_array)

    index = 0
    cov = np.zeros((N_data, N_data))

    i_values = np.array(range(0, N_data * N_data)).astype(int) % N_data
    i_values = i_values.astype(int)
    j_values = np.array(range(0, N_data * N_data)).astype(int) / N_data
    j_values = j_values.astype(int)

    cov[i_values, j_values] = Gaussian_array
    
    np.save(f'{output_dir}/shear_planck/covariance.npy', cov)

    return


def main():
    
    parser = argparse.ArgumentParser(description='Run covariance calculation')
    parser.add_argument('output_dir', type=str, help='The output_dir to use')
    args = parser.parse_args()

    print(f'Run covariance and save it in {args.output_dir}/shear_planck')
    covariance_simple(args.output_dir)

          
if __name__ == "__main__":
    main()
    