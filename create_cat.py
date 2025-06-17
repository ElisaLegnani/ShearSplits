import numpy as np
import pickle
import pandas as pd
import fitsio
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
import sys

np.random.seed(1738)
sys.path.append('/global/u1/e/elegnani/sompz/')
from sompz import CellMap


def assign_balrog_zbins():

    sompz_data_dir = '/global/cfs/cdirs/des/jmyles/sompz_data/v0.50/'

    print('Running assign Balrog to zbins\n')

    cm = CellMap.read(f'{sompz_data_dir}data_cellmap_10e6.h5')

    # balrog_file = f'{sompz_data_dir}deep_balrog.pkl' # Was pickled with an older pandas version --> saved to csv
    # sample = pickle.load(open(balrog_file, 'rb'), encoding='latin1')
    balrog_file = '/global/cfs/cdirs/des/elisa/IA_decisiontree/sompz_data/deep_balrog.csv'
    sample = pd.read_csv(balrog_file)
    print(f'\nLoaded Balrog entries: {len(sample)}')
    print(f'Unique Balrog sources: {sample["ID"].nunique()}\n')

    flux_columns = ['METACAL/flux_i', 'METACAL/flux_r', 'METACAL/flux_z']
    fluxtypes = ['unsheared', 'sheared_1p', 'sheared_1m', 'sheared_2p', 'sheared_2m']
    for fluxtype in fluxtypes:
        sample_renamed = sample.rename(
            columns={
                col.replace('METACAL', fluxtype): col
                for col in flux_columns + [col.replace('flux', 'flux_ivar') for col in flux_columns]
            }
        )
        sample[f'cell_wide_{fluxtype}'] = cm.assign_wide(sample_renamed)
        print(f'Assigned wide SOM using {fluxtype} fluxes')
    print('Completed wide SOM assignment')

    sample['cell_deep'] = cm.assign_deep(sample)
    print('Completed deep SOM assignment')

    with open(f'{sompz_data_dir}tomo_bins_wide_modal_even.pkl', 'rb') as f:
        tomo_bins = pickle.load(f, encoding='latin1')

    sample['bin'] = sample['cell_wide_unsheared'].apply(
        lambda idx: next((num for num, cells in tomo_bins.items() if idx in cells), -1)
    )
    print('\nAssigned tomographic bins')

    data_pzc = np.load(f'{sompz_data_dir}pzc.npy')
    sample['pzc'] = sample['cell_deep'].apply(lambda idx: data_pzc[idx])
    print('Assigned photometric redshift cell (pzc)')

    return sample


def merge_bagpipes():

    bagpipes_dir = '/global/cfs/cdirs/des/elisa/y3_physical_ia/bagpipes/'

    print('Running merge Bagpipes runs\n')

    bagpipes_lengths = []
    data_frames = []

    for i in range(1, 6):
        bagpipes = fitsio.read(f'{bagpipes_dir}bagpipes_out_BigRun{i}.fits')
        bagpipes_lengths.append(len(bagpipes))
        print(f'Length of Bagpipes catalog {i}: {len(bagpipes)}')
        data_frames.append(pd.DataFrame(bagpipes))

    bagpipes_merged = pd.concat(data_frames, ignore_index=True)
    assert len(bagpipes_merged) == sum(bagpipes_lengths), 'Merged length mismatch'

    bagpipes_extra = fitsio.read(f'{bagpipes_dir}bagpipes_out_BigRunPatch.fits')
    print(f'Length of Bagpipes extra COSMOS catalog: {len(bagpipes_extra)}\n')
    df_extra = pd.DataFrame(bagpipes_extra.byteswap().newbyteorder())

    bagpipes_merged = pd.concat([bagpipes_merged, df_extra], ignore_index=True)
    bagpipes_merged.drop_duplicates(subset=['ID'], inplace=True)

    print(f'Length of merged Bagpipes catalogs: {len(bagpipes_merged)}')

    return bagpipes_merged


def merge_balrog_bagpipes(balrog, bagpipes):
    
    print('Running merge Balrog with Bagpipes\n')
    
    balrog_coords = SkyCoord(ra=balrog['RA'].values*u.degree, dec=balrog['DEC'].values*u.degree)
    bagpipes_coords = SkyCoord(ra=bagpipes['RA'].values*u.degree, dec=bagpipes['DEC'].values*u.degree)
    idx, d2d, d3d = match_coordinates_sky(bagpipes_coords, balrog_coords)

    max_sep = 0.05 * u.arcsec # Smaller than 1.5'' - I want 1! obect per single balrog source
    matched_indices = [i for i, sep in enumerate(d2d) if sep < max_sep]
    matched_balrog_indices = idx[matched_indices]

    assert len(matched_balrog_indices) <= balrog['ID'].nunique(), 'More matches than unique Balrog sources'
    assert len(matched_indices) <= len(bagpipes), 'More matches than Bagpipes sources'
    
    print(f'Length of Balrog single sources matching by coordinates with Bagpipes: {len(matched_indices)}')

    matching_balrog = balrog.iloc[matched_balrog_indices]
    
    # Add all the multiple (injected) sorces to the single sources
    matching_all = balrog[balrog['ID'].isin(matching_balrog['ID'])]

    print(f'Length of Balrog data matching with Bagpipes: {len(matching_all)}')

    merged = pd.merge(bagpipes, matching_all, on='ID', how='inner')

    assert len(merged) == len(matching_all), 'Merged length mismatch'
    
    return merged


def mesh_average(quantity,indexx,indexy,steps,count):
    
    m = np.zeros((steps,steps)) # Revised version, was -1 before
    np.add.at(m,(indexx,indexy),quantity)
    m /= count
    return m


def assign_loggrid(x, y, xmin, xmax, xsteps, ymin, ymax, ysteps):
    """
    Returns x and y indices of data (x,y) on a log-spaced grid that runs from [xy]min to [xy]max in [xy]steps
    """
    
    logstepx = np.log10(xmax/xmin)/xsteps
    logstepy = np.log10(ymax/ymin)/ysteps
    
    indexx = (np.log10(x/xmin)/logstepx).astype(int)
    indexy = (np.log10(y/ymin)/logstepy).astype(int)
    
    indexx = np.maximum(indexx,0)
    indexx = np.minimum(indexx, xsteps-1)
    indexy = np.maximum(indexy,0)
    indexy = np.minimum(indexy, ysteps-1)
    
    return indexx,indexy


snmin=10
snmax=300
sizemin=0.5
sizemax=5
steps=20

def apply_loggrid(x, y, grid, xmin=snmin, xmax=snmax, xsteps=steps, ymin=sizemin, ymax=sizemax, ysteps=steps):
    
    indexx,indexy = assign_loggrid(x, y, xmin, xmax, xsteps, ymin, ymax, ysteps)
    res = np.zeros(len(x))
    res = grid[indexx,indexy]
    return res


def apply_nzs_weighting(mcal_s2n_r, mcal_T_r, mcal_Tpsf):
    """
    Generates smooth gridded shear weights and response weights from the simulations
    """
    
    # Use lookup tables generated from DES Y3 data
    grids_dir = '/global/cfs/cdirs/des/elisa/IA_decisiontree/y3_shape_grids/'
    w = np.genfromtxt(f'{grids_dir}y3_shape_w_grid_03_31_20.txt')
    r = np.genfromtxt(f'{grids_dir}y3_shape_response_grid_03_31_20.txt')

    shear_weight    = apply_loggrid(mcal_s2n_r, mcal_T_r/mcal_Tpsf, w) # Applies a shear weight that was gridded (how well we measure shape)
    smooth_response = apply_loggrid(mcal_s2n_r, mcal_T_r/mcal_Tpsf, r) # Smoothes the response
    
    overlap_weight = shear_weight * smooth_response # Total weight
    print('Shear weights adjusted')
    print('Response weights adjusted')
        
    return overlap_weight


columns_to_drop = [
    'U', 'G', 'R', 'I', 'Z', 'J', 'H', 'Ks',
    'UERR', 'GERR', 'RERR', 'IERR', 'ZERR', 'JERR', 'HERR', 'KsERR',
    'Chi^2',
    'Av_low', 'Av_best', 'Av_high',
    'fwhm_low', 'fwhm_best', 'fwhm_high',
    'massf_low', 'massf_best', 'massf_high',
    'metal_low', 'metal_best', 'metal_high',
    'tmax_low', 'tmax_best', 'tmax_high',
    'z_low', 'z_best', 'z_high',
    'forme_low', 'forme_best', 'forme_high',
    'sfr_low', 'sfr_best', 'sfr_high',
    'nsfr_low', 'nsfr_best', 'nsfr_high',
    'mass__low', 'mass__best', 'mass__high',
    'tform_low', 'tform_best', 'tform_high',
    'tquen_low', 'tquen_best', 'tquen_high',
    'sfh_low', 'sfh_best', 'sfh_high',
    'photo_low', 'photo_best', 'photo_high',
    'spect_low', 'spect_best', 'spect_high',
    'uvj_low', 'uvj_best', 'uvj_high',
    'chisq_low', 'chisq_best', 'chisq_high',
    'dust__low', 'dust__best', 'dust__high',
    'BDF_FLUX_DERED_CALIB_G', 'BDF_FLUX_DERED_CALIB_I', 'BDF_FLUX_DERED_CALIB_R',
    'BDF_FLUX_DERED_CALIB_U', 'BDF_FLUX_DERED_CALIB_Z',
    'BDF_FLUX_ERR_DERED_CALIB_G', 'BDF_FLUX_ERR_DERED_CALIB_I', 'BDF_FLUX_ERR_DERED_CALIB_R',
    'BDF_FLUX_ERR_DERED_CALIB_U', 'BDF_FLUX_ERR_DERED_CALIB_Z',
    'BDF_MAG_DERED_CALIB_G', 'BDF_MAG_DERED_CALIB_I', 'BDF_MAG_DERED_CALIB_R',
    'BDF_MAG_DERED_CALIB_U', 'BDF_MAG_DERED_CALIB_Z',
    'DEC_y', 'FLAGS', 'MASK_FLAGS', 'RA_y',
    'BDF_FLUX_DERED_CALIB_H', 'BDF_FLUX_DERED_CALIB_J', 'BDF_FLUX_DERED_CALIB_K',
    'BDF_FLUX_ERR_DERED_CALIB_H', 'BDF_FLUX_ERR_DERED_CALIB_J', 'BDF_FLUX_ERR_DERED_CALIB_K',
    'BDF_MAG_DERED_CALIB_H', 'BDF_MAG_DERED_CALIB_J', 'BDF_MAG_DERED_CALIB_K',
    'DEC_jhk', 'FLAGS_jhk', 'MASK_FLAGS_jhk', 'RA_jhk',
    'FIELD'
]

