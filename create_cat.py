import numpy as np
import pickle
import pandas as pd
import fitsio
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
import sys

np.random.seed(1738)
sys.path.append("/global/u1/e/elegnani/sompz/")
from sompz import CellMap


def assign_balrog_zbins():

    sompz_data_dir = "/global/cfs/cdirs/des/jmyles/sompz_data/v0.50/"

    print("Running assign Balrog to zbins")

    cm = CellMap.read(f"{sompz_data_dir}data_cellmap_10e6.h5")

    # balrog_file = f"{sompz_data_dir}deep_balrog.pkl"
    # sample = pickle.load(open(balrog_file, "rb"), encoding="latin1")
    balrog_file = "/global/cfs/cdirs/des/elisa/IA_decisiontree/sompz_data/deep_balrog.csv"
    sample = pd.read_csv(balrog_file)
    print(f"\nLoaded Balrog entries: {len(sample)}")
    print(f'Unique Balrog sources: {sample["ID"].nunique()}\n')

    flux_columns = ["METACAL/flux_i", "METACAL/flux_r", "METACAL/flux_z"]
    fluxtypes = ["unsheared", "sheared_1p", "sheared_1m", "sheared_2p", "sheared_2m"]
    for fluxtype in fluxtypes:
        sample_renamed = sample.rename(
            columns={
                col.replace("METACAL", fluxtype): col
                for col in flux_columns + [col.replace("flux", "flux_ivar") for col in flux_columns]
            }
        )
        sample[f"cell_wide_{fluxtype}"] = cm.assign_wide(sample_renamed)
        print(f"Assigned wide SOM using {fluxtype} fluxes")
    print("Completed wide SOM assignment")

    sample["cell_deep"] = cm.assign_deep(sample)
    print("Completed deep SOM assignment")

    with open(f"{sompz_data_dir}tomo_bins_wide_modal_even.pkl", "rb") as f:
        tomo_bins = pickle.load(f, encoding="latin1")

    sample["bin"] = sample["cell_wide_unsheared"].apply(
        lambda idx: next((num for num, cells in tomo_bins.items() if idx in cells), -1)
    )
    print("Assigned tomographic bins")

    data_pzc = np.load(f"{sompz_data_dir}pzc.npy")
    sample["pzc"] = sample["cell_deep"].apply(lambda idx: data_pzc[idx])
    print("Assigned photometric redshift cell (pzc)")

    return sample


def merge_bagpipes():

    bagpipes_dir = "/global/cfs/cdirs/des/elisa/y3_physical_ia/bagpipes/"

    print("Running merge Bagpipes runs")

    bagpipes_lengths = []
    data_frames = []

    for i in range(1, 6):
        bagpipes = fitsio.read(f"{bagpipes_dir}bagpipes_out_BigRun{i}.fits")
        bagpipes_lengths.append(len(bagpipes))
        print(f"Length of Bagpipes catalogue {i}: {len(bagpipes)}")
        data_frames.append(pd.DataFrame(bagpipes))

    bagpipes_merged = pd.concat(data_frames, ignore_index=True)
    assert len(bagpipes_merged) == sum(bagpipes_lengths), "Merged length mismatch"

    bagpipes_extra = fitsio.read(f"{bagpipes_dir}bagpipes_out_BigRunPatch.fits")
    print("Length of Bagpipes extra COSMOS catalogue: {len(bagpipes_extra)}\n")
    df_extra = pd.DataFrame(bagpipes_extra.byteswap().newbyteorder())

    bagpipes_merged = pd.concat([bagpipes_merged, df_extra], ignore_index=True)
    bagpipes_merged.drop_duplicates(subset=["ID"], inplace=True)

    print(f"Length of merged Bagpipes catalogues: {len(bagpipes_merged)}\n")

    return bagpipes_merged


def merge_balrog_bagpipes(balrog, bagpipes):
    
    print('Running merge Balrog with Bagpipes')
    
    balrog_coords = SkyCoord(ra=balrog['RA'].values*u.degree, dec=balrog['DEC'].values*u.degree)
    bagpipes_coords = SkyCoord(ra=bagpipes['RA'].values*u.degree, dec=bagpipes['DEC'].values*u.degree)
    idx, d2d, d3d = match_coordinates_sky(bagpipes_coords, balrog_coords)

    max_sep = 0.05 * u.arcsec # Smaller than 1.5'' - I want 1! obect per single balrog source
    matched_indices = [i for i, sep in enumerate(d2d) if sep < max_sep]
    matched_balrog_indices = idx[matched_indices]

    assert len(matched_balrog_indices) <= balrog['ID'].nunique(), "More matches than unique Balrog sources"
    assert len(matched_indices) <= len(bagpipes), "More matches than Bagpipes sources"
    
    print(f'\nLength of Balrog single sources matching by coordinates with Bagpipes: {len(matched_indices)}')

    matching_balrog = balrog.iloc[matched_balrog_indices]
    
    # Add all the multiple (injected) sorces to the single sources
    matching_all = balrog[balrog['ID'].isin(matching_balrog['ID'])]

    print(f'Length of Balrog data matching with Bagpipes: {len(matching_all)}\n')

    merged = pd.merge(bagpipes, matching_all, on='ID', how='inner')

    assert len(merged) == len(matching_all), "Merged length mismatch"
    
    return merged

