import numpy as np
import pickle
import pandas as pd
import fitsio
import astropy.units as units
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
