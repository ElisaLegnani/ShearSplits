import numpy as np
import fitsio
import shutil

np.random.seed(1738)

def get_shape_noise(w2e2_1, w2e2_2, w, w2): 

    a1 = np.sum(w2e2_1)
    a2 = np.sum(w2e2_2)
    a = np.sum(w) ** 2
    b = np.sum(w2)

    return np.sqrt((a1 / a + a2 / a) * (a / b) / 2)


def get_density(w, w2):

    a = np.sum(w) ** 2
    b = np.sum(w2)

    area = 4143 * 60.0 * 60.0  # deg^2 to arcmin^2
    # ra, dec = get_coord(data)
    # area = get_area(ra, dec) * 60 * 60

    return a / b / area
    

def get_log_normal_shift(zbin):
    # from Oliver's config file
    log_normal_shift = [
        0.00452546001972,
        0.00885240483805,
        0.0191781225951,
        0.0328746634851,
    ]
    return log_normal_shift[zbin]


def write_nz_file(nz1, nz2, output_dir):

    nz1 = nz1 / 100
    nz2 = nz2 / 100
    # /100: it's the normalization of the original n(z)s in DESshear.fits

    dt = {"names": ["BIN1", "BIN2"], "formats": [">f8", ">f8"]}
    NZ = np.zeros(len(nz1), dtype=dt)
    NZ["BIN1"] = nz1
    NZ["BIN2"] = nz2
    # Only used by cosmosis for n(z)s and theta - no need to add density, shape_noise

    output_file = f'{output_dir}/DESshear.fits'
    shutil.copy('cosmosis_runs/DESshear.fits', output_file)
    with fitsio.FITS(output_file, 'rw') as fits:
        fits[4].write(NZ)

    return 0
    