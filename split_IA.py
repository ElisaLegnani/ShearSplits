import numpy as np

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