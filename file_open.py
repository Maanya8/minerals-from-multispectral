import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt


def cuprite_data():
    # Load the .mat file
    data = sio.loadmat("Cuprite/Cuprite_data_R188.mat")


    # print(data.keys())

    cube = data['Y']  # or data['HSI'], etc.
    # print(cube.shape)  # expect something like (250, 191, 188) or (188, 250*191)


    if cube.ndim == 2:
        n_bands, n_pixels = cube.shape
        # you need the original height/width — usually documented in the repo's README
        height, width = 250, 190
        cube_3d = cube.T.reshape(height, width, n_bands, order='F')

    # plt.imshow(cube_3d[:, :, 20], cmap='gray')
    # plt.title('Band 30 - checking reshape order')
    # plt.colorbar()
    # plt.show()

    # plt.plot(cube_3d[134, 106, :], label='bright spot')  # pick coordinates from a bright blob
    # plt.plot(cube_3d[163, 21, :], label='dark ridge')     # pick coordinates from the upper-left dark area
    # plt.legend()
    # plt.show()

    return cube_3d