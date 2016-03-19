import glob2 as glob
import astropy.io.fits as ft
import numpy as np
import os
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import galaxy as gl
import pandas as pd
import warnings

warnings.filterwarnings('ignore')
a = gl.Galaxy('/home/fmxustc/Desktop/type1cut/J122342.81+581446.1_r.fits')
a.truncate()
a.find_pollutions()
a.eliminate_pollution()
a.calculate_parameters()
# cxx = 1.0129027e-02
# cyy = 3.6632743e-02
# cxy = 1.6219877e-02
# r = 3.5
# data = np.array(ft.open('galaxy.fits')[0].data)
# bg = np.zeros(data.size).reshape(data.shape)
# for i in np.arange(data.shape[0]):
#     for j in np.arange(data.shape[1]):
#         if cxx*((i-data.shape[0]/2)**2)+cyy*((j-data.shape[1]/2)**2)+cxy*(j-data.shape[0]/2)*(i-data.shape[1]/2) <= r**2:
#             bg[i][j] = 1
# if os.path.exists('qwe.fits'):
#     os.system('rm qwe.fits')
# ft.writeto('qwe.fits', bg)
# os.system('ds9 qwe.fits')
# for k in bg.flattern():
#     print(bg.)


