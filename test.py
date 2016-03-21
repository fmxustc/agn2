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
a = gl.Galaxy('/home/fmxustc/Desktop/type1cut/J083732.70+284218.7_r.fits')
a.truncate()
a.find_pollutions()
a.eliminate_pollution()
a.calculate_parameters()
print(a.galaxy_information)



