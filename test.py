import glob2 as glob
import astropy.io.fits as ft
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import galaxy as gl
import pandas as pd


# data = ft.open('test.fits')[0].data
# ct = data.shape[0]/2
#
# box = np.zeros(ct+1)
#
#
# for i in np.arange(data.shape[1]):
#     for j in np.arange(data.shape[0]):
#         dy = abs(i-ct)
#         dx = abs(j-ct)
#         box[max(dx, dy)]+=data[i][j]
# # s=np.sum
# # inbox = np.copy(box)
# for k in np.arange(1,ct+1):
#     box[k]+=box[k-1]
#
# circle = np.array([np.sum([
#     data[ct-i, ct-i:ct+i],
#     data[ct+i, ct-i:ct+i],
#     data[ct-i:ct+i, ct-i],
#     data[ct-i:ct+i, ct+i],
# ]) for i in np.arange(1,ct)])
# print(sum(circle)+data[ct][ct])
# print(sum(box))
# # plt.figure(1)
# plt.plot([circle[i]/box[i] for i in range(len(box)-1)])
# # plt.plot([1]*len(box))
# plt.show()

a = gl.Galaxy('/home/fmxustc/Desktop/type1cut/J075525.29+391109.9_g.fits')
a.truncate()
# a.show_truncate_image()
a.find_pollutions()
a.eliminate_pollution()
a.show_galaxy_treated_image()
# a.show_treated_image()

# data = ft.open('/home/fmxustc/Desktop/tt/J000605.59-092007.0_g.fits')[0].data
# bg = np.array(data[200:350, 200:350]).flatten()
#
#
# md = np.median(bg)
# mn = np.mean(bg)
# std = np.std(bg)
# var = np.var(bg)
# print(np.max(bg),np.min(bg),std, mn)
#
# a,b,c = plt.hist(bg, 50,normed=1, facecolor='green', alpha=0.75)
# fl = mlab.normpdf(b,mn,std)
# plt.plot(b,fl,'r--')
# plt.grid()
# plt.show()Threshold



