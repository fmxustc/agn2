import glob2 as glob
import astropy.io.fits as ft
import numpy as np
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

a = gl.Galaxy('/home/fmxustc/Desktop/type1cut/J120257.82+045045.1_g.fits')
a.truncate()
a.find_pollutions()
a.eliminate_pollution()

# b = a.initial_information
# c = a.truncate_information
# print(3*np.median(b['data'])-3*np.mean(b['data']),2.5*np.median(b['data'])-1.5*np.mean(b['data']))
# print(3*np.median(c['data'])-3*np.mean(c['data']),2.5*np.median(c['data'])-1.5*np.mean(c['data']))
# a = pd.read_table('result.txt', header=None)
# a= np.loadtxt('result.txt')

