import astropy.io.fits as ft
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd


class Galaxy(object):

    def __init__(self, file):
        # initial data and attributes
        self.__hdu = ft.open(file)
        self.__header = self.__hdu[0].header
        self.__initialData = np.copy(self.__hdu[0].data)
        self.__hdu.close()
        self.__initialCenter = self.__initialData.shape[0]/2
        self.__boxRadius = min(250, self.__initialCenter)
        self.__initialLoop = []
        # truncate data and attributes
        self.__loopRatio = []
        self.__truncateData = []
        if os.path.exists('truncate.fits'):
            os.system('rm truncate.fits')
        if os.path.exists('result.txt'):
            os.system('rm result.txt')
        self.__truncateFile = 'truncate.fits'
        self.__truncateRadius = 0
        # pollutions attributes
        self.__pollutionFile = 'result.txt'
        self.__pollutionDataFrame = None
        # processed data and galaxy real attributes
        self.__galaxyData = None
        self.__galaxySeries = None
        self.__galaxyRadius = 0
        self.__galaxyPollutions = None
        self.__structuralParameters = {
            'G': None,
            'M': None,
            'A': None,
            'C': None
        }
        # some flags
        self.__truncated = False

    @property
    def initial_information(self):
        return {
            'hdu': self.__hdu,
            'header': self.__header,
            'data': self.__initialData,
            'center': self.__initialCenter
        }

    def truncate(self):
        if self.__truncated:
            print('you have run this function')
            return
        self.__initialLoop = np.array([np.sum([
            self.__initialData[self.__initialCenter-i, self.__initialCenter-i:self.__initialCenter+i],
            self.__initialData[self.__initialCenter+i, self.__initialCenter-i:self.__initialCenter+i],
            self.__initialData[self.__initialCenter-i:self.__initialCenter+i, self.__initialCenter-i],
            self.__initialData[self.__initialCenter-i:self.__initialCenter+i, self.__initialCenter+i],
        ]) for i in np.arange(self.__boxRadius)])
        self.__initialLoop[0] = self.__initialData[self.__initialCenter][self.__initialCenter]
        self.__loopRatio = [self.__initialLoop[k]/sum(self.__initialLoop[:k]) for k in np.arange(1,self.__initialLoop.size)]
        for lr in self.__loopRatio:
            if lr < 1e-03:
                self.__truncateRadius = self.__loopRatio.index(lr)
                break
        if not self.__truncateRadius:
            raise SystemError
        self.__truncateData = self.__initialData[
                    self.__initialCenter-self.__truncateRadius:self.__initialCenter+self.__truncateRadius,
                    self.__initialCenter-self.__truncateRadius:self.__initialCenter+self.__truncateRadius]
        ft.writeto(self.__truncateFile, self.__truncateData)
        self.__truncated = True
        return

    @property
    def truncate_information(self):
        return {
            'size': self.__boxRadius,
            'loop': self.__initialLoop,
            'ratio': self.__loopRatio,
            'scale': self.__truncateRadius,
            'data': self.__truncateData,
            'file': self.__truncateFile
        }

    def find_pollutions(self):
        if not self.__truncated:
            print('you don\'t have the truncate image')
            return
        os.system('sextractor '+self.__truncateFile)
        if not os.path.exists(self.__truncateFile):
            print('your SExtractor software may thread some errors')
            return
        self.__pollutionDataFrame = pd.DataFrame(np.vstack((np.loadtxt(self.__pollutionFile), np.zeros(6))), columns=[
            'MAG_AUTO', 'X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE'])
        return

    @property
    def pollution_information(self):
        return {
            'file': self.__pollutionFile,
            'info': self.__pollutionDataFrame
        }

    def eliminate_pollution(self):
        tr = self.__truncateRadius
        self.__galaxySeries = self.__pollutionDataFrame[
            (self.__pollutionDataFrame.X_IMAGE < tr+10) & (self.__pollutionDataFrame.X_IMAGE > tr-10) &
            (self.__pollutionDataFrame.Y_IMAGE < tr+10) & (self.__pollutionDataFrame.Y_IMAGE > tr-10)]
        if self.__galaxySeries is None:
            print('can\'t get the galaxy\'s info')
            return
        elif len(self.__galaxySeries) > 1:
            print('it occurs two or more galaxy nuclei')
            return
        self.__pollutionDataFrame = self.__pollutionDataFrame.drop(self.__galaxySeries.index)
        gs = self.__galaxySeries
        gr = self.__galaxyRadius = max(gs['A_IMAGE'].values*3*np.sin(gs['THETA_IMAGE'].values),
                                       gs['A_IMAGE'].values*3*np.cos(gs['THETA_IMAGE'].values),
                                       gs['B_IMAGE'].values*3*np.sin(gs['THETA_IMAGE'].values),
                                       gs['B_IMAGE'].values*3*np.cos(gs['THETA_IMAGE'].values))
        self.__galaxyPollutions = self.__pollutionDataFrame[
            (abs(self.__pollutionDataFrame.X_IMAGE-gs['X_IMAGE'].values) < gr) &
            (abs(self.__pollutionDataFrame.Y_IMAGE-gs['Y_IMAGE'].values) < gr)
        ]
        print(gs['A_IMAGE']*3*np.cos(gs['THETA_IMAGE']))
        return

    def show_loop_ratio(self):
        if os.path.exists('truncate.fits'):
            plt.title('loopRatio')
            plt.plot([1e-03]*len(self.__loopRatio))
            plt.plot(self.__loopRatio)
            plt.show()
        else:
            print('You may not run the function truncate()')
        return

    def show_truncate_image(self):
        if self.__truncated:
            os.system('ds9 '+self.__truncateFile)
        else:
            print('You may not run the function truncate()')
        return
