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
        self.__boxSize = 250
        self.__initialLoop = []
        # truncate data and attributes
        self.__loopRatio = []
        self.__truncateData = []
        if os.path.exists('truncate.fits'):
            os.system('rm truncate.fits')
        if os.path.exists('result.txt'):
            os.system('rm result.txt')
        self.__truncateFile = 'truncate.fits'
        self.__truncateSize = 0
        # pollutions attributes
        self.__pollutionDataFrame = None
        # processed data and galaxy real attributes
        self.__galaxyData = None
        self.__galaxySize = 0
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
        ]) for i in np.arange(self.__boxSize)])
        self.__initialLoop[0] = self.__initialData[self.__initialCenter][self.__initialCenter]
        self.__loopRatio = [self.__initialLoop[k]/sum(self.__initialLoop[:k]) for k in np.arange(1,self.__initialLoop.size)]
        for lr in self.__loopRatio:
            if lr < 1e-03:
                self.__truncateSize = self.__loopRatio.index(lr)
                break
        if not self.__truncateSize:
            raise SystemError
        self.__truncateData = self.__initialData[
                    self.__initialCenter-self.__truncateSize:self.__initialCenter+self.__truncateSize,
                    self.__initialCenter-self.__truncateSize:self.__initialCenter+self.__truncateSize]
        ft.writeto(self.__truncateFile, self.__truncateData)
        self.__truncated = True
        return

    @property
    def truncate_information(self):
        return {
            'size': self.__boxSize,
            'loop': self.__initialLoop,
            'ratio': self.__loopRatio,
            'scale': self.__truncateSize,
            'data': self.__truncateData,
            'file': self.__truncateFile
        }

    def eliminate_pollutions(self):
        if not self.__truncated:
            print('you don\'t have the truncate image')
            return
        os.system('sextractor truncate.fits ')
        if not os.path.exists('result.txt'):
            print('your SExtractor software may thread some errors')
            return
        self.__pollutionDataFrame = pd.DataFrame(np.loadtxt('result.txt'), columns=[
            'MAG_AUTO', 'X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE'])
        return

    @property
    def pollution_information(self):
        return {
            'dataFrame': self.__pollutionDataFrame
        }

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
