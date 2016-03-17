import astropy.io.fits as ft
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import subprocess


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
        self.__SExOutput = ''
        self.__pollutionFile = 'result.txt'
        self.__pollutionDataFrame = None
        # background attributes
        self.__backgroundRMS = None
        self.__backgroundThreshold = None
        self.__backgroundMean = None
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
        self.__loopRatio = [self.__initialLoop[k]/sum(self.__initialLoop[:k])
                            for k in np.arange(1, self.__initialLoop.size)]
        for lr in self.__loopRatio:
            if lr < 1.5e-03:
                self.__truncateRadius = self.__loopRatio.index(lr)
                break
        if not self.__truncateRadius:
            raise AttributeError('can\'t get the radius from truncate data')
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
            raise KeyError('you don\'t have the truncate image')
        subprocess.call('sextractor '+self.__truncateFile, shell=True)
        self.__SExOutput = subprocess.Popen('sextractor '+self.__truncateFile,
                                            shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if not os.path.exists(self.__truncateFile):
            raise FileExistsError('your SExtractor software may thread some errors')
        self.__pollutionDataFrame = pd.DataFrame(np.vstack((np.loadtxt(self.__pollutionFile), np.zeros(6))), columns=[
            'MAG_AUTO', 'X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE'])
        bg = str(self.__SExOutput.stdout.readlines()[14]).split()
        self.__backgroundMean = bg[2]
        self.__backgroundRMS = bg[4]
        self.__backgroundThreshold = bg[7]

        return

    @property
    def pollution_information(self):
        return {
            'file': self.__pollutionFile,
            'info': self.__pollutionDataFrame,
            'treated': self.__galaxyPollutions
        }

    @property
    def background_information(self):
        return {
            'mean': self.__backgroundMean,
            'rms': self.__backgroundRMS,
            'threshold': self.__backgroundThreshold
        }

    def eliminate_pollution(self):
        tr = self.__truncateRadius
        pdf = self.__pollutionDataFrame
        a = pdf[abs(np.sin(pdf.THETA_IMAGE)) > abs(np.cos(pdf.THETA_IMAGE))]
        a['RADIUS'] = a.A_IMAGE*3.5*abs(np.sin(a.THETA_IMAGE))
        b = pdf[abs(np.sin(pdf.THETA_IMAGE)) <= abs(np.cos(pdf.THETA_IMAGE))]
        b['RADIUS'] = b.A_IMAGE*3.5*abs(np.cos(b.THETA_IMAGE))
        self.__pollutionDataFrame = pd.concat([a, b])
        self.__pollutionDataFrame = self.__pollutionDataFrame.sort_index()
        pdf = self.__pollutionDataFrame
        self.__galaxySeries = pdf[(abs(pdf.X_IMAGE-tr) < 10) & (abs(pdf.Y_IMAGE-tr) < 10)]
        gs = self.__galaxySeries
        if gs is None:
            raise ValueError('can\'t get the galaxy\'s info')
        elif len(gs) > 1:
            raise ValueError('it occurs two or more galaxy nuclei')
        self.__pollutionDataFrame = self.__pollutionDataFrame.drop(self.__galaxySeries.index)
        pdf = self.__pollutionDataFrame
        self.__galaxyPollutions = pdf[
            (abs(pdf.X_IMAGE-gs.X_IMAGE.values) < pdf.RADIUS+gs.RADIUS.values) &
            (abs(pdf.X_IMAGE-gs.X_IMAGE.values) < pdf.RADIUS+gs.RADIUS.values) &
            pdf.MAG_AUTO < -1]
        print(self.__galaxyPollutions)
        for pollution in self.__galaxyPollutions:
            pass
        return

    def show_loop_ratio(self):
        if os.path.exists('truncate.fits'):
            plt.title('loopRatio')
            plt.plot([1e-03]*len(self.__loopRatio))
            plt.plot(self.__loopRatio)
            plt.show()
        else:
            raise ChildProcessError('You may not run the function truncate()')
        return

    def show_truncate_image(self):
        if self.__truncated:
            os.system('ds9 '+self.__truncateFile)
        else:
            raise ChildProcessError('You may not run the function truncate()')
        return
