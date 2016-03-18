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
        self.__truncateFile = 'truncate.fits'
        self.__truncateRadius = 0
        # pollutions attributes
        self.__SExOutput = ''
        self.__pollutionFile = 'result.txt'
        self.__pollutionDataFrame = None
        self.__pollutionTreatedData = []
        self.__pollutionTreatedFile = 'treated.fits'
        # background attributes
        self.__backgroundRMS = None
        self.__backgroundThreshold = None
        self.__backgroundMean = None
        # processed data and galaxy real attributes
        self.__galaxyData = []
        self.__galaxyFile = 'galaxy.fits'
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

    # properties
    @property
    def initial_information(self):
        return {
            'hdu': self.__hdu,
            'header': self.__header,
            'data': self.__initialData,
            'center': self.__initialCenter
        }

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

    @property
    def pollution_information(self):
        return {
            'file': self.__pollutionFile,
            'info': self.__pollutionDataFrame,
            'treated': {
                'info': self.__galaxyPollutions,
                'data': self.__pollutionTreatedData,
                'file': self.__pollutionTreatedFile
            }
        }

    @property
    def background_information(self):
        return {
            'mean': self.__backgroundMean,
            'rms': self.__backgroundRMS,
            'threshold': self.__backgroundThreshold
        }

    # data process methods
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
        self.__truncateData = np.copy(self.__initialData[
                    self.__initialCenter-self.__truncateRadius:self.__initialCenter+self.__truncateRadius,
                    self.__initialCenter-self.__truncateRadius:self.__initialCenter+self.__truncateRadius])
        if os.path.exists(self.__truncateFile):
            os.system('rm '+self.__truncateFile)
        ft.writeto(self.__truncateFile, self.__truncateData)
        self.__truncated = True
        return

    def find_pollutions(self):
        if not self.__truncated:
            raise FileExistsError('you don\'t have the truncate image')
        if os.path.exists(self.__pollutionFile):
            os.system('rm '+self.__pollutionFile)
        subprocess.call('s'+'extractor '+self.__truncateFile, shell=True, executable='/bin/bash')
        self.__SExOutput = subprocess.Popen('s'+'extractor '+self.__truncateFile,
                                            shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if not os.path.exists(self.__truncateFile):
            raise FileExistsError('your SExtractor software may thread some errors')
        self.__pollutionDataFrame = pd.DataFrame(np.vstack((np.loadtxt(self.__pollutionFile), np.zeros(6))), columns=[
            'MAG_AUTO', 'X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE'])
        bg = str(self.__SExOutput.stdout.readlines()[14]).split()
        self.__backgroundMean = float(bg[2])
        self.__backgroundRMS = float(bg[4])
        self.__backgroundThreshold = float(bg[7])
        return

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
            (pdf.MAG_AUTO.values < gs.MAG_AUTO.values*0.2)]
        self.__pollutionTreatedData = np.copy(self.__truncateData)
        pollutions = self.__galaxyPollutions.iterrows()
        size = self.__pollutionTreatedData.shape[0]
        for pollution in pollutions:
            p = pollution[1]
            x_min = max(p.X_IMAGE-p.RADIUS, 0)
            x_max = min(p.X_IMAGE+p.RADIUS, size-1)
            y_min = max(p.Y_IMAGE-p.RADIUS, 0)
            y_max = min(p.Y_IMAGE+p.RADIUS, size-1)
            ptf = self.__pollutionTreatedData[y_min:y_max, x_min:x_max]
            self.__pollutionTreatedData[y_min:y_max, x_min:x_max] = \
                np.random.normal(self.__backgroundMean, self.__backgroundRMS, (ptf.shape[0], ptf.shape[1]))
        if os.path.exists(self.__pollutionTreatedFile):
            os.system('rm '+self.__pollutionTreatedFile)
        ft.writeto(self.__pollutionTreatedFile, self.__pollutionTreatedData)
        self.__galaxyData = np.copy(self.__pollutionTreatedData[
                            gs.Y_IMAGE.values-gs.RADIUS.values:gs.Y_IMAGE.values+gs.RADIUS.values,
                            gs.X_IMAGE.values-gs.RADIUS.values:gs.X_IMAGE.values+gs.RADIUS.values])
        if os.path.exists(self.__galaxyFile):
            os.system('rm '+self.__galaxyFile)
        ft.writeto(self.__galaxyFile, self.__galaxyData)
        return

    # data visualization methods
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

    def show_treated_image(self):

        os.system('ds9 '+self.__pollutionTreatedFile)
        return

    def show_galaxy_treated_image(self):
        if not os.path.exists(self.__galaxyFile):
            raise FileExistsError('the galaxy file is not exist')
        os.system('ds9 '+self.__galaxyFile)
        return
