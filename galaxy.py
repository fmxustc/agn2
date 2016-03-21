import astropy.io.fits as ft
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
import json


class Galaxy(object):

    def __init__(self, file):
        # initial data and attributes
        self.__initialFile = file
        self.__hdu = ft.open(file)
        self.__header = self.__hdu[0].header
        self.__initialData = np.copy(self.__hdu[0].data)
        self.__hdu.close()
        self.__initialCentroid = self.__initialData.shape[0]/2
        self.__boxRadius = min(250, self.__initialCentroid-1)
        self.__initialLoop = None
        # truncate data and attributes
        self.__loopRatio = None
        self.__truncateData = None
        self.__truncateFile = 'truncate.fits'
        self.__truncateRadius = 0
        # pollutions attributes
        self.__SExOutput = ''
        self.__pollutionFile = 'result.txt'
        self.__pollutionDataFrame = None
        self.__pollutionTreatedData = None
        self.__pollutionTreatedFile = 'treated.fits'
        # background attributes
        self.__backgroundRMS = None
        self.__backgroundThreshold = None
        self.__backgroundMean = None
        # processed data and galaxy real attributes
        self.__galaxyData = None
        self.__galaxyFile = 'galaxy.fits'
        self.__galaxySeries = None
        self.__galaxyPollutions = None
        self.__galaxySurfaceBrightness = None
        self.__galaxyMeanSurfaceBrightness = None
        self.__galaxyStructuralParameter = {
            'G': None,
            'M': None,
            'A': None,
            'C': None
        }
        # some flags
        self.__truncated = False
        self.__found = False
        self.__treated = False
        # calculating extra things
        self.__isInGalaxy = None
        self.__galaxyMoment = None

    # properties
    @property
    def initial_information(self):
        return {
            'file': self.__initialFile,
            'hdu': self.__hdu,
            'header': self.__header,
            # 'data': self.__initialData,
            'centroid': self.__initialCentroid
        }

    @property
    def truncate_information(self):
        return {
            'size': self.__boxRadius,
            'loop': self.__initialLoop,
            'ratio': self.__loopRatio,
            'scale': self.__truncateRadius,
            # 'data': self.__truncateData,
            'file': self.__truncateFile
        }

    @property
    def pollution_information(self):
        return {
            'file': self.__pollutionFile,
            'info': self.__pollutionDataFrame,
            'treated': {
                'info': self.__galaxyPollutions,
                # 'data': self.__pollutionTreatedData,
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

    @property
    def galaxy_information(self):
        return {
            'info': self.__galaxySeries,
            'file': self.__galaxyFile,
            'parameter': self.__galaxyStructuralParameter,
            'surface_brightness': self.__galaxySurfaceBrightness
        }

    # data process methods
    def truncate(self):
        if self.__truncated:
            print('you have run this function')
            return
        self.__initialLoop = np.array([np.sum([
            self.__initialData[self.__initialCentroid-i, self.__initialCentroid-i:self.__initialCentroid+i+1],
            self.__initialData[self.__initialCentroid+i, self.__initialCentroid-i:self.__initialCentroid+i+1],
            self.__initialData[self.__initialCentroid-i:self.__initialCentroid+i+1, self.__initialCentroid-i],
            self.__initialData[self.__initialCentroid-i:self.__initialCentroid+i+1, self.__initialCentroid+i],
        ])-np.sum([
            self.__initialData[self.__initialCentroid-i][self.__initialCentroid+i],
            self.__initialData[self.__initialCentroid-i][self.__initialCentroid-i],
            self.__initialData[self.__initialCentroid+i][self.__initialCentroid+i],
            self.__initialData[self.__initialCentroid+i][self.__initialCentroid-i]
        ]) for i in np.arange(self.__boxRadius)])
        self.__initialLoop[0] = self.__initialData[self.__initialCentroid][self.__initialCentroid]
        self.__loopRatio = [self.__initialLoop[k]/sum(self.__initialLoop[:k])
                            for k in np.arange(1, self.__initialLoop.size)]
        for lr in self.__loopRatio:
            if lr < 1.5e-03:
                self.__truncateRadius = self.__loopRatio.index(lr)
                break
        if not self.__truncateRadius:
            raise AttributeError('can\'t get the radius from truncate data')
        self.__truncateData = np.copy(self.__initialData[
                    self.__initialCentroid-self.__truncateRadius:self.__initialCentroid+self.__truncateRadius,
                    self.__initialCentroid-self.__truncateRadius:self.__initialCentroid+self.__truncateRadius])
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
        self.__pollutionDataFrame = pd.DataFrame(np.vstack((np.loadtxt(self.__pollutionFile), np.zeros(9))), columns=[
            'MAG_AUTO', 'X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 'CXX_IMAGE', 'CYY_IMAGE', 'CXY_IMAGE'])
        bg = str(self.__SExOutput.stdout.readlines()[14]).split()
        self.__backgroundMean = float(bg[2])
        self.__backgroundRMS = float(bg[4])
        self.__backgroundThreshold = float(bg[7])
        self.__found = True
        return

    def eliminate_pollution(self):
        if not self.__found:
            raise KeyError('you may not use SExtractor find the pollutions')
        tr = self.__truncateRadius
        pdf = self.__pollutionDataFrame
        a = pdf[abs(np.sin(np.deg2rad(pdf.THETA_IMAGE))) > abs(np.cos(np.deg2rad(pdf.THETA_IMAGE)))]
        a['RADIUS'] = a.A_IMAGE*3.5*abs(np.sin(np.deg2rad(a.THETA_IMAGE)))
        b = pdf[abs(np.sin(np.deg2rad(pdf.THETA_IMAGE))) <= abs(np.cos(np.deg2rad(pdf.THETA_IMAGE)))]
        b['RADIUS'] = b.A_IMAGE*3.5*abs(np.cos(np.deg2rad(b.THETA_IMAGE)))
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
        if not np.sum(self.__pollutionDataFrame.values):
            self.__treated = True
            return
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
                            gs.Y_IMAGE.values-gs.RADIUS.values:gs.Y_IMAGE.values+gs.RADIUS.values+1,
                            gs.X_IMAGE.values-gs.RADIUS.values:gs.X_IMAGE.values+gs.RADIUS.values+1])
        if os.path.exists(self.__galaxyFile):
            os.system('rm '+self.__galaxyFile)
        ft.writeto(self.__galaxyFile, self.__galaxyData)
        # change the coordinates of x and y
        self.__galaxySeries.X_IMAGE -= (self.__pollutionTreatedData.shape[1]+1-self.__galaxyData.shape[1])/2+1
        self.__galaxySeries.Y_IMAGE -= (self.__pollutionTreatedData.shape[0]+1-self.__galaxyData.shape[0])/2+1
        self.__treated = True
        return

    def calculate_parameters(self):
        if not self.__treated:
            raise KeyError('you haven\'t treated the image')
        self.__galaxyStructuralParameter['A'] = np.sum(abs(self.__galaxyData-np.rot90(self.__galaxyData, 2)))/(2*np.sum(abs(self.__galaxyData)))
        self.__isInGalaxy = np.zeros(self.__galaxyData.shape)
        self.__galaxyMoment = np.zeros(self.__galaxyData.shape)
        _x = self.__galaxySeries.X_IMAGE.values
        _y = self.__galaxySeries.Y_IMAGE.values
        cxx = self.__galaxySeries.CXX_IMAGE.values
        cyy = self.__galaxySeries.CYY_IMAGE.values
        cxy = self.__galaxySeries.CXY_IMAGE.values
        r = 3.5
        for y in np.arange(self.__isInGalaxy.shape[0]):
            for x in np.arange(self.__isInGalaxy.shape[1]):
                self.__galaxyMoment[y][x] = self.__galaxyData[y][x]*((x-_x)**2+(y-_y)**2)
                if cxx*(x-_x)**2+cyy*(y-_y)**2+cxy*(x-_x)*(y-_y) <= r**2:
                    self.__isInGalaxy[y][x] = 1
        arg = np.argsort(self.__galaxyData, axis=None)[::-1]
        f = np.array([self.__galaxyData[arg[i]//self.__galaxyData.shape[0]][arg[i] % self.__galaxyData.shape[0]]
                      for i in np.arange(self.__galaxyData.size)
                      if self.__isInGalaxy[arg[i]//self.__galaxyData.shape[0]][arg[i] % self.__galaxyData.shape[0]]])
        m = np.array([self.__galaxyMoment[arg[i]//self.__galaxyData.shape[0]][arg[i] % self.__galaxyData.shape[0]]
                      for i in np.arange(self.__galaxyData.size)
                      if self.__isInGalaxy[arg[i]//self.__galaxyData.shape[0]][arg[i] % self.__galaxyData.shape[0]]])
        self.__galaxyStructuralParameter['G'] = sum([(2*l-f.size-1)*f[l]/(f.mean()*f.size*(f.size-1)) for l in np.arange(f.size)])
        ff = np.copy(f)
        for k in np.arange(1, f.size):
            ff[k] += ff[k-1]
            if ff[k] > f.sum()*0.2:
                self.__galaxyStructuralParameter['M'] = np.log10(np.sum(m[:k])/m.sum())
                break
        length = int(self.__galaxySeries.RADIUS.values)
        self.__galaxySurfaceBrightness = np.copy(self.__initialLoop[:length])
        # self.__galaxyMeanSurfaceBrightness = np.copy(self.__galaxySurfaceBrightness)
        print(length, self.__galaxySurfaceBrightness.size)
        found = False
        for i in np.arange(length-1):
            # self.__galaxyMeanSurfaceBrightness[i+1] += self.__galaxyMeanSurfaceBrightness[i]
            # self.__galaxyMeanSurfaceBrightness[i] /= (2*i+1)**2
            self.__galaxySurfaceBrightness[i+1] /= 8*(i+1)
            if self.__galaxySurfaceBrightness[i] < 2*self.__backgroundMean and not found:
                print(i)
                self.__galaxyStructuralParameter['C'] = np.sum(self.__initialLoop[:int(i*0.3)])/np.sum(self.__initialLoop[:i])
                found = True
        # self.__galaxyMeanSurfaceBrightness[length-1] /= (2*length-1)**2
        # eta = np.array([self.__galaxyMeanSurfaceBrightness[i]/self.__galaxySurfaceBrightness[i] for i in np.arange(length)])
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
