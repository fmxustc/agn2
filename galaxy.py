import astropy.io.fits as ft
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import subprocess


class Galaxy(object):

    def __init__(self, file):
        # initial data and attributes
        self.__initialFile = file
        self.__hdu = ft.open(file)
        self.__header = self.__hdu[0].header
        self.__initialData = np.copy(self.__hdu[0].data)
        self.__hdu.close()
        self.__initialCentroid = self.__initialData.shape[0]/2
        self.__boxRadius = min(200, self.__initialCentroid-1)
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
            if lr < 1e-03:
                self.__truncateRadius = self.__loopRatio.index(lr)
                break
        if not self.__truncateRadius:
            self.__truncateRadius = np.argmin(self.__loopRatio[:200])
            # raise AttributeError('can\'t get the radius from truncate data')
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
        subprocess.call('sex '+self.__truncateFile, shell=True, executable='/bin/bash')
        self.__SExOutput = subprocess.Popen('sex '+self.__truncateFile,
                                            shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if not os.path.exists(self.__truncateFile):
            raise FileExistsError('your SExtractor software may thread some errors')
        tr = open(self.__pollutionFile, 'r')
        tw = open('pollution.csv', 'w')
        tw.write('MAG_AUTO,X_IMAGE,Y_IMAGE,A_IMAGE,B_IMAGE,THETA_IMAGE,RADIUS,DIST\n')
        ct = self.__truncateRadius
        for line in tr.readlines():
            elements = line.split()
            for elem in elements:
                tw.write(elem+',')
            cos = abs(np.cos(np.deg2rad(float(elements[5]))))
            sin = abs(np.sin(np.deg2rad(float(elements[5]))))
            if cos > sin:
                tw.write(str(float(elements[3])*cos*6))
            else:
                tw.write(str(float(elements[3])*sin*6))
            tw.write(','+str(abs(float(elements[1])-ct)+abs(float(elements[2])-ct)))
            tw.write('\n')
        tr.close()
        tw.close()
        self.__pollutionFile = 'pollution.csv'
        bg = str(self.__SExOutput.stdout.readlines()[14]).split()
        self.__backgroundMean = float(bg[2])
        self.__backgroundRMS = float(bg[4])
        self.__backgroundThreshold = float(bg[7])
        self.__found = True
        return

    def eliminate_pollution(self):
        if not self.__found:
            raise KeyError('you may not use SExtractor find the pollutions')
        pdf = pd.read_csv(self.__pollutionFile)
        pdf = pdf.sort_index(by='DIST')
        pdf.index = range(len(pdf.index))
        gs = pdf.ix[0]
        pdf = pdf.drop(pdf.index[0])
        gp = pdf[(pdf.MAG_AUTO < gs.MAG_AUTO*0.2) &
                 (pdf.A_IMAGE/pdf.B_IMAGE < 2.5)]
        pollutions = gp.iterrows()
        self.__pollutionTreatedData = np.copy(self.__truncateData)
        edge = self.__pollutionTreatedData.shape[0]
        for pollution in pollutions:
            p = pollution[1]
            x_min = max(p.X_IMAGE-p.RADIUS, 0)
            x_max = min(p.X_IMAGE+p.RADIUS, edge-1)
            y_min = max(p.Y_IMAGE-p.RADIUS, 0)
            y_max = min(p.Y_IMAGE+p.RADIUS, edge-1)
            bg = np.random.normal(self.__backgroundMean, self.__backgroundRMS, int(np.pi*p.A_IMAGE*p.B_IMAGE*20))
            _x = p.X_IMAGE
            _y = p.Y_IMAGE
            a = p.A_IMAGE
            b = p.B_IMAGE
            th = np.deg2rad(p.THETA_IMAGE)
            cxx = np.cos(th)**2/a**2+np.sin(th)**2/b**2
            cyy = np.sin(th)**2/a**2+np.cos(th)**2/b**2
            cxy = 2*np.sin(th)*np.cos(th)*(1/a**2-1/b**2)
            r = 4
            cnt = 0
            for i in np.arange(y_min, y_max):
                for j in np.arange(x_min, x_max):
                    if cxx*(j-_x)**2+cyy*(i-_y)**2+cxy*(i-_y)*(j-_x) <= r**2:
                        self.__pollutionTreatedData[i][j] = bg[cnt]
                        cnt += 1
        if os.path.exists(self.__pollutionTreatedFile):
            os.system('rm '+self.__pollutionTreatedFile)
        ft.writeto(self.__pollutionTreatedFile, self.__pollutionTreatedData)
        self.__galaxyData = np.copy(self.__pollutionTreatedData[
                            max(gs.Y_IMAGE-gs.RADIUS, 0):min(gs.Y_IMAGE+gs.RADIUS+1, edge),
                            max(gs.X_IMAGE-gs.RADIUS, 0):min(gs.X_IMAGE+gs.RADIUS+1, edge)])
        if self.__galaxyData.shape[0]-self.__galaxyData.shape[1]:
            m = np.min(self.__galaxyData.shape)
            self.__galaxyData = self.__galaxyData[self.__galaxyData.shape[0]-m:self.__galaxyData.shape[0],
                                                  self.__galaxyData.shape[1]-m:self.__galaxyData.shape[1]]
        if os.path.exists(self.__galaxyFile):
            os.system('rm '+self.__galaxyFile)
        ft.writeto(self.__galaxyFile, self.__galaxyData)
        gs.X_IMAGE = np.argmax(self.__galaxyData) % self.__galaxyData.shape[0]
        gs.Y_IMAGE = np.argmax(self.__galaxyData) // self.__galaxyData.shape[0]
        self.__pollutionDataFrame = pdf
        self.__galaxySeries = gs
        self.__galaxyPollutions = gp
        self.__treated = True
        return

    def calculate_parameters(self):
        if not self.__treated:
            raise KeyError('you haven\'t treated the image')
        self.__galaxyStructuralParameter['A'] = np.sum(abs(self.__galaxyData-np.rot90(self.__galaxyData, 2)))/(2*np.sum(abs(self.__galaxyData)))
        self.__isInGalaxy = np.zeros(self.__galaxyData.shape)
        self.__galaxyMoment = np.zeros(self.__galaxyData.shape)
        _x = self.__galaxySeries.X_IMAGE
        _y = self.__galaxySeries.Y_IMAGE
        a = self.__galaxySeries.A_IMAGE
        b = self.__galaxySeries.B_IMAGE
        th = np.deg2rad(self.__galaxySeries.THETA_IMAGE)
        cxx = np.cos(th)**2/a**2+np.sin(th)**2/b**2
        cyy = np.sin(th)**2/a**2+np.cos(th)**2/b**2
        cxy = 2*np.sin(th)*np.cos(th)*(1/a**2-1/b**2)
        r = 4
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
        self.__galaxySurfaceBrightness = np.zeros(self.__galaxyData.shape[0])
        for i in np.arange(self.__galaxyData.shape[0]):
            for j in np.arange(self.__galaxyData.shape[1]):
                dist = max(abs(i-self.__galaxySeries.Y_IMAGE), abs(j-self.__galaxySeries.X_IMAGE))
                self.__galaxySurfaceBrightness[dist] += self.__galaxyData[i][j]
        found = False
        for i in np.arange(min(self.__galaxySeries.RADIUS, self.__galaxyData.shape[0])-1):
            self.__galaxySurfaceBrightness[i + 1] /= 8 * (i + 1)
            if self.__galaxySurfaceBrightness[i] < 5 * self.__backgroundRMS and not found:
                self.__galaxyStructuralParameter['C'] = np.sum(self.__initialLoop[:int(i * 0.3)]) / np.sum(
                    self.__initialLoop[:i])
                found = True
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
