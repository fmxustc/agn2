import astropy.io.fits as ft
import numpy as np
import os
import matplotlib.pyplot as plt


class Galaxy(object):

    def __init__(self, file):
        # initial data and attributes
        self._hdu = ft.open(file)
        self._header = self._hdu[0].header
        self._initialData = np.copy(self._hdu[0].data)
        self._hdu.close()
        self._initialCenter = self._initialData.shape[0]/2
        self._boxSize = 250
        self._initialLoop = np.array([np.sum([
            self._initialData[self._initialCenter-i, self._initialCenter-i:self._initialCenter+i],
            self._initialData[self._initialCenter+i, self._initialCenter-i:self._initialCenter+i],
            self._initialData[self._initialCenter-i:self._initialCenter+i, self._initialCenter-i],
            self._initialData[self._initialCenter-i:self._initialCenter+i, self._initialCenter+i],
        ]) for i in np.arange(self._boxSize)])
        self._initialLoop[0] = self._initialData[self._initialCenter][self._initialCenter]
        # truncate data and attributes
        self._loopRatio = []
        self._truncateData = []
        if os.path.exists('truncate.fits'):
            os.system('rm truncate.fits')
        self._truncateFile = 'truncate.fits'
        # processed data and galaxy real attributes
        self._galaxyData = None
        self._galaxySize = 0
        self._structuralParameters = {
            'G': None,
            'M': None,
            'A': None,
            'C': None
        }

    def truncate(self):
        self._loopRatio = [self._initialLoop[k]/sum(self._initialLoop[:k]) for k in np.arange(1,self._initialLoop.size)]
        for lr in self._loopRatio:
            if lr < 1e-03:
                self._galaxySize = self._loopRatio.index(lr)
                break
        if not self._galaxySize:
            raise SystemError
        self._truncateData = self._initialData[
                    self._initialCenter-self._galaxySize:self._initialCenter+self._galaxySize,
                    self._initialCenter-self._galaxySize:self._initialCenter+self._galaxySize]
        ft.writeto(self._truncateFile, self._truncateData)

        return

    def show_loop_ratio(self):
        if os.path.exists('truncate.fits'):
            plt.title('loopRatio')
            plt.plot([1e-03]*len(self._loopRatio))
            plt.plot(self._loopRatio)
            plt.show()
        else:
            print('You may not run the function truncate()')
        return

    def show_truncate_image(self):
        if os.path.exists('truncate.fits'):
            os.system('ds9 '+self._truncateFile)
        else:
            print('You may not run the function truncate()')
        return
