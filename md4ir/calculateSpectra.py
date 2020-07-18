"""
Calculate vibrational spectra from any 1 or 3-column variable
This code is part of md4ir, available at:
https://gitlab.ethz.ch/paenurke/md4ir
"""

__all__ = [ 'calcSpec' ]

import numpy as np
import scipy.signal as sps
from . import fileHandling


def calcSpec(name, file, scale, timestep, time_start, time_end, wn_start, wn_end, method):
    """"Calculates the power or IR spectrum (depends on the input type)
    """
    # Set up variables
    files = [i.split(',') for i in file]
    if name == 'dummy':
        name = [str(i[0][:-4]) for i in files] # assuming .xxx extension
    else:
        name = [name] # define as list for easier handling later

    # Run routines
    
    for n in np.arange(0,len(files)):

        # Calculate the spectrum
        if len(files[n]) > 1: # average spectra if more than one files given per -f flag
            signal = 0
            for f in files[n]:
                data = fileHandling.read_cols(f)
                spec = _Calc(data, scale, timestep, time_start, time_end, wn_start, wn_end)
                spec.gen_spectrum(method, norm = False)
                signal += spec.signal
            spec.signal = signal / np.amax(signal)
        else:
            data = fileHandling.read_cols(files[n][0])
            spec = _Calc(data, scale, timestep, time_start, time_end, wn_start, wn_end)
            spec.gen_spectrum(method)

        # Write the spectrum into a file
        with open(name[n]+'_spectrum.txt','w',encoding='utf8') as f:
            for i in np.arange(0,len(spec.signal)):
                    f.write("%f %f\n" % (spec.wavenumbers[i],  spec.signal[i]))
            f.close()


# Classes

class _Calc:
    
    def __init__(self, data, scale, timestep, t_start, t_end, wn_start, wn_end):
        """
        data: list of lists of the variable (dipoles, coordinates, etc.)
        scale: scaling of the wavenumbers
        timestep: timestep of the simulation in fs
        t_start: start time of the evaluation
        t_end: end time of the evaluation in ps
        wn_start: start-point of the spectrum in cm-1
        wn_end: end-point of the spectrum in cm-1
        """
        self.timestep = timestep * 1e-15 #conversion to s
        self.t_start = t_start
        self.t_end = t_end
        self.data = [val for idx,val in enumerate(data) if idx >= t_start/timestep*1000 and idx <= t_end/timestep*1000]
        self.wn_start = wn_start
        self.wn_end = wn_end
        self.wavenumbers = None
        self.signal = None
        self.Hz_2invcm = 3.335641e-11
        self.scale = scale
        
    def dx_dt(self):
        """
        Returns the time derivative of the values in the array
        """
        length = len(self.data[0])
        dx_dt = np.zeros((len(self.data)-1,length))
        for i in np.arange(0,len(self.data)-1):
            for j in np.arange(0,length):
                dx_dt[i][j] = (float(self.data[i+1][j])-float(self.data[i][j]))/self.timestep
        return dx_dt

    def autocorr_numpy(self):
        """
        Calculates the autocorrelation function with standard numpy tools
        """
        dx_dt = self.dx_dt().transpose()
        signal = np.zeros(len(dx_dt[0])-1, dtype=complex)
        for d in range(0,len(dx_dt)):
            acf = np.correlate(dx_dt[d],dx_dt[d],mode='full')
            halfLength = int((len(acf)+1)/2)
            acf = acf[halfLength:]
            bh = sps.blackmanharris(len(acf))
            acf_bh = (acf.transpose()*bh).transpose()
            signal_temp = np.fft.fft(acf_bh, axis=0)
            signal += signal_temp
        return signal, acf_bh

    def gen_spectrum(self, method = 'numpy', norm = True):
        """
        Calculates the spectrum in cm-1 frame
        """
        # Calculate the FFT
        if method == 'numpy': # other methods can be added by elif statements
            signal, acf_bh = self.autocorr_numpy()
        # Calculate the wavenumbers and treat the data
        wavenumbers = np.fft.fftfreq(len(acf_bh), d=self.timestep)
        wavenumbers = wavenumbers * self.Hz_2invcm * self.scale # convert from Hz to cm-1 and scale
        signal = signal[wavenumbers > self.wn_start] # take signals corresponding to wavenumbers above wn_end
        wavenumbers = wavenumbers[wavenumbers > self.wn_start] # take wavenumbers above wn_end
        signal = signal[wavenumbers < self.wn_end] # take signals corresponding to wavenumbers below wn_end
        self.wavenumbers = wavenumbers[wavenumbers < self.wn_end] # take wavenumbers below wn_end
        self.signal = abs(signal) # take the absolute
        if norm == True:
            self.signal = self.signal / np.amax(self.signal) # normalize to 1
