import numpy as np
import pyfits
from .helper_functions import chan_freq

import logging
logger = logging.getLogger(__name__)


def tophits_writer(spectra_slice, hit_indices, header, format='txt'):
    return None


def log_writer(filename, info_str):
    return None


def fits_writer(filename, header, fitsdata):
    return None

class GeneralWriter:
    """ """
    def __init__(self, filename='', mode='a'):
        with open(filename, mode) as myfile:
            self.filehandle = myfile
            self.filename = filename
        return None

    def close(self):
        if self.filehandle.closed:
            pass
        else:
            self.filehandle.close()

    def open(self, mode='a'):
        if self.filehandle.closed:
            with open(self.filename, mode) as myfile:
                self.filehandle = myfile
        elif self.filehandle.mode == mode:
            return
        else:
            self.close()
            with open(self.filename, mode) as myfile:
                self.filehandle = myfile

    def is_open(self):
        return not self.filehandle.closed

    def writable(self):
        if self.is_open() and (('w' in self.filehandle.mode) or ('a' in self.filehandle.mode)):
            return True
        else:
            return False

    def write(self, info_str, mode='a'):
        if (not 'w' in mode) and (not 'a' in mode):
            mode = 'a'
        if not self.writable():
            with open(self.filename, mode) as myfile:
                myfile.write(info_str)
                self.filehandle = myfile
        else:
            self.filehandle.write(info_str)

    def start_over(self):
        self.open('w')
        self.write('')
        self.open('a')



class FileWriter(GeneralWriter):
    """ """
    def __init__(self, filename):
        GeneralWriter.__init__(self, filename)
        self.tophit_count = 0

    def report_header(self, header):
        info_str = 'MJD: %18.12f\tRA: %18.12f\tDEC: %18.12f\tDELTAT: %18.12f\tDOPPLER: %18.12f\n'%(header['MJD'], header['RA'], header['DEC'], header['DELTAT'], header['DOPPLER'])
        self.write(info_str)
        self.write('--------------------------\n')


    def report_tophit(self, max_val, ind, ind_tuple, spec_slice, header):
        if not self.tophit_count:
            self.report_header(header)
        tdwidth =  len(max_val.maxsnr)
        info_str = 'Top Hit #%d: \t'%(self.tophit_count + 1)
        self.tophit_count += 1
        info_str += 'Drift Rate: %10.6f\t'%max_val.maxdrift[ind]
        self.write(info_str)
        info_str = 'Uncorrected Frequency: %14.6f\t'%chan_freq(header, ind, tdwidth, 0)
        info_str += 'Corrected Frequency: %14.6f\t'%chan_freq(header, ind, tdwidth, 1)
        info_str += 'Index: %d\n'%ind
        self.write(info_str)
        freq_start = chan_freq(header, ind_tuple[0], tdwidth, 0)
        freq_end = chan_freq(header, ind_tuple[1]-1, tdwidth, 0)
        info_str = 'Freqs: (%14.6f ... %14.6f), DELTAF: %12.8f (MHz)\n'%(freq_start, freq_end, header['DELTAF'])
        self.write(info_str)
        for i in range(0, spec_slice.shape[0]):
            info_str = ''
            for j in range(0, spec_slice.shape[-1]):
                info_str += '%14.6f '%spec_slice[i, j]
            info_str += '\n'
            self.write(info_str)
        self.write('\n')
        return self


class LogWriter(GeneralWriter):
    """ """
    def report_candidate(self, info_str):
        return None

    def info(self, info_str):
        self.write(info_str + '\n')
        return None



