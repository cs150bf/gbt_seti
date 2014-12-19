import pyfits
import os
import numpy as np
import math
import copy
from pkg_resources import resource_filename
from . import barycenter

import logging
logger = logging.getLogger(__name__)

class FITSHandle:
    """ """
    def __init__(self, filename=None, size_limit = 1024.0, split_dir='.'):
        if filename and os.path.isfile(filename):
            self.filename = filename
            self.filestat = os.stat(filename)
            filesize = self.filestat.st_size/(1024.0**2)
            if filesize > size_limit:
                logger.info("The file is of size %f MB, exceeding our size limit %f MB. Split needed..."%(filesize, size_limit))
                self.split = True
                self.fits_list = self.split_fits(filename, split_dir=split_dir, size_limit=size_limit)
            else:
                logger.debug("File size %f MB within range %f MB, okay..."%(filesize, size_limit))
                self.split = False
                fits_obj = FITS(filename)
                self.fits_list = [fits_obj]
            self.split_filenames = str([fits.filename if fits.status else '' for fits in self.fits_list])
            if self.split_filenames == '[]':
                logger.error("Failed to initiate a FITS instance, aborting...")
                self.status = False
                return None
            else:
                self.status = True
        else:
            logger.error("FITS file %s doesn\'t exists, please check!"%filename)
            self.status = False
            return None

    def get_info(self):
        return ""

    @staticmethod
    def split_fits(filename=None, split_dir='.', size_limit = 1024.0):
        fits_list = []
        filestat = os.stat(filename)
        filesize = filestat.st_size/(1024.0**2)
        if filesize <= size_limit:
            logger.error("This file is only %f MB, smaller than our size limit %f MB, no split."%(filesize, size_limit))
            return []
        try:
            bighdulist = pyfits.open(filename, memmap=True)
            first_row = bighdulist[0]
            header = first_row.header
        except:
            logger.error("Error encountered when trying to open FITS file %s"%filename)
            return []

        fn = filename[filename.rfind('/')+1:filename.rfind('.fits')]
        deltaf = header['DELTAF']
        fftlen = header['NAXIS1']
        fcntr = header['FCNTR']
        frange = [fcntr - fftlen*deltaf/2, fcntr + fftlen*deltaf/2]

        nfiles_min = int(math.ceil(filesize/size_limit))
        new_width_max = fftlen/nfiles_min
        new_width = 2**math.floor(np.log2(new_width_max))
        nfiles = int(math.ceil(fftlen/new_width))
        new_files = []
        new_fcntrs = []
        new_filenames = []
        indices = []
        new_primary_header = copy.deepcopy(header)
        to_create = []
        for i in range(0, nfiles):
            new_filename = split_dir + '/' + fn + '_%d'%i + '.fits'
            new_filenames.append(new_filename)
            new_fcntr_tmp = frange[0] + deltaf * new_width * (i + 0.5)
            new_fcntrs.append(new_fcntr_tmp)
            new_primary_header['FCNTR'] = new_fcntr_tmp
            ind = (i*new_width, min(fftlen, (i+1)*new_width))
            indices.append(ind)
            if os.path.isfile(new_filename):
                logger.error("file %s already existed!"%new_filename)
                to_create.append(False)
                continue
            to_create.append(True)
            data = first_row.data[0][ind[0]:ind[1]]
            prihdu = pyfits.PrimaryHDU([data], header = new_primary_header)
            prihdu.writeto(new_filename)
            logger.info("Created new file: %s"%new_filename)

        for i, ohdu in enumerate(bighdulist[1:]):
            logger.debug("Dealing with row %d"%i)
            new_header = copy.deepcopy(ohdu.header)
            for j, new_filename in enumerate(new_filenames):
                if not to_create[j]:
                    continue
                new_header['FCNTR'] = new_fcntrs[j]
                ind = indices[j]
                data = ohdu.data[0][ind[0]:ind[1]]
                pyfits.append(new_filename, [data], new_header)

        for new_filename in new_filenames:
            fits_obj = FITS(new_filename)
            fits_list.append(fits_obj)

        return fits_list



class FITS:
    """ """
    def __init__(self, filename=None, size_limit = 1024.0):
        if filename and os.path.isfile(filename):
            self.filename = filename
            self.filestat = os.stat(filename)
            filesize = self.filestat.st_size/(1024.0**2)
            if filesize > size_limit:
                logger.info("The file is of size %f MB, exceeding our size limit %f MB. Aborting..."%(filesize, size_limit))
                return None
            try:
                hdulist = pyfits.open(filename, memmap=True)
                #self.info = hdulist.info()
                first_row = hdulist[0]
                header = first_row.header
                hdulist.close()
            except:
                logger.error("Error encountered when trying to open FITS file %s"%filename)
                self.status = False
                return None
            self.fftlen = header['NAXIS1']
            self.tsteps_valid = len(hdulist)
            self.tsteps = int(math.pow(2, math.ceil(np.log2(math.floor(self.tsteps_valid)))))
            self.obs_length = self.tsteps_valid * header['DELTAT']
            self.tdwidth = self.fftlen + 8*self.tsteps
            self.drift_rate_resolution = (1e6 * header['DELTAF']) / self.obs_length
            self.nom_max_drift = self.drift_rate_resolution * self.tsteps_valid
            self.header = barycenter.correct(header, self.obs_length)
            logger.info('barycenter done for fits file %s! baryv: %f'%(filename, self.header['baryv']))
            # some default values
            self.original_vals= {'tsteps_valid': self.tsteps_valid, 'tsteps': self.tsteps,
                                 'tdwidth': self.tdwidth, 'fftlen':self.fftlen}
            self.compressed_t = False
            self.compressed_f = False
            self.status = True
            return None

    def load_data(self, max_search_rate=None, bw_compress_width=None, logwriter=None):
        hdulist = pyfits.open(self.filename, memmap=True)
        spectra = np.empty((self.tsteps_valid, self.fftlen), dtype=np.float32)

        if max_search_rate and self.nom_max_drift > max_search_rate:
            logger.info("nominal max drift rate greater than allowed.... decimating.")
            self.compressed_t = True
            decimate_factor = math.floor(self.nom_max_drift / max_search_rate)
            self.tsteps_valid = int((self.tsteps_valid - decimate_factor)/decimate_factor)
            self.tsteps = int(math.pow(2, math.ceil(np.log2(math.floor(self.tsteps_valid)))))
            self.tdwidth = self.fftlen + 8*self.tsteps
            self.nom_max_drift = self.drift_rate_resolution * self.tsteps_valid
            for i in range(0, self.tsteps_valid+1):
                print "loading row %d..."%(i*decimate_factor)
                np.copyto(spectra[i], hdulist[i*decimate_factor].data[0])
                for k in range(1, decimate_factor):
                    # ???
                    logger.debug("Decimation: adding row %d"%(i*decimate_factor + k))
                    spectra[i, :] += hdulist[i*decimate_factor + k].data[0]
                spectra = spectra[0:self.tsteps_valid, :]
        else:
            for i in range(0, self.tsteps_valid):
                print "loading row %d"%i
                np.copyto(spectra[i], hdulist[i].data[0])

        compressed_spectra = None
        if bw_compress_width:
            compressed_spectra = matrix_compression(spectra, compression_width = bw_compress_width, axis=1, method='max')
        if compressed_spectra is None:
            compressed_spectra = spectra
        else:
            self.compressed_f = True
            self.fftlen = compressed_spectra.shape[-1]
            self.tdwidth = self.fftlen + 8*self.tsteps

        drift_indexes = self.load_drift_indexes()
        return compressed_spectra, drift_indexes


    def load_drift_indexes(self):
        n = int(np.log2(self.tsteps))
        di_array = np.genfromtxt(resource_filename('dedoppler', 'drift_indexes/drift_indexes_array_%d.txt'%n), delimiter='\t')
        ts2 = self.tsteps/2
        drift_indexes = di_array[self.tsteps_valid - 1 - ts2, 0:self.tsteps_valid]
        return drift_indexes

    def get_info(self):
        return ""



def matrix_compression(matrix_original, compression_width=1, axis=0, method='add'):
    """
    matrix_original: assumed to be a well-formed 2-d array
    compression_width: an integer
    axis: 0 or 1, if 0 sum along column, if 1 then sum along row
    """
    matrix_shape = matrix_original.shape
    if len(matrix_shape) != 2:
        logger.error("Sorry, this function can only handle 2-dim arrays, aborting...")
        return None
    compression_width = int(max(1, math.floor(compression_width)))
    if compression_width > 1 and compression_width < matrix_shape[axis]:
        target_ind = range(0,  int(matrix_shape[axis]), compression_width)
        if method=='max':
            matrix_compressed = np.empty_like(matrix_original[:, target_ind] if axis else matrix_original[target_ind, :])
            for i, val in enumerate(target_ind):
                temp = np.amax(matrix_original[:, val:val+compression_width] if axis else matrix_original[val:val+compression_width], axis=axis)
                np.copyto( matrix_compressed[:, i] if axis else matrix_compressed[i, :], temp)
        elif method=='add':
            matrix_compressed = matrix_original[:, target_ind] if axis else  matrix_original[target_ind, :]
            for i in target_ind:
                ind1 = int(i/compression_width)
                if axis: # sum along rows
                    np.copyto(matrix_compressed[:, ind1], matrix_original[:, i])
                else:
                    np.copyto(matrix_compressed[ind1, :], matrix_original[i, :])
                for j in range(1, compression_width):
                    ind2 = i + j
                    if (not axis) and ind2 < matrix_shape[axis]:
                        matrix_compressed[ind1, :] += matrix_original[ind2, :]
                    elif ind2 < matrix_shape[axis]:
                        matrix_compressed[:, ind1] += matrix_original[:, ind2]
        else:
            logger.error("Method unkonwn/unimplmeneted, returning None...")
            return None
        logger.info("Compression done, width: %d, axis: %d"%(compression_width, axis))
    else:
        logger.info("No compression to be performed, returning None...")
        matrix_compressed = None
    return matrix_compressed
