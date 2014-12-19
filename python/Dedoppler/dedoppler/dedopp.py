import numpy as np
import pyfits #fitsio
import sys
import os
from numba import double, jit, autojit
from . import fits_wrapper
from . import file_writers
from .helper_functions import chan_freq

import logging
logger = logging.getLogger(__name__)


class max_vals:
    def __init__(self):
        self.maxsnr = None
        self.maxdrift = None
        self.maxsmooth = None
        self.maxid = None


class DedopplerTask:
    """ """
    def __init__(self, fitsfile, min_drift, max_drift, snr, bw, rfithresh=25.0, rfiwindow = 2, split_dir='/tmp'):
        self.min_drift = min_drift
        self.max_drift = max_drift
        self.snr = snr
        self.bw = bw
        self.rfithresh = rfithresh
        self.rfiwindow = rfiwindow
        self.fits_handle = fits_wrapper.FITSHandle(fitsfile, split_dir=split_dir)
        if (self.fits_handle is None) or (self.fits_handle.status is False):
            logger.error("FITS file error, aborting...")
            self.status = False
            return None
        logger.info(self.fits_handle.get_info())
        logger.info("A new Dedoppler Task instance created!")
        self.status = True
        fn = fitsfile[fitsfile.rfind('/')+1:fitsfile.rfind('.fits')]
        self.logwriter = file_writers.LogWriter('%s/%s.log'%(split_dir, fn))
        self.filewriter = file_writers.FileWriter('%s/%s.dat'%(split_dir, fn))
        return None

    def get_info(self):
        info_str = "FITS file: %s\nSplit FITS files: %s\ndrift rates (min, max): (%f, %f)\nSNR: %f\nbw: %f\n"%(self.fits_handle.filename,
            self.fits_handle.split_filenames, self.min_drift, self.max_drift,
            self.snr, self.bw)
        return info_str

    def search(self):
        logger.info("Start searching...")
        logger.debug(self.get_info())
        for target_fits in self.fits_handle.fits_list:
            self.search_fits(target_fits)

    def search_fits(self, fits_obj):
        logger.info("Start searching for %s"%fits_obj.get_info())
        self.logwriter.info("Start searching for %s"%fits_obj.get_info())
        spectra, drift_indexes = fits_obj.load_data(bw_compress_width = self.bw, logwriter=self.logwriter)
        tsteps = fits_obj.tsteps
        tsteps_valid = fits_obj.tsteps_valid
        tdwidth = fits_obj.tdwidth
        fftlen = fits_obj.fftlen
        nframes = tsteps_valid

        spectrum = np.empty_like(spectra[0])
        spectrum_sum = np.empty_like(spectrum)

        # allocate array for dedopplering 
        # init dedopplering array to zero 
        tree_dedoppler = np.zeros(tsteps * tdwidth,
                                  dtype=np.float32)

        # allocate array for holding original 
        tree_dedoppler_original = np.empty_like(tree_dedoppler)

        #/* build index mask for in-place tree doppler correction */
        ibrev = np.zeros(tsteps, dtype=np.int32)

        for i in range(0, tsteps):
            ibrev[i] = bitrev(i, int(np.log2(tsteps)))

        for i in range(0, nframes):
            sind = tdwidth*i + tsteps*4
            cplen = fftlen
            np.copyto(tree_dedoppler[sind: sind + cplen], spectra[i])
            #//load end of current spectra into left hand side of next spectra 
            sind = i * tdwidth
            np.copyto(tree_dedoppler[sind: sind + tsteps*4],
                  spectra[i, fftlen-(tsteps*4):fftlen])


        #/* allocate array for negative doppler rates */
        tree_dedoppler_flip = np.empty_like(tree_dedoppler)

        #/* malloc stat mask if we need to */
        #/* zero out stat mask */
        stat_mask = np.zeros(tdwidth, dtype='uint8')

        for i in range(0, tsteps*4):
            stat_mask[i] = 1
        for i in range(tdwidth - tsteps*4, tdwidth):
            stat_mask[i] = 1


        max_val = max_vals()
        if max_val.maxsnr == None:
            max_val.maxsnr = np.zeros(tdwidth, dtype=np.float32)
        if max_val.maxdrift == None:
            max_val.maxdrift = np.zeros(tdwidth, dtype=np.float32)
        if max_val.maxsmooth == None:
            max_val.maxsmooth = np.zeros(tdwidth, dtype='uint8')
        if max_val.maxid == None:
            max_val.maxid = np.zeros(tdwidth, dtype='uint64')


        #/* populate original array */
        np.copyto(tree_dedoppler_original, tree_dedoppler)

        print "Running interference rejection\n"
        logger.debug('tree_dedoppler shape - %s'%str(tree_dedoppler.shape))
        rfi_mask = rfirej(tree_dedoppler, tdwidth, nframes, tsteps, self.rfiwindow, self.rfithresh)

        #/* populate neg doppler array */
        np.copyto(tree_dedoppler_flip, tree_dedoppler)

        logger.info("Doppler correcting forward...")
        taylor_flt(tree_dedoppler, tsteps * tdwidth, tsteps)
        logger.info("done...")
        if (tree_dedoppler == tree_dedoppler_original).all():
             logger.error("taylor_flt has no effect?")
        else:
             logger.debug("tree_dedoppler changed")

        # candidate search!
        max_val.maxsnr = np.zeros(tdwidth, dtype = np.float32)
        max_val.maxdrift = np.zeros(tdwidth, dtype = np.float32)


        for k in range(0, tsteps_valid):
            indx  = (ibrev[drift_indexes[k]] * tdwidth)
            #/* SEARCH POSITIVE DRIFT RATES */
            spectrum = tree_dedoppler[indx: indx+tdwidth]
            mean_val, stddev = comp_stats(spectrum, tdwidth, stat_mask)

            #/* normalize */
            for i in range(0, tdwidth):
                spectrum[i] = (spectrum[i] - mean_val)/stddev

            m = 0 # hack
            if m==0:
                specstart = (tsteps*4)
                specend = tdwidth - (tsteps * 6)
            else:
                specstart = (tsteps*2)
                specend = tdwidth - (tsteps * 6)

            channel = m # hack
            drift_rate = k*fits_obj.drift_rate_resolution
            logger.info('drift_rate: %f'%drift_rate)
            if drift_rate < self.max_drift and drift_rate > self.min_drift:
                n_candi, max_val = candsearch(spectrum, specstart, specend, self.snr, \
                                      drift_rate, fits_obj.header, fftlen, tdwidth, channel, max_val, 0)
                info_str = "found %d candidates at drift rate %15.15f\n"%(n_candi, drift_rate)
                logger.info(info_str)
                self.logwriter.info(info_str)
                self.filewriter = tophitsearch(tree_dedoppler_original, max_val, tsteps, nframes, tdwidth, fits_obj.header, split_dir = '/queencow/dedoppler_test', rfithresh=self.rfithresh, logwriter=self.logwriter, filewriter=self.filewriter)

        #/* copy back original array */
        np.copyto(tree_dedoppler, tree_dedoppler_flip)

        #/* Flip matrix across X dimension to search negative doppler drift rates */
        FlipX(tree_dedoppler_flip, tdwidth, tsteps)

        print "Doppler correcting reverse...\n"
        taylor_flt(tree_dedoppler_flip, tsteps * tdwidth, tsteps)
        print "done...\n"

        for k in range(0, tsteps_valid):
            indx  = ibrev[drift_indexes[k]] * tdwidth

            #/* SEARCH NEGATIVE DRIFT RATES */
            spectrum = tree_dedoppler_flip[indx: indx + tdwidth]
            mean_val, stddev = comp_stats(spectrum, tdwidth, stat_mask) 


            #/* normalize */
            for i in range(0, tdwidth):
                spectrum[i] = (spectrum[i] - mean_val)/stddev

            if m==0:
                specstart = (tsteps*4)
                specend = tdwidth - (tsteps * 6)
            else:
                specstart = tsteps * 4
                specend = tdwidth - (tsteps * 4)

            drift_rate = -1 * k * fits_obj.drift_rate_resolution
            logger.debug("Drift rate: %f"%drift_rate)
            if drift_rate < self.max_drift and drift_rate > self.min_drift:
                 n_candi, max_val = candsearch(spectrum, specstart, specend, self.snr, \
                                drift_rate, fits_obj.header, \
                                fftlen, tdwidth, channel, max_val, 1)
                 info_str = "found %d candidates at drift rate %15.15f\n"%(n_candi, drift_rate)
                 logger.info(info_str)
                 self.logwriter.info(info_str)
                 self.filewriter = tophitsearch(tree_dedoppler_original, max_val, tsteps, nframes, tdwidth, fits_obj.header, split_dir = '/queencow/dedoppler_test', rfithresh=self.rfithresh, logwriter=self.logwriter, filewriter=self.filewriter)



#  ======================================================================  #
#  This function bit-reverses the given value "inval" with the number of   #
#  bits, "nbits".    ----  R. Ramachandran, 10-Nov-97, nfra.               #
#  python version ----  H. Chen   Modified 2014                            #
#  ======================================================================  #
@jit
def bitrev(inval, nbits):
    if nbits <= 1:
        ibitr = inval
    else:
        ifact = 1
        for i in range(1, nbits):
           ifact *= 2
        k = inval
        ibitr = (1 & k) * ifact
        for i in range(2, nbits+1):
            k /= 2
            ifact /= 2
            ibitr += (1 & k) * ifact
    return ibitr

#  ======================================================================  #
#  This function bit-reverses the given value "inval" with the number of   #
#  bits, "nbits".                                                          #
#  python version ----  H. Chen   Modified 2014                            #
#  reference: stackoverflow.com/questions/12681945                         #
#  ======================================================================  #
def bitrev2(inval, nbits, width=32):
    b = '{:0{width}b}'.format(inval, width=width)
    ibitr = int(b[-1:(width-1-nbits):-1], 2)
    return ibitr


@jit
def AxisSwap(inbuf, outbuf, nchans, NTSampInRead):
    #long int    j1, j2, indx, jndx;
    for j1 in range(0, NTSampInRead):
        indx  = j1 * nchans
        for j2 in range(nchans-1, -1, -1):
            jndx = j2 * NTSampInRead + j1
            outbuf[jndx]  = inbuf[indx+j2]


@jit
def FlipBand(outbuf, nchans, NTSampInRead):
    temp = np.zeros(nchans*NTSampInRead, dtype=np.float32)

    indx  = (nchans - 1);
    for i in range(0, nchans):
        jndx = (indx - i) * NTSampInRead
        kndx = i * NTSampInRead
        np.copyto(temp[jndx: jndx+NTSampInRead], outbuf[kndx + NTSampInRead])
    #memcpy(outbuf, temp, (sizeof(float)*NTSampInRead * nchans));
    outbuf = temp
    return


@jit
def FlipX(outbuf, xdim, ydim):
    temp = np.empty_like(outbuf[0:xdim])
    logger.debug("temp array dimension: %s"%str(temp.shape))

    for j in range(0, ydim):
        revi = xdim - 1
        indx = j * xdim
        np.copyto(temp, outbuf[indx:indx+xdim])
        np.copyto(outbuf[indx: indx+xdim], temp[::-1])
    return


#  ======================================================================  #
#  This is a function to Taylor-tree-sum a data stream. It assumes that    #
#  the arrangement of data stream is, all points in first spectra, all     #
#  points in second spectra, etc...  Data are summed across time           #
#                     Original version: R. Ramachandran, 07-Nov-97, nfra.  #
#                     Modified 2011 A. Siemion float/64 bit addressing     #
#                     Modified 2014 H. Chen python version                 #
#  outbuf[]       : input array (float), replaced by dedispersed data      #
#                   at the output                                          #
#  mlen           : dimension of outbuf[] (long int)                       #
#  nchn           : number of frequency channels (long int)                #
#                                                                          #
#  ======================================================================  #
@autojit
def taylor_flt(outbuf, mlen, nchn):
    '''
    Parameters:
        outbuf       : input array (float), replaced by dedispersed data 
                       at the output
        mlen         : dimension of outbuf[] (long int) 
        nchn         : number of frequency channels (long int)
    '''
    nsamp = (mlen/nchn) - (2*nchn)
    npts = nsamp + nchn
    ndat1 = nsamp + 2.0 * nchn
    nstages = int(np.log2(nchn))
    nmem = 1.0
    for istages in range(0, nstages):
        nmem  *= 2.0
        nsec1  = nchn/nmem
        nmem2  = nmem - 2.0
        for isec in range(0, int(nsec1)):
            ndelay = -1.0
            koff = isec * nmem
            for ipair in range(0, int(nmem2+1), 2):
                ioff1 = (bitrev(ipair, istages+1) + koff) * ndat1
                i2 = (bitrev(ipair+1, istages+1) + koff) * ndat1
                ndelay += 1
                ndelay2 = (ndelay + 1)
                nfin = (npts + ioff1)
                for i1 in range(int(ioff1), int(nfin)):
                    itemp = outbuf[i1] + outbuf[i2+ndelay]
                    #logger.debug("changing values...")
                    outbuf[i2] = outbuf[i1] + outbuf[i2+ndelay2]
                    outbuf[i1] = itemp
                    i2 += 1
    return


@jit
def comp_stats(vec, veclen, ignore):
    #//compute mean and stddev of floating point vector vec, ignoring elements in ignore != 0
    tmean = 0
    tstddev = 0
    valid_points=0

    for i in range(0, veclen):
        if ignore[i] == 0:
            tmean = tmean + vec[i]
            tstddev = tstddev + (vec[i] * vec[i])
            valid_points += 1

    tstddev = pow((tstddev - ((tmean * tmean)/valid_points))/(valid_points - 1), 0.5)
    tmean = tmean / (valid_points)

    return tmean, tstddev


@jit
def rfirej(tree_dedoppler, tdwidth, nframes, tsteps, rfiwindow, rfithresh):
    spectrum = np.zeros(tdwidth)
    spectrum_sum = np.zeros(tdwidth)
    stat_mask = np.zeros(tdwidth, dtype='uint8')
    rfi_mask = np.zeros(tdwidth, dtype='uint8')
    channelmedians = np.zeros(tdwidth)

    for i in range(0, tsteps*4):
        stat_mask[i] = 1
    for i in range(tdwidth - (tsteps*4), tsteps*4):
        stat_mask[i] = 1

    print "summing spectra\n"
    #/* sum all spectra */
    for i in range(0, nframes):
        spectrum_sum += tree_dedoppler[i * tdwidth : (i+1)*tdwidth]

    print "computing stats\n"
    mean, stddev = comp_stats(spectrum_sum, tdwidth, stat_mask)
    print "normalizing\n"
    #/* normalize */
    spectrum_sum = (spectrum_sum - mean)/stddev

    j = 0
    print "building rfimask\n"
    #/* build RFI mask based on 0 Hz/sec spectrum */
    for i in range(0, tdwidth):
        if spectrum_sum[i] > rfithresh:
            rfi_mask[i] = 1
            j += 1
    print "excluded %d channels with large 0 Hz/sec signals\n"%j

    #/* median filter */
    channelmedians= np.median(tree_dedoppler.reshape((tsteps, tdwidth))[0:nframes, :], axis=0)

    for i in range(0, nframes):
        #/* copy each frame */
        np.copyto(spectrum, tree_dedoppler[i * tdwidth: (i+1) * tdwidth])

        #/* divide every point by channel median */
        nonzero_indices = channelmedians.nonzero()
        spectrum[nonzero_indices] = spectrum[nonzero_indices]/channelmedians[nonzero_indices]

        #/* copy the spectrum back into the array */
        np.copyto(tree_dedoppler[i * tdwidth: (i+1) * tdwidth], spectrum)

    logger.debug("after rfirej, tree_dedoppler shape: %s"%str(tree_dedoppler.shape))

    return rfi_mask


#@autojit
@jit
def candsearch(spectrum, specstart, specend, candthresh, drift_rate, \
               header, fftlen, tdwidth, channel, max_val, reverse):
    logger.debug('Start searching for drift rate: %f'%drift_rate)
    j = 0
    for i in (spectrum[specstart:specend] > candthresh).nonzero()[0] + specstart:
        k =  (tdwidth - 1 - i) if reverse else i
        info_str = 'Candidate found at SNR %f! %s\t'%(spectrum[i], '(reverse)' if reverse else '')
        info_str += 'Spectrum index: %d, Drift rate: %f\t'%(i, drift_rate)
        info_str += 'Uncorrected frequency: %f\t'%chan_freq(header,  k, tdwidth, 0)
        info_str += 'Corrected frequency: %f'%chan_freq(header, k, tdwidth, 1)
        logger.info(info_str)
        j += 1
        used_id = j
        if spectrum[i] > max_val.maxsnr[k]:
            max_val.maxsnr[k] = spectrum[i]
            max_val.maxdrift[k] = drift_rate
            max_val.maxid[k] = used_id

    return j, max_val


#@jit
#@autojit
def tophitsearch(tree_dedoppler_original, max_val, tsteps, nframes, tdwidth, header, split_dir='/queencow/dedoppler_test', rfithresh = 400, logwriter=None, filewriter=None):
    maxsnr = max_val.maxsnr
    logger.debug("original matrix size: %d\t(%d, %d)"%(len(tree_dedoppler_original), tsteps, tdwidth))
    tree_orig = tree_dedoppler_original.reshape((tsteps, tdwidth))
    logger.debug("tree_orig shape: %s"%str(tree_orig.shape))
    for i in ((maxsnr > 0) & (maxsnr < rfithresh)).nonzero()[0]:
        lbound = max(0, i - tsteps)
        ubound = min(tdwidth, i + tsteps)
        skip = 0
        for j in ((maxsnr[lbound:i] > maxsnr[i]) & (maxsnr[lbound:i] < rfithresh)).nonzero()[0]:
            skip = 1
        for j in ((maxsnr[i+1:ubound] > maxsnr[i]) & (maxsnr[i+1:ubound] < rfithresh)).nonzero()[0]:
            skip = 1

        if skip:
            pass
            logger.debug("SNR not big enough... %f pass... index: %d"%(maxsnr[i], i))
        else:
            info_str = "Top hit found! SNR: %f ... index: %d"%(maxsnr[i], i)
            logger.debug(info_str)
            if logwriter:
                logwriter.info(info_str)
                #logwriter.report_tophit(max_val, i, header)
            logger.debug("slice of spectrum...size: %s"%str(tree_orig[0:nframes, lbound:ubound].shape))
            if filewriter:
                filewriter = filewriter.report_tophit(max_val, i, (lbound, ubound), tree_orig[0:nframes, lbound:ubound], header)
            else:
                np.save(split_dir + '/spec_id_%d.npy'%i, tree_orig[0:nframes, lbound:ubound])
            llbound = max(0, i - int(0.03/header['DELTAF']))
            uubound = min(tdwidth, i + int(0.03/header['DELTAF']))
            np.save(split_dir + '/spec_id_big_%d.npy'%i, tree_orig[0:nframes, llbound:uubound])
        return filewriter
