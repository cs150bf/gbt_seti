import logging
import dedoppler
import numpy as np
from matplotlib.pyplot import *

fitsfilename = '/queencow/dedoppler_test/test56553.034965_6.fits'
snr = 10

logging.basicConfig(
        format='%(relativeCreated)5d %(name)-15s %(levelname)-8s %(message)s', 
        level = logging.DEBUG)

logger = logging.getLogger(__name__)

fits_handle = dedoppler.fits_wrapper.FITSHandle(filename = fitsfilename, split_dir='/queencow/dedoppler_test')

fits_obj = fits_handle.fits_list[0]
spectra, drift_indexes = fits_obj.load_data(bw_compress_width=5)

tsteps = fits_obj.tsteps
tsteps_valid = fits_obj.tsteps_valid
tdwidth = fits_obj.tdwidth
fftlen = fits_obj.fftlen
nframes = tsteps_valid

tree_dedoppler = np.zeros(tsteps * tdwidth, dtype=np.float32)

ibrev = np.zeros(tsteps, dtype=np.int32)


for i in range(0, tsteps):
    ibrev[i] = dedoppler.dedopp.bitrev(i, int(np.log2(tsteps)))



for i in range(0, nframes):
    sind = tdwidth*i + tsteps*4
    cplen = fftlen
    np.copyto(tree_dedoppler[sind: sind + cplen], spectra[i])
    sind = i * tdwidth
    np.copyto(tree_dedoppler[sind: sind + tsteps*4], spectra[i, fftlen-(tsteps*4):fftlen])


if not (tree_dedoppler.reshape((tsteps, tdwidth))[0:tsteps_valid, tsteps*4:tdwidth-tsteps*4] == spectra).all():
    logger.error("Something went wrong during the array copying!")

tree_dedoppler_flip = np.empty_like(tree_dedoppler)
np.copyto(tree_dedoppler_flip, tree_dedoppler)

rfi_mask = dedoppler.dedopp.rfirej(tree_dedoppler, tdwidth, nframes, tsteps, 2, 100)

dedoppler.dedopp.taylor_flt(tree_dedoppler, tsteps*tdwidth, tsteps)
maxsnr = np.zeros(tdwidth, dtype=np.float32)


stat_mask = np.zeros(tdwidth, dtype='uint8')
for i in range(0, tsteps*4):
    stat_mask[i] = 1
for i in range(tdwidth - tsteps*4, tsteps*4):
    stat_mask[i] = 1

k = 1
indx = ibrev[drift_indexes[k]]*tdwidth
spectrum = tree_dedoppler[indx : indx+tdwidth]
mean_val, stddev = dedoppler.dedopp.comp_stats(spectrum, tdwidth, stat_mask)

for i in range(0, tdwidth):
    spectrum[i] = (spectrum[i] - mean_val)/stddev

specstart = tsteps*4
specend = tdwidth - tsteps*6

drift_rate = k * fits_obj.drift_rate_resolution

logger.info("drift rate: %f"%drift_rate)

for i in range(specstart, specend):
    if spectrum[i] > snr:
       print "Index: %d\tSNR: %f"%(i, spectrum[i])

'''
for i in range(0, 32):
    plot(range(32320*i, 32320*(i+1)), spectrum[32320*i: 32320*(i+1)])
    title('%d-th'%i)
    show()
'''

logger.info("Doppler correcting reverse...")
dedoppler.dedopp.FlipX(tree_dedoppler_flip, tdwidth, tsteps)
dedoppler.dedopp.taylor_flt(tree_dedoppler_flip, tdwidth, tsteps)

drift_rate = -1 * k * fits_obj.drift_rate_resolution

logger.info("drift rate: %f"%drift_rate)

for i in range(specstart, specend):
    if spectrum[i] > snr:
       print "Index: %d\tSNR: %f"%(i, spectrum[i])


