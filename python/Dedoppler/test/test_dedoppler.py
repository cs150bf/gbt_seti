import logging
import dedoppler

fitsfilename = '/queencow/dedoppler_test/test56553.034965_6.fits'

logging.basicConfig(
        format='%(relativeCreated)5d %(name)-15s %(levelname)-8s %(message)s', 
        level = logging.DEBUG)
mydedopp = dedoppler.dedopp.Dedoppler(fitsfile=fitsfilename, max_drift = 10, min_drift = 0.0, bw = 1, snr=10, split_dir='/queencow')

mydedopp.search()

