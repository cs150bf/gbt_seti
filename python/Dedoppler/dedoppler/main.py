import sys
import logging
from . import dedopp

def main():
    from optparse import OptionParser

    p = OptionParser()
    p.set_usage('dedoppler <FULL_PATH_TO_FITS_FILE>')
    p.set_description(__doc__)

    p.add_option('-m', '--min_drift', dest='min_drift', type='float', default=0.0,
            help='Set the minimum drift rate to search. Unit: Hz/sec. Default:0.0')
    p.add_option('-M', '--max_drift', dest='max_drift', type='float', default=10.0,
            help='Set the drift rate to search. Unit: Hz/sec. Default: 10.0')
    p.add_option('-s', '--snr', dest='snr', type='float', default=10.0, help='SNR threshold. Unit: ?. Default: 10.0')
    p.add_option('-b', '--bw', dest='bw', type='float', default=1, help='Specify the amount of \'compression\' to be done in frequency domain to search for more \'spread out\' signals. Unit:?. Default: ?')
    p.add_option('-r', '--rfithresh', dest='rfithresh', type='float', default=25.0, help='Specify the RFI threshold. Unit:?. Default: 25.0')
    p.add_option('-p', '--path', dest='split_dir', type='str', default='/tmp',
            help='In the case that the input FITS file size is too big to handle at once, we\'ll need to split it into smaller FITS files. This option specify where to put those FITS files. Default: /tmp ')
    p.add_option('-o', '--output', dest='out', type='str', default='', help='')
    p.add_option('-w', '--width', dest='slice_width', type='int', default=512, help='')
    p.add_option('-l', '--loglevel', dest='loglevel', type='str', default='debug', help='Specify log level')

    opts, args = p.parse_args(sys.argv[1:])

    if len(args)!=1:
        print 'Please specify a FITS file \nExiting.'
        sys.exit()
    else:
        fitsfile = args[0]

    logging.basicConfig(
        format='%(relativeCreated)5d %(name)-15s %(levelname)-8s %(message)s', level = logging.DEBUG)

    mydedopptask = dedopp.DedopplerTask(fitsfile = fitsfile, max_drift = opts.max_drift, min_drift = opts.min_drift, snr = opts.snr, bw = opts.bw, rfithresh = opts.rfithresh, split_dir = opts.split_dir)

    mydedopptask.search()




if __name__=='__main__':
    main()
