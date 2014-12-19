#!/usr/bin/env python

import sys
import logging
from dedoppler.plot_utils import *

def main():
    from optparse import OptionParser

    p = OptionParser()
    p.set_usage('plot_results <FULL_PATH_TO_TEXT_FILE>')
    p.set_description(__doc__)

    opts, args = p.parse_args(sys.argv[1:])

    if len(args)!=1:
        print 'Please specify an input file \nExiting.'
        sys.exit()
    else:
        filename = args[0]


    logging.basicConfig(format='%(relativeCreated)5d %(name)-15s %(levelname)-8s %(message)s', level = logging.DEBUG)

    plot_tophits(filename)


if __name__=='__main__':
    main()
