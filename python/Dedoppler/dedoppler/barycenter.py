import struct
from subprocess import Popen, PIPE, call
import numpy as np
import math
import tempfile

import logging
logger = logging.getLogger(__name__)

SECPERDAY = 86400.0

def doppler(freq_observed, voverc):
    '''
    /* This routine returns the frequency emitted by a pulsar */
    /* (in MHz) given that we observe the pulsar at frequency */
    /* freq_observed (MHz) while moving with radial velocity */
    /* (in units of v/c) of voverc wrt the pulsar. */
    '''
    return freq_observed * (1.0 + voverc)


def read_resid_rec(outfile, firsttime=1, use_ints=0):
    '''
    /* This routine reads a single record (i.e. 1 TOA) from */
    /* the file resid2.tmp which is written by TEMPO.       */
    /* It returns 1 if successful, 0 if unsuccessful.       */
    '''

    #// The default Fortran binary block marker has changed
    #// several times in recent versions of g77 and gfortran.
    #// g77 used 4 bytes, gfortran 4.0 and 4.1 used 8 bytes
    #// and gfortrans 4.2 and higher use 4 bytes again.
    #// So here we try to auto-detect what is going on.
    #// The current version should be OK on 32- and 64-bit systems

    if firsttime:
        llraw = outfile.read(8)
        ddraw = outfile.read(8)
        #print llraw, ddraw
        ll = struct.unpack('<q', llraw)[0]
        dd = struct.unpack('<d', ddraw)[0]
        if 1:
            print "(long long) index = %d  (MJD = %17.10f)\n"%(ll, dd)
        if (ll != 72 or dd < 40000.0 or dd > 70000.0):
            outfile.seek(0)
            ii = struct.unpack('<i', outfile.read(4))[0]
            dd = struct.unpack('<d', outfile.read(8))[0]
            if True:
                print "(int) index = %d    (MJD = %17.10f)\n"%(ii, dd)
            if (ii == 72 and (dd > 40000.0 and dd < 70000.0)):
                use_ints = 1
            else:
                print "\nError:  Can't read the TEMPO residuals correctly!\n"
                return None, None, None
        outfile.seek(0)
        firsttime = 0

    if use_ints:
        print "Using ints..."
        ii = struct.unpack('<i', outfile.read(4))[0]
    else:
        ll = struct.unpack('<q', outfile.read(8))[0]
    #//  Now read the rest of the binary record
    d = struct.unpack('<9d', outfile.read(8*9))
    if 1: 
        print "Barycentric TOA = %17.10f\n"%d[0]
        print "Postfit residual (pulse phase) = %f\n"%d[1]
        print "Postfit residual (seconds) = %f\n"%d[2]
        print "Orbital phase = %f\n"%d[3]
        print "Barycentric Observing freq = %f\n"%d[4]
        print "Weight of point in the fit = %f\n"%d[5]
        print "Timing uncertainty = %g\n"%d[6]
        print "Prefit residual (seconds) = %f\n"%d[7]
        print "??? = %f\n\n"%d[8]
    toa = d[0]
    obsf = d[4]
    if use_ints:
        ii = struct.unpack('i', outfile.read(4))[0]
        return ii, toa, obsf
    else:
        ll = struct.unpack('q', outfile.read(8))[0]
        return ll, toa, obsf



def correct(header, obs_length): #topotimes, N, ra_str, dec_str, obs, ephem):
    '''
    /* This routine uses TEMPO to correct a vector of           */
    /* topocentric times (in *topotimes) to barycentric times   */
    /* (in *barytimes) assuming an infinite observation         */
    /* frequency.  The routine also returns values for the      */
    /* radial velocity of the observation site (in units of     */
    /* v/c) at the barycentric times.  All three vectors must   */
    /* be initialized prior to calling.  The vector length for  */
    /* all the vectors is 'N' points.  The RA and DEC (J2000)   */
    /* of the observed object are passed as strings in the      */
    /* following format: "hh:mm:ss.ssss" for RA and             */
    /* "dd:mm:ss.ssss" for DEC.  The observatory site is passed */
    /* as a 2 letter ITOA code.  This observatory code must be  */
    /* found in obsys.dat (in the TEMPO paths).  The ephemeris  */
    /* is either "DE200" or "DE405".                            */
    '''
    topotimes =  [header['MJD'], header['MJD'] + obs_length/SECPERDAY]
    N = len(topotimes)
    ra_str = dec2hms(header['RA'])
    dec_str = dec2hms(header['DEC'])
    obs = "GB"
    ephem = "DE405"
    fobs = 1000.0
    command = np.zeros(100, dtype='uint8')
    barytimes = np.zeros(N)
    voverc = np.zeros(N)

    print 'topotimes: ', topotimes

    #/* Write the free format TEMPO file to begin barycentering */

    tempdir = '.'
    #tempdir = tempfile.mkdtemp(prefix='barycenter')
    outfile = tempfile.NamedTemporaryFile(prefix='bary', 
            suffix='.tmp', dir=tempdir, delete=False)
    barytmp_path = outfile.name
    outfile.write("C  Header Section\n"
                    + "  HEAD                    \n"
                    + "  PSR                 bary\n"
                    + "  NPRNT                  2\n"
                    + "  P0                   1.0 1\n"
                    + "  P1                   0.0\n"
                    + "  CLK            UTC(NIST)\n"
                    + "  PEPOCH           %19.13f\n"%(topotimes[0])
                    + "  COORD              J2000\n"
                    + "  RA                    %s\n"%(ra_str)
                    + "  DEC                   %s\n"%(dec_str)
                    + "  DM                   0.0\n"
                    + "  EPHEM                 %s\n"%(ephem)
                    + "C  TOA Section (uses ITAO Format)\n"
                    + "C  First 8 columns must have + or -!\n"
                    + "  TOA\n")

    #/* Write the TOAs for infinite frequencies */

    for i in range(0, N):
        outfile.write("topocen+ %19.13f  0.00     0.0000  0.000000  %s\n"%(topotimes[i], obs))
    outfile.write("topocen+ %19.13f  0.00     0.0000  0.000000  %s\n"%(
                    topotimes[N - 1] + 10.0 / SECPERDAY, obs))
    outfile.write("topocen+ %19.13f  0.00     0.0000  0.000000  %s\n"%(
                    topotimes[N - 1] + 20.0 / SECPERDAY, obs))
    outfile.close()

    #/* Call TEMPO */

    #/* Check the TEMPO *.tmp and *.lis files for errors when done. */

    command = "tempo " + barytmp_path +" > " + tempdir + "/tempoout_times.tmp"
    print command
    if call(command, shell=True) ==-1:
        print "\nError calling TEMPO in barycenter.c!\n"
        return None, None

    #/* Now read the TEMPO results */

    temporaryfile = "resid2.tmp"
    outfile = open(temporaryfile, "rb")

    #/* Read the barycentric TOAs for infinite frequencies */

    for i in range(0, N):
        tmp, barytimes[i], obsf = read_resid_rec(outfile)
    outfile.close()

    '''
    /* rename("itoa.out", "itoa1.out"); */
    /* rename("bary.tmp", "bary1.tmp"); */
    /* rename("bary.par", "bary1.par"); */

    /* Write the free format TEMPO file to begin barycentering */
    '''
    outfile = tempfile.NamedTemporaryFile(prefix='bary', suffix='.tmp', 
            dir=tempdir, delete=False)
    barytmp_path = outfile.name
    outfile.write("C  Header Section\n"
                    + "  HEAD                    \n"
                    + "  PSR                 bary\n"
                    + "  NPRNT                  2\n"
                    + "  P0                   1.0 1\n"
                    + "  P1                   0.0\n"
                    + "  CLK            UTC(NIST)\n"
                    + "  PEPOCH           %19.13f\n"%(topotimes[0])
                    + "  COORD              J2000\n"
                    + "  RA                    %s\n"%(ra_str)
                    + "  DEC                   %s\n"%(dec_str)
                    + "  DM                   0.0\n"
                    + "  EPHEM                 %s\n"%(ephem)
                    + "C  TOA Section (uses ITAO Format)\n"
                    + "C  First 8 columns must have + or -!\n"
                    + "  TOA\n")

    #/* Write the TOAs for finite frequencies */

    for i in range(0, N):
        outfile.write("topocen+ %19.13f  0.00  %9.4f  0.000000  %s\n"%(topotimes[i], fobs, obs))
    outfile.write("topocen+ %19.13f  0.00  %9.4f  0.000000  %s\n"%(
                   topotimes[N - 1] + 10.0 / SECPERDAY, fobs, obs))
    outfile.write("topocen+ %19.13f  0.00  %9.4f  0.000000  %s\n"%(
                    topotimes[N - 1] + 20.0 / SECPERDAY, fobs, obs))
    outfile.close()

    #/* Call TEMPO */

    #/* Insure you check the file tempoout.tmp for  */
    #/* errors from TEMPO when complete.            */

    command = "tempo " + barytmp_path +" > " + tempdir + "/tempoout_vels.tmp"
    if call(command, shell=True) == -1:
        print "\nError calling TEMPO in barycenter.c!\n"
        return None, None


    #/* Now read the TEMPO results */
    temporaryfile = "resid2.tmp"
    outfile = open(temporaryfile, "rb")

    #/* Determine the radial velocities using the emitted freq */
    for i in range(0, N):
        tmp, tmp2, femit = read_resid_rec(outfile)
        voverc[i] = femit / fobs - 1.0
    outfile.close()

    #/* Cleanup the temp files */ 

    call("rm tempo.lis", shell=True)
    call("rm tempoout_times.tmp", shell=True)
    call("rm tempoout_vels.tmp", shell=True)
    call("rm resid2.tmp", shell=True)
    #call("rm bary.tmp", shell=True)
    call("rm matrix.tmp", shell=True)
    call("rm bary.par", shell=True)
    call("rm bary*.tmp", shell=True)

    #return barytimes, voverc
    #header['barytimes'] = barytimes
    #header['voverc'] = voverc

    header['baryv'] = voverc[0]
    header['barya'] = (voverc[0] - voverc[1])/obs_length
    return header

# edited from test_psrfits.c
def dec2hms(in_double, sflag=False):
    sign = 1;
    if in_double < 0.0:
        sign = -1
        in_double = -in_double
    h = int(math.floor(in_double))
    in_double = in_double - h
    in_double *= 60.0
    m = int(math.floor(in_double))
    in_double -= m
    in_double *= 60.0
    s = in_double;
    if sign==1 and sflag:
        ptr='+'
    elif sign==-1:
        ptr='-'
    else:
        ptr=''
    ptr = ptr + "%2.2d:%2.2d:%07.4f"%(h, m, s)
    return ptr

