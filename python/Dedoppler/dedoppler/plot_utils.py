import numpy as np
import re
from subprocess import Popen, PIPE
from threading import Thread
from datetime import datetime
import copy
from glob import glob
import gc
import os
from matplotlib.pyplot import *
from matplotlib import gridspec
from matplotlib import ticker
from matplotlib.offsetbox import TextArea, AnnotationBbox
from scipy import stats
import scipy

import logging
logger = logging.getLogger(__name__)

HFONT = {'fontname':'Helvetica'}
TFONT = {'fontname': 'Times'}

def load_tophitsheader(filehandle):
    header = None
    regex0 = r"MJD:\s+(?P<MJD>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\s+RA:\s+(?P<RA>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\s+DEC:\s+(?P<DEC>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\s+DELTAT:\s+(?P<DELTAT>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\s+DOPPLER: (?P<DOPPLER>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)?"
    r0 = re.compile(regex0)
    while True:
        aline = filehandle.readline()
        if not aline:
            break
        elif aline[0]=='\n':
            pass
        elif aline[0]=='-':
            pass
        elif aline[0]=='M':
            amatch = r0.match(aline)
            if amatch:
                header = amatch.groupdict()
                for key in header.keys():
                    header[key] = np.float32(header[key])
                logger.debug(aline)
                logger.info(str(header))
            return header
    return header


def load_tophits(filehandle):
    regexp1 = r"Top Hit #(?P<top_hit>[-+]?\d+):\s+Drift Rate:\s+(?P<drift_rate>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\s+Uncorrected Frequency:\s+(?P<ufreq>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\s+Corrected Frequency:\s+(?P<cfreq>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\s+Index:\s+(?P<index>[-+]?\d+)?"
    #dtype1 =[('tophit', np.int), ('driftrate', np.float32), ('ufreq', np.float32), ('cfreq', np.float32), ('index', np.int32)]
    regexp2 = r"Freqs:\s+\(\s+(?P<ufreqs_start>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\s+\.\.\.\s+(?P<ufreqs_end>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\),\s+DELTAF:\s+(?P<DELTAF>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)?"
    #dtype2 =  [('cfreqstart', np.float32), ('cfreqend', np.float32), ('deltaf', np.float32)]
    r1 = re.compile(regexp1)
    r2 = re.compile(regexp2)

    opts1 = None
    opts2 = None
    specs = None
    while True:
        aline = filehandle.readline()
        if not aline: break
        if aline[0] == '\n':
            if opts1 and opts2: break
            pass
        elif aline[0:3]=='---':
            pass
        elif aline[0:3]=='Top':
            #opts1 = np.fromregex(aline, regexp1, dtype=dtype1)
            opts1 = r1.match(aline).groupdict()
            for key in opts1.keys():
                opts1[key] = np.float32(opts1[key])
            logger.debug(aline)
            logger.info(str(opts1))
        elif opts1 and aline[0:5]=='Freqs':
            #opts2 = np.fromregex(aline, regexp2, dtype=dtype2)
            opts2 = r2.match(aline).groupdict()
            for key in opts2.keys():
                opts2[key] = np.float32(opts2[key])
            logger.debug(aline)
            logger.info(str(opts2))
        elif opts1 and opts2 and (not aline[0] in 'TF'):
            specs_row = np.fromstring(aline, dtype=np.float32, sep=' ')
            if specs is None: 
                specs = specs_row
            else:
                specs = np.vstack((specs, specs_row))
        else:
            pass
    return opts1, opts2, specs


def threshold(spec, threshmin = 1.0):
    nspec = np.zeros(spec.shape)
    #lspec = np.zeros(spec.shape)
    nspec = stats.threshold(spec, threshmin = threshmin, newval = threshmin)
    #lspec = np.log10(nspec) - np.log10(threshmin)
    return nspec


def logscale(spec, threshmin = 1.0):
    lspec = np.zeros(spec.shape)
    lspec = np.log10(spec) - np.log10(threshmin)
    return lspec


def integrate(specs):
    ispec = np.zeros(specs.shape[-1])
    nspecs = specs.shape[0]
    for i in range(0, nspecs):
        ispec += specs[i]
    return ispec


def myim(specs, freqs, ispec, header, dec_rate, vmin=0, vmax=None):
    global TFONT
    global HFONT
    DELTAT = header['DELTAT']
    if vmax is None:
        vmax = np.median(specs)*2
    figure(figsize=(32, 14), dpi=200)
    matplotlib.rcParams.update({'font.size': 32, 
        'font.family': 'serif', 
        'font.serif': ['Palatino', 'Times', 'New Century Schoolbook']})
    gs = gridspec.GridSpec(2, 1, height_ratios = [3, 1])

    ax1 = subplot(gs[0])
    ax1.imshow(specs, aspect='auto', vmin=vmin, vmax=vmax, cmap='binary')
    textstr = 'Telescope: GBT\nMJD: %f\nRA: %f\nDEC: %f\n'%(header['MJD'], 
               header['RA'], header['DEC']) 
    textstr += r'$F_{cntr}$' + ': %f (Uncorredcted)\n'%(header['ufreq'])
    textstr += r'$F_{cntr}$' + ': %f (Corredcted)\n'%(header['cfreq'])
    #textstr += 'Doppler: %f\n'%header['DOPPLER']
    textstr += 'Decimation rate (plotting): %d'%dec_rate
    #ax1.text(8, 16, textstr, bbox={'alpha':0.5, 'pad':10})
    offsetbox = TextArea(textstr, minimumdescent=False)
    ab = AnnotationBbox(offsetbox, (0.3, 0.95), xycoords='axes fraction')              
    ax1.add_artist(ab)
    nspecs = specs.shape[0]
    ytick_step = max(2, int(round(nspecs*header['DELTAT']/25))*5)
    yts = np.arange(0, specs.shape[0], ytick_step/DELTAT)
    ytlabels = np.arange(0, len(yts))*ytick_step
    ax1.yaxis.set_ticks(yts)
    ax1.yaxis.set_ticklabels(ytlabels)
    ax1.yaxis.grid(color='0.8', linestyle='dashed')
    nfreqs = specs.shape[-1]
    nm_ind = nfreqs/2
    if nm_ind == int(nm_ind):
        cfreq = (freqs[nm_ind-1] + freqs[nm_ind])/2
    else:
        cfreq = freqs[(nfreqs-1)/2]
    xtick_step = int(max(2, int(round(1e6*nfreqs*header['DELTAF']/700)*100)/(1e6*header['DELTAF'])))
    nx = int(nm_ind/xtick_step)
    nm_ind = int(nm_ind)
    xts = range(nm_ind - nx*xtick_step, nfreqs, xtick_step)
    xtscale = int(-(np.log10((freqs[xts[int(len(xts)/2)+1]] - cfreq))))+1
    xtlabels = [ "%.0f" % f for f in (freqs[xts] - cfreq)*(10**xtscale)] 
    ax1.xaxis.set_ticks(xts)
    ax1.xaxis.set_ticklabels([])
    ylabel('Time (seconds)', **TFONT)
    print 'yticklabels', yts, ytlabels

    #subplot(2, 1, 2)
    #ax2 = subplot2grid((4,4), (3,0), colspan=4)
    ax2 = subplot(gs[1])
    subplots_adjust(left=0.1, bottom=0.1, hspace=0)
    ax2.plot(ispec)
    mi = min(ispec)
    ma = max(ispec)
    adj = (ma - mi)/4
    ylim([mi-adj, ma+adj])
    xlim([0, len(ispec)])
    ax2.xaxis.set_ticks(np.array(xts))
    ax2.xaxis.set_ticklabels(xtlabels)
    #ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    print 'xticklabels', xts, xtlabels, xtscale
    ax2.yaxis.set_ticks([])
    ax2.yaxis.set_ticklabels([])
    xlabel('Frequency (' + r'$F_{cntr}$' + '+/- Hz)', **TFONT)
    annotate(r'$\times10^%d$'%(6-xtscale), xy=(0.92, 0.05), xycoords='figure fraction', fontsize=35)
    suptitle('VOYAGER 1, drift rate: %f'%header['drift_rate'], fontsize=32, **TFONT)
    #show()

def myimshow(specs, freqs, ispec, header, dec_rate, vmin=0, vmax=None):
    myim(specs, freqs, ispec, header, dec_rate, vmin=vmin, vmax=vmax)
    show()

def myimsave(specs, freqs, ispec, header, dec_rate, vmin=0, vmax=None, filename='test.pdf'):
    myim(specs, freqs, ispec, header, dec_rate, vmin=vmin, vmax=vmax)
    savefig(filename)
    #show()

def plot_tophits(filename):
    fn_stem = filename[:filename.rfind('.dat')]
    dec_rate = 1
    with open(filename, 'r') as f:
        header = load_tophitsheader(f)
        DELTAT = header['DELTAT']

        while True:
            opts1, opts2, specs = load_tophits(f)
            #header['DELTAF'] = opts2['deltaf']
            header.update(opts1)
            header.update(opts2)
            nframes = specs.shape[0]
            xtick_step = int(1e-2/(header['DELTAF']*dec_rate))
            ytick_step = max(2, int(round(nframes*header['DELTAT']/25))*5)

            ofreqs = np.arange(header['ufreqs_start'], header['ufreqs_end'], header['DELTAF'])
            freqs = ofreqs[0::dec_rate]
            if len(freqs) != specs.shape[-1]:
                logger.error("Frequency array dimension doesn\t match that of the specs!")
                freqs = np.arange(header['ufreq']- 0.5*header['DELTAF']*specs.shape[-1], header['ufreq'] + 0.5*header['DELTAF']*specs.shape[-1], header['DELTAF'])

            logger.debug('Integrating spectra...')
            ispec = integrate(specs)
            #vmin = 0.15
            #vmax = 1.0
            vmin = np.min(specs)
            vmax = np.max(specs)
            #myimshow(specs, freqs, ispec, header, dec_rate, vmin=vmin, vmax=vmax)
            filename = fn_stem + '_tophit%d.pdf'%(opts1['top_hit'])
            myimsave(specs, freqs, ispec, header, dec_rate, vmin=vmin, vmax=vmax, filename = filename)

            #np.save(fn_stem + '_tophit%d'%(opts1['top_hit']), specs)

