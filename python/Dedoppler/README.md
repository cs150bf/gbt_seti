Dedoppler 
=============================


&nbsp;
[TOC]

&nbsp;


## Summary

Based on `rawdopplersearch.c` [`gbt_seti/src/rawdopplersearch.c`](http://github.com/casper-astro/gbt_seti/blob/master/src/rawdopplersearch.c).

- Moved to Python
- Pre-calculated `drift_indexes` (no need for one extra round of `taylor_flt`)
- Matrix `compression` along either direction
- Output text file of certain formats (and some parser functions)
- ?

&nbsp;

-------------------

## Installation

```bash
$ python setup.py install
```

### Dependencies

- numba
- pyfits
- **tempo**


&nbsp;

--------------------------

## Usage


### Expected Inputs

- **FITS** files
- HDU: header, data
    - header: keywords `FCNTR`, `DELTAF`, `DELTAT`, etc.
    - data: *shape* **(fftlen, )**   &nbsp; (consistant throughout the file)
- what else...?

### Command Line

> **`$ dedoppler <FULL_PATH_TO_INPUT_FITS_FILE> [OPTIONS]`**
>
> Use `$ dedoppler -h` to view usage details.
>
> &nbsp;
> 
> Parameters:
> - `max_drift`:
> - `min_drift`:
> - `snr`:
> - `bw`:
> 
> (and so on...)



#### Example:

```bash
$ dedoppler /queencow/dedoppler_test/VOYAGER1.8448.437500_56553.034965.fits -p /queencow/dedoppler_test > dedopp.log &
```

This will take `/queencow/dedoppler_test/VOYAGER1.8448.437500_56553.034965.fits` as input (and in this particular case it will discover that this file is too big to handle all at once, so it will first partition it into smaller FITS files and save them into the directory specified by option **`-p`**, and then proceed with dedopplering for each small FITS files). Everything else was set to default values.

**TO DO**: at this moment the logging scheme is very, very poorly implemented...

#### Sample Outputs

See `/queencow/dedoppler_test/*.log`, `/queencow/dedoppler_test/*.dat` for search results and see `/queencow/dedoppler_test/*.pdf` for some plots.

&nbsp;

> MJD: 56553.034967189247	RA:   257.063400000000	DEC:    12.179900000000	DELTAT:     0.330301440000	DOPPLER:     0.000000000000
>
> `---------------------`
>
> Top Hit #1: 	Drift Rate:   0.062781	Uncorrected Frequency:    8419.278882	Corrected Frequency:    8418.606974	Index: 593508
>
> Freqs: (   8419.278107 ...    8419.279654), DELTAF:   0.00000303 (MHz)
>
> [a matrix....]
>
>
> Top Hit #2: 	Drift Rate:   0.125561	Uncorrected Frequency:    8419.278882	Corrected Frequency:    8418.606974	Index: 593508
>
> Freqs: (   8419.278107 ...    8419.279654), DELTAF:   0.00000303 (MHz)
>
> [matrix (part of the spectra)]
>
>
> Top Hit #3:...

&nbsp;

### Python Executable (script)

- `logging` configuration required
- in `dedoppler/bin` directory
- one script for searching, the other for plotting the **top hits**


> **`$ python dedopp.py  <FULL_PATH_TO_INPUT_FITS_FILE> [OPTIONS]`**
>
> **`$ ./dedopp.py  <FULL_PATH_TO_INPUT_FITS_FILE> [OPTIONS]`**

&nbsp;

> **`$ ./plot_results.py <FULL_PATH_TO_TOP_HITS_FILE>`**
>
> A bunch of \*.pdf files will show up in the same directory as the input file.
> 
> **Example:**
>
> **`$ ./plot_results.py /queencow/dedoppler_test/VOYAGER1.8448.437500_56553.034965_6.dat`**
>

&nbsp;


### Use as a package

```python
> import dedoppler
```


&nbsp;
--------------------------

