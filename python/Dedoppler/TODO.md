TO DO List for Dedoppler
==============================

- [ ] More speed optimization (numba.jit, vectorization, numpy, etc.)
- [] Program: parameters
    - currently only half-implemented...
    - [] matrix compression along time or frequency ...
        - [x] basic functions
        - [○] how to correctly calculate new frequency resolution, etc.?
- [○] logging, logger, logs!
    - log files?
- [○] Error checking!
    - Check all class and methods
    - Particularly those having to do with file manipulation
        - [○] write permission?
        - [x] file already existed?
- [○] Output file format control
    - Text file, Numpy saved file, other format?
    - Top hits reporting... e.g. how big a sub-table to output?
- [○] tempo? tempo2?
    - really not sure about this... `tempo` deals with quite a few temporary files...
- [○] Plotting
    - The annotation needs some more work
    - drift rate?
- [○] Style
    - line wrapping
    - standarize variable names?

