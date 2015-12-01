# optics_photochem
Various scripts for reading in, parsing, and analyzing marine optical and photochemical data. Some dependencies may reside in the repo https://github.com/jamesrco/dependencies-useful-scripts

Currently contains:

1. [JAZ_data_read.m](https://github.com/jamesrco/optics_photochem/blob/master/JAZ_data_read.m) and dependencies, a MATLAB script to read in and store to file UV-VIS spectra from an Ocean Optics JAZ device. Includes some code for basic analysis as well, such as calculation of integrated photon fluxes for spectral bands of photochemical interest such as UVB, UVA, PAR, etc. You'll have to configure the script and [JAZ_wavelengths.csv](https://github.com/jamesrco/optics_photochem/blob/master/JAZ_wavelengths.csv) for your own JAZ device if it doesn't have a CCD with the same specifications as the one for which the script was written (200-850 nm).
