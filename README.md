# Goodman Spectrograph Exposure Time Calculator (ETC)

Python tools and configuration files to predict Signal-to-Noise of Goodman spectra

---

## Overview

This repository contains Python code to produce plots of the expected signal-to-noise for Goodman 
observations.

---

## Requirements

The code that performs the calculations requires NumPy, AstroPy, SciPy, and Matplotlib

---

## Structure

The repository contains four subdirectories:
* bin - contains the executable to run the code
* eff - contains the available Goodman efficiency curves
* sky - contains a set of sky models for various Lunar phases
* spec - contains model spectra templates of sources
* results - this directory is created by the software for the output (figures, FITS files, and tables)

---

## Usage

```text
Usage: ./bin/goodman_etc.py [options]

Options:
  -grating <value>         Select grating (e.g., 400, 600, 930, 1200)
  -mode <value>            Select the grating observing mode. The following modes are supported:
                            400: M1 or M2
                            600: UV, Blue, Mid, Red
                            930: M1, M2, M3, M4, M5, M6
                           1200: M0, M1, M2, M3, M4, M5, M6, M7
  -bin_sky <value>         Set pixel binning along the spatial axis
  -bin_disp <value>        Set pixel binning along the dispersion axis
  -slit <value>            Set the slit width in arcsec.
  -gain <slow|fast>        This sets the readout speed, gain and readout noise for the ETC.  The
                           "slow" setting corresponds to 344kHz ATTN3 readout.  The "fast" setting
                           corresponds to the 750kHz ATTN0 readout.
  -vmag <value>            The V magnitude of the source on the Vega system
  -sptype <value>          The spectral type of the object from the Pickles UVKLIB spectral library
                           (Pickles, PASP, 110, 863, 1998).  An AstroPy table of the allowable spectral
                           types is located in <path>/spec/spectab.txt
  -metal <weak|solar|rich> The metallicity of the source, as listed in the file spectab.txt
  -lp <value>              The number of days since new moon.  This is used to determine the sky
                           model.  Sky model templates have been made for the following:
                            0 - New moon (Darkest)
                            3 - Crescent (Dark)
                            7 - 1st/3rd quarter (Gray)
                           10 - Gibbous (Bright)
                           14 - Full (Brightest)
-see <value>               The seeing for the observations
-expt <value>              The exposure time in seconds
-num <value>               The number of exposures
-out <yes|no>              Create output FITS files (1-D spectra) of the expected counts and S/N.  This also
                           creates an AstroPy table with the results.  Note that the results have been converted
                           from electrons to ADU by using the gain value.  This allows for predicted spectra and
                           S/N to be compared to existing data.
```

---

## Output

Upon successful completion of running the code, the software automatically produces a PNG file that contains four subplots.
These plots are:<br> 
(Upper Left) S/N plot for 1 x ExpTime spectrum<br> 
(Upper Right) S/N plot for N x Exptime spectra<br>
(Lower Left) Predicted electrons per spectral bin for source plus sky and noise spectra for 1 x Exptime spectrum<br> 
(Lower Right) Predicted electrons per spectral bin for source plus sky and noise spectra for N x ExpTime spectra<br> 

If the "-out" argument is "yes", the software will additionally produce an AstroPy table that lists wavelength, source counts, 
sky counts, and S/N for each wavelength. The numbers in the table are determined by dividing the electrons detected by the 
gain (electrons per ADU) selected.  This was selected to facilitate comparisons with real data.  This option also write two
1-D spectra FITS files for the source and the calculated S/N per spectral bin.

---

## Contributions or comments

Contributions are welcome as are coments to improve the code to provide additional
functionality.

---

## Author

Sean Points  
NSFs National Optical-Infrared Astronomy Research Laboratory (NOIRLab)<br>
Cerro Tololo Inter-American Observatory (CTIO)</br>
sean.points@noirlab.edu
