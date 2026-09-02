#!/usr/bin/env python
# coding: utf-8

"""Calculate S/N for Goodman spectra

Cerro Tololo Inter-American Observatory

Synopsis
--------

Read user input to calculate signal-to-noise for Goodman observations

Command Line Usage
------------------

::

    Basic usage::
    
    goodman_etc.py -grating [400,600,930,1200] 
    -mode [M0...M7; UV, Blue, Mid, Red] 
    -bin_sky [1,2] 
    -bin_disp [1,2] 
    -slit [0.45, 0.6, 0.8, 1.0, 1.2] 
    -gain [low,high] 
    -vmag [V_mag] 
    -sptype [Spectral Type] 
    -metal [weak, solar, rich]
    -ld [0,3,7,10,14] 
    -see [Seeing]
    -expt [Exp Time]

    Options
    ^^^^^^^

    -grating 
    Choose grating for observations.  The currently supported gratings are
    400, 600, 930, and 1200

    -mode 
    Select the grating mode. The following modes are supported:
    400: M1 or M2
    600: UV, Blue, Mid, Red
    930: M1, M2, M3, M4, M5, M6
    1200: M0, M1, M2, M3, M4, M5< M6, M7

    -bin_sky
    Pixel binning along the spatial axis

    -bin_disp
    Pixel binning along the dispersion axis

    -slit
    The slit width in arcsec.  You should choose one of the installed 
    long slits.

    -gain
    The gain setting.  This can be "slow" or "fast".  The "slow" setting
    corresponds to 344kHz ATTN readout (Gain 1.48e/ADU; RdNoise 3.89e) 
    and the "fast" setting corresponds to the 750kHz ATTN0 readout (Gain
    3.77; RdNoise 8.99e)

    -vmag
    The V magnitude of the source on the Vega system

    -sptype
    The spectral type of the object for Solar metallicity

    -metal
    The metallicity used by the Pickels templates

    -lp
    The lunar day.  This is used to determine the sky model.
    0 - New moon (Darkest)
    3 _ Crescent (Dark)
    7 - 1st/3rd quarter (Gray)
    10 - Gibbous (Bright)
    14 - Full (Brightest)

    -see 
    The seeing for the observations

    -expt
    The exposure time in seconds

Description
-----------

goodman_etc.py is written to provide basic signal-to-noise for the Goodman
Red camera.  

Version History
---------------

2026-08-24 sdp
    Coding begun

Examples
--------

"""
import os
from astropy.table import Table
from astropy.io import fits
import astropy.units as u
from astropy.constants import h, c
import matplotlib.pyplot as plt
import warnings
from astropy.units import UnitsWarning
from scipy.special import erf
from scipy.interpolate import interp1d
import numpy as np


# ------------------------------------------------------------
# Ignore warnings
# ------------------------------------------------------------
warnings.simplefilter("ignore", UnitsWarning)

# ------------------------------------------------------------
# Setup directories
# ------------------------------------------------------------
exec_dir = os.getcwd()
SPECDIR = exec_dir + '/spec/'
SKYDIR = exec_dir + '/sky/'
EFFDIR = exec_dir + '/eff/'

# ------------------------------------------------------------
# Telescope data
# ------------------------------------------------------------

# Telescope diameter [cm]
diam = 4.1 * 100.0

# Central obstruction [cm]
central_obs = 0.9365 * 100.0

# Telescope collecting area [cm^2]
sa = (
    np.pi * (diam / 2.0)**2
    - np.pi * (central_obs / 2.0)**2
)

area = sa * u.cm**2

#print("Telescope collecting area =", sa, "cm^2")

pix_scale = 0.15 * u.arcsec / u.pixel

# ============================================================
# Get gain and readut noise
# ============================================================
def get_gain(readout):

    if readout == "slow":
        gain = 1.48
        readnoise = 3.89

    if readout == "fast":
        gain = 3.77
        readnoise = 8.99

    return gain, readnoise

# ============================================================
# Define slit loss function
# ============================================================

def slit_loss(slit, seeing):

    Trans = erf(
        (slit / 2.0) /
        ((2.0**0.5) * (seeing / 2.35482))
    )

    Loss = 1.0 - Trans

    return Trans, Loss

# ============================================================
# Find spectral template from Pickles library
# ============================================================

def get_spec_template(spec_type, metallicity):
    spec_type = spec_type.lower()
    metallicity = metallicity.lower()
    speclist = SPECDIR + 'spectab.txt'
    spectab = Table.read(speclist, format='ascii.fixed_width_two_line')
    spectab["SpType"] = [str(x).lower() for x in spectab["SpType"]]
    spectab["Metal"] = [str(x).lower() for x in spectab["Metal"]]

    match = (spectab["SpType"] == spec_type) & (spectab["Metal"] == metallicity)

    cap_sptype = spec_type.upper()

    if not any(match):
        filename = "none.txt"
    else:
        filename = spectab["Specfile"][match][0]

    return filename

# ============================================================
# Get Goodman throughput efficiency
# ============================================================

def get_efficiency(grism, obsmode):
    obsmode = obsmode.lower()
    eff_file = f"eff_{grism}{obsmode}.txt"
    #print(eff_file)

    return eff_file

# ============================================================
# Add a leeading 0 to an integer string
# ============================================================

def add_leading_zero_fixed_length(integer_value, length):
    # Convert the integer to a string and prepend zeros to achieve the fixed length
    str_value = f"{integer_value:0{length}}"
    return str_value

# ============================================================
# Get sky model
# ============================================================

def get_sky_model(moon_phase):

    temp_moon_phase = int(moon_phase)
    if temp_moon_phase < 2:
        print("Requested lunar phase:", temp_moon_phase)
        moon_phase = 0
        print("Using lunar phase:", moon_phase)
    elif temp_moon_phase >= 2 and temp_moon_phase < 5:
        print("Requested lunar phase:", temp_moon_phase)
        moon_phase = 3
        print("Using lunar phase:", moon_phase)
    elif temp_moon_phase >= 5 and temp_moon_phase < 9:
        print("Requested lunar phase:", temp_moon_phase)
        moon_phase = 7
        print("Using lunar phase:", moon_phase)
    elif temp_moon_phase >=9 and temp_moon_phase < 12:
        print("Requested lunar phase:", temp_moon_phase)
        moon_phase = 10
        print("Using lunar phase:", moon_phase)
    elif temp_moon_phase >=12:
        print("Requested lunar phase:", temp_moon_phase)
        moon_phase = 14
        print("Using lunar phase:", moon_phase)
    lp = str(moon_phase)
    llp = len(lp)
    if llp < 2:
       moon_phase = add_leading_zero_fixed_length(moon_phase, 2)

    temp_sky = f"sky_lp{moon_phase}_am12.fits"
    sky_model = SKYDIR + temp_sky

    return moon_phase, sky_model

# ============================================================
# Calculate dispersion per pixel
# ============================================================

def get_dispersion(grism, bin_dispersion):

    if grism == '400':
        disp = 1.0 * bin_dispersion
    elif grism == '600':
        disp = 0.65 * bin_dispersion
    elif grism == '930':
        disp = 0.42 * bin_dispersion
    elif grism == '1200':
        disp = 0.31 * bin_dispersion
    else:
        disp = 0.0

    return disp

# ============================================================
# Read efficiency curve
# ============================================================

def read_efficiency_curve(eff_file):

    eff_tab = Table.read(
        eff_file,
        format='ascii.fixed_width_two_line'
    )

    eff_tab["WAVE"].unit = u.Angstrom
    wave_eff = eff_tab["WAVE"].quantity
    eff = np.asarray(
        eff_tab["EFF"],
        dtype=float
    )

    return wave_eff, eff

# ============================================================
# Get wavelength range
# ============================================================

def get_wavelength_range(wave_eff):

    wmin = np.min(wave_eff)
    wmax = np.max(wave_eff)

    return wmin, wmax
    
# ============================================================
# Make source spectrum
# ============================================================

def make_source_spectrum(spec_file, vmag, eff_file, trans, disp, exptime, sa):
    src_tab = Table.read(spec_file)
    src_tab["WAVELENGTH"].unit = u.Angstrom
    src_tab["FLUX"].unit = (
        u.erg
        / u.cm**2
        / u.s
        / u.Angstrom
    )
    #
    # Pickels magnitudes scaled to V=0
    # 
    vref = 0
    flux_ratio = 10**(0.4 * (vmag - vref))

    src_tab["NFLUX"] = (
        src_tab["FLUX"] / flux_ratio
    )

    wave_eff, eff = read_efficiency_curve(eff_file)
    #print(wave_eff)
    #print(eff)

    wmin, wmax = get_wavelength_range(wave_eff)
    print("Wave_min", wmin)
    print("Wave_max", wmax)

    
    wave_src = src_tab["WAVELENGTH"].quantity
    flux_src = src_tab["NFLUX"].quantity
    spec_mask = (
        (wave_src >= wmin)
        & (wave_src <= wmax)
    )

    wave_src_obsmode = wave_src[spec_mask]
    flux_src_obsmode = flux_src[spec_mask]

    f_interp = interp1d(
        wave_src_obsmode.to_value(u.Angstrom),
        flux_src_obsmode.to_value(
            u.erg
            / u.cm**2
            / u.s
            / u.Angstrom
        ),
        kind='linear',
        bounds_error=True
    )

    flux_resampled = (
        f_interp(
            wave_eff.to_value(u.Angstrom)
        )
        * u.erg
        / u.cm**2
        / u.s
        / u.Angstrom
    )

    #print(wave_eff)
    #print(flux_resampled)

    # Wavelength bin width
    delta_lambda = disp * u.Angstrom
    #print(delta_lambda)

    # Flux in each dispersion bin
    flux_per_bin = (
        flux_resampled 
        * delta_lambda
    )

    # Photon energy
    photon_energy = (
       h * c / wave_eff
    ).to(u.erg)

    # Photon flux
    # photons / cm^2 / s / bin (dispersion per pixel)
    photon_flux = (
       flux_per_bin 
       / photon_energy
    ).value * u.ph / (u.s * u.cm**2)
    print(photon_flux[0])

    # Apply slit transmission
    photon_flux_slit = (
        photon_flux
        * trans
    )

    # Apply total efficiency
    detected_photon_flux = (
        photon_flux_slit
        * eff
    )

    # Apply telescope collecting area
    # photons / s / bin
    detected_photon_rate = (
        detected_photon_flux
        * area
    )
    print(detected_photon_rate[0])

    # Total photons detected
    # photons / bin
    detected_photons = (
        detected_photon_rate
        * exptime
        * u.s
    )
    print(detected_photons[0])
    

    return wave_eff, detected_photons

# ============================================================
# Make noise spectrum
# ============================================================

def make_noise_spectrum(source_spec, sky_model, eff_file, rdnoise, disp, exptime, sa, pix_scale, bin_dispersion, bin_spatial, seeing):

    with fits.open(sky_model, mode="update") as hdul:

        for hdu in hdul:
            if not isinstance(hdu, fits.BinTableHDU):
                continue

            for i in range(1, len(hdu.columns) + 1):

                key = f"TUNIT{i}"
                if key not in hdu.header:
                    continue
                unit = hdu.header[key]

                if unit == "ph/s/m2/micron/arcsec2":
                    hdu.header[key] = "ph/s/m2/um/arcsec2"
                elif unit == "1":
                    del hdu.header[key]
        hdul.flush()

    sky_tab = Table.read(sky_model)

    pixel_scale = bin_spatial * pix_scale

    flux=sky_tab["flux"]
    wave=sky_tab["lam"]

    # Convert nm to Angstrom
    new_wave = wave * 10
    # Convert Ph/m^2/s/micron/arcsec^2 to Ph/cm^2/s/A/pix
    new_flux = flux * 1e-8 * pixel_scale**2

    sky_tab["WAVELENGTH"] = new_wave
    sky_tab["NFLUX"] = new_flux

    sky_tab["WAVELENGTH"].unit = u.Angstrom
    sky_tab["NFLUX"].unit = (
        u.ph
        / u.cm**2
        / u.s
        / u.Angstrom
        / u.pixel
    )

    see = seeing*u.arcsec
    # Get number of pixels used for sky extraction and read noise
    npix = 2 * (see / pixel_scale)
    npix = npix.to_value()
    #print(npix)

    # 
    wave_eff, eff = read_efficiency_curve(eff_file)
    wmin, wmax = get_wavelength_range(wave_eff)
    
    wave_sky = sky_tab["WAVELENGTH"].quantity
    flux_sky = sky_tab["NFLUX"].quantity
    sky_mask = (
        (wave_sky >= wmin)
        & (wave_sky <= wmax)
    )

    wave_sky_obsmode = wave_sky[sky_mask]
    flux_sky_obsmode = flux_sky[sky_mask]

    s_interp = interp1d(
        wave_sky_obsmode.to_value(u.Angstrom),
        flux_sky_obsmode.to_value(
            u.ph
            / u.cm**2
            / u.s
            / u.Angstrom
            / u.pixel
        ),
        kind='linear',
        bounds_error=True
    )

    sky_resampled = (
        s_interp(wave_eff.to_value(u.Angstrom)
        )
        * u.ph
        / u.cm**2
        / u.s
        / u.Angstrom
        / u.pixel
    )

    delta_lambda = disp * u.Angstrom

    sky_per_bin = (
        sky_resampled
        * delta_lambda
    )

    # Sky already in photons 
    # Don't convert to energy

    # Sky fills the slit
    # No slit loss

    # Apply total efficiency
    detected_sky_flux = (
       sky_per_bin
       * eff
    )

    print("Detected sky flux:", detected_sky_flux[0])

    # Apply telescope collecting area
    # photons / s / bin / pixel
    detected_sky_rate = (
        detected_sky_flux
        * area
    )
    print("Detected sky rate:", detected_sky_rate[0])

    # Total sky photons detected
    # photons / bin / pixel
    detected_sky = (
        detected_sky_rate
        * exptime
        * u.s
    )

    detected_sky = detected_sky.to_value()

    source_cnts = source_spec.to_value()
    sky_cnts = detected_sky * npix
    print("Sky cnts", sky_cnts[0])
    read_noise = rdnoise

    #print("Sky Counts")
    #print(sky_cnts)

    #print(source_cnts)
    #print(sky_cnts)
    #print(npix)
    #print(rdnoise)

    variance = (
        source_cnts 
        + sky_cnts 
        + npix * read_noise**2
    )

    noise_spectrum = np.sqrt(variance)
    print("Noise spec:", noise_spectrum[0])

    return noise_spectrum, sky_cnts

# ============================================================
# Make 4 panel plots
# ============================================================

def plot_4spectrum(
    s2n_spec,
    source_spec,
    sky_spec,
    noise_spec,
    wave_eff,
    grism,
    obsmode,
    spec_type,
    vmag,
    seeing,
    exptime,
    nexp,
    slit,
    moon_phase,
    src1,
    src2,
    exec_dir,
):

    
    #print("Moon", moon_phase, type(moon_phase))

    # ------------------------------------------------------------
    # Create 2 x 2 figure
    # ------------------------------------------------------------

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        2,
        2,
        figsize=(12, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 1]},
    )

    # ------------------------------------------------------------
    # Top-left panel: S/N
    # ------------------------------------------------------------

    ax1.plot(
        wave_eff.to_value(u.Angstrom),
        s2n_spec,
        linewidth=1.0,
    )

    ax1.set_ylabel("S/N per binned spectral pixel")
    ax1.grid(alpha=0.3)

    # ------------------------------------------------------------
    # Top-right panel: S/N
    # ------------------------------------------------------------

    ax2.plot(
        wave_eff.to_value(u.Angstrom),
        np.sqrt(nexp) * s2n_spec,
        linewidth=1.0,
    )

    ax2.set_ylabel("S/N per binned spectral pixel")
    ax2.grid(alpha=0.3)

    # ------------------------------------------------------------
    # Bottom-left panel: Source, Sky, Noise
    # ------------------------------------------------------------

    ax3.plot(
        wave_eff.to_value(u.Angstrom),
        source_spec.to_value(source_spec.unit),
        linewidth=1.0,
        label="Source",
    )

    ax3.plot(
        wave_eff.to_value(u.Angstrom),
        sky_spec,
        linewidth=1.0,
        label="Sky",
    )

    ax3.plot(
        wave_eff.to_value(u.Angstrom),
        noise_spec,
        linewidth=1.0,
        label="Noise",
    )

    ax3.set_yscale("log")
    ax3.set_xlabel("Wavelength (Å)")
    ax3.set_ylabel("e- per binned spectral pixel")
    ax3.grid(alpha=0.3)
    ax3.legend()

    # ------------------------------------------------------------
    # Bottom-right panel: Source, Sky, Noise
    # ------------------------------------------------------------

    ax4.plot(
        wave_eff.to_value(u.Angstrom),
        nexp * source_spec.to_value(source_spec.unit),
        linewidth=1.0,
        label="Source",
    )

    ax4.plot(
        wave_eff.to_value(u.Angstrom),
        nexp * sky_spec,
        linewidth=1.0,
        label="Sky",
    )

    ax4.plot(
        wave_eff.to_value(u.Angstrom),
        np.sqrt(nexp) * noise_spec,
        linewidth=1.0,
        label="Noise",
    )

    ax4.set_yscale("log")
    ax4.set_xlabel("Wavelength (Å)")
    ax4.set_ylabel("e- per binned spectral pixel")
    ax4.grid(alpha=0.3)
    ax4.legend()

    # ------------------------------------------------------------
    # Title and observation information
    # ------------------------------------------------------------

    obsmode = obsmode.upper()
    sptype = spec_type.upper()

    # ------------------------------------------------------------
    # Column titles
    # ------------------------------------------------------------

    fig.text(
        0.27, 0.96,
        src1,
        ha="center",
        va="top",
        fontsize=14,
    )

    fig.text(
        0.73, 0.96,
        src2,
        ha="center",
        va="top",
        fontsize=14,
    )

    info = (
        f"ObsMode: {grism}{obsmode}\n"
        f"SpectralType: {sptype}\n"
        f"VMag: {vmag}\n"
        f"ExpTime: {exptime}s\n"
        f"SlitWidth: {slit}\"\n"
        f"Seeing: {seeing}\"\n"
        f"Lunar Day: {moon_phase}"
    )

    ax1.text(
        0.02,
        0.96,
        info,
        transform=ax1.transAxes,
        fontsize=10,
        verticalalignment="top",
        horizontalalignment="left",
        bbox=dict(
            boxstyle="round",
            facecolor="white",
            alpha=0.8,
            edgecolor="black",
        ),
    )

    # ------------------------------------------------------------
    # Save plot
    # ------------------------------------------------------------

    obsmode = obsmode.lower()

    plotname = (
        f"plot_{spec_type}_{vmag}mag_"
        f"{grism}{obsmode}_lp{moon_phase}.png"
    )

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    plt.savefig(
        exec_dir + "/" + plotname,
        dpi=300,
        bbox_inches="tight",
    )

    # plt.show()
    plt.close(fig)

# ============================================================
# Main steering program
# ============================================================

def steer(argv):
    
    # Setup defaults
    grism = '400'
    obsmode = 'm1'
    bin_spatial = 2
    bin_dispersion = 2
    slit = 1.0
    readout = 'slow'
    vmag = 15.0
    spec_type = 'A0V'
    metallicity = 'solar'
    moon_phase = '7'
    seeing = 1.0
    exptime = 60.0
    nexp = 1

    i = 1
    while i < len(argv):
        if argv[i] == '-h':
            print(__doc__)
            return
        elif argv[i] == '-grating':
            i += 1
            grism = argv[i]
        elif argv[i] == '-mode':
            i += 1
            obsmode = argv[i]
        elif argv[i] == '-bin_sky':
            i += 1
            bin_spatial = int(argv[i])
        elif argv[i] == '-bin_disp':
            i += 1
            bin_dispersion = int(argv[i])
        elif argv[i] == '-slit':
            i += 1
            slit = float(argv[i])
        elif argv[i] == '-gain':
            i += 1
            readout = argv[i]
        elif argv[i] == '-vmag':
            i += 1
            vmag = float(argv[i])
        elif argv[i] == '-sptype':
            i += 1
            spec_type = argv[i]
        elif argv[i] == '-metal':
            i += 1
            metallicity = argv[i]
        elif argv[i] == '-lp':
            i += 1
            moon_phase = argv[i]
        elif argv[i] == '-see':
            i += 1
            seeing = float(argv[i])
        elif argv[i] == '-expt':
            i += 1
            exptime = float(argv[i])
        elif argv[i] =='-num':
            i += 1
            nexp = int(argv[i])
        i += 1    

    # Get gain and readout noise
    if readout == "slow" or readout == "fast":
        #print(f"Readout \"{readout}\"")
        gain, rdnoise = get_gain(readout)
    else: 
        print(f"Readout must be \"slow\" or \"fast\". Readout of {readout} is not supported.  Exiting.")
        sys.exit()
    print(f"Calculations will be made using a gain of {gain} e/ADU and a readout noise of {rdnoise}e.")


    # Get spectral template
    if metallicity != 'solar' and metallicity != 'weak' and metallicity != 'rich':
        print(f"Metallicity must be \"solar\", \"weak\", or \"rich\".  Metallicity of \"{metallicity}\" is not supported.  Exiting.")
        sys.exit()

    spec_file = get_spec_template(spec_type, metallicity)
    cap_sptype = spec_type.upper()
    
    if spec_file == 'none.txt':
        print(f"No template of spectral type {cap_sptype} and {metallicity} metallicity found.  Exiting")
        sys.exit()
    print(f"Using spectral template {spec_file} for a {cap_sptype} star with {metallicity} metallicity")
    spec_file = SPECDIR + spec_file

    # Get throughput curve for grating and observing mode
    eff_file = get_efficiency(grism, obsmode)
    cap_obsmode = obsmode.upper()
    temp_eff_file = eff_file
    eff_file = EFFDIR + eff_file
    if not os.path.isfile(eff_file):
       print(f"No efficiency curve for {grism}{cap_obsmode} is found. Exiting.")
       sys.exit()

    print(f"Using efficiency curve {temp_eff_file}.")

    # Get lunar phase
    lp = int(moon_phase)
    #print(lp)
    if lp < 0 or lp > 15:
       print(f"Lunar phase is days past full moon and must be between 0 and 15. Lunar phase of {lp} is out of range. Exiting.")
       sys.exit()
    
    moon_phase, sky_model = get_sky_model(moon_phase)
    
    if not os.path.isfile(sky_model):
        print(f"Sky model {sky_model} does not exist.  Exiting.")
        sys.exit
    temp_sky = os.path.basename(sky_model)
    
    print(f"Using sky model {temp_sky}.")

    # Get slit loss
    trans, loss = slit_loss(slit, seeing)
    print(f"For a slit width of {slit}\" and seeing of {seeing}\", the calculated slit loss is {loss:.2f}.")
    print(f"Slit transmission: {trans:.2f}")

    # Get dispersion per pixel
    disp = get_dispersion(grism, bin_dispersion)
    if disp < 0.1:
        print(f"No dispersion information for {grism}l/mm grating.  Exiting.")
        sys.exit()
    print(f"Using the {grism}l/mm grating and spectral binning of {bin_dispersion}.")
    print(f"The dispersion is {disp} angstrom/pixel.")   

    wave_eff, source_spec = make_source_spectrum(spec_file, vmag, eff_file, trans, disp, exptime, sa)


    noise_spec, sky_spec = make_noise_spectrum(source_spec, sky_model, eff_file, rdnoise, disp, exptime, sa, pix_scale, bin_dispersion, bin_spatial, seeing)

    s2n = source_spec / noise_spec

    src1 = f"Single {exptime}s Exposure {grism}{obsmode}"
    src2 = rf"{nexp} $\times$ {exptime}s Exposures {grism}{obsmode}"

    plot_4spectrum(s2n, source_spec, sky_spec, noise_spec,
        wave_eff, grism, obsmode, spec_type, vmag, seeing, exptime,
        nexp, slit, moon_phase, src1, src2, exec_dir)



if __name__ == "__main__":
     import sys
     if len(sys.argv) > 1:
         steer(sys.argv)
     else:
         print(__doc__)
