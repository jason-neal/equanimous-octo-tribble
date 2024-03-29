#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function

import argparse
import os
import urllib

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as sci
from astropy.io import fits
from astropy.modeling import fitting, models

from gooey import Gooey, GooeyParser


def _download_spec(fout):
    """Download a spectrum from my personal web page."""
    import requests
    spec = fout.rpartition('/')[-1]
    url = 'http://www.astro.up.pt/~dandreasen/{0!s}'.format(spec)
    spectrum = requests.get(url)
    with open(fout, 'w') as file:
        file.write(spectrum.content)


class Cursor:
    """Get a crosshair at the cursor's position.

    The code is from here:
    http://matplotlib.org/examples/pylab_examples/cursor_demo.html
    """

    def __init__(self, ax):
        self._ax = ax
        self.lx = ax.axhline(color='b', lw=2, alpha=0.7)
        self.ly = ax.axvline(color='b', lw=2, alpha=0.7)

    def mouse_move(self, event):
        if not event.inaxes:
            return
        x, y = event.xdata, event.ydata
        self.lx.set_ydata(y)
        self.ly.set_xdata(x)
        plt.draw()


def ccf_astro(spectrum1, spectrum2, rvmin=0, rvmax=200, drv=1):
    """Make a CCF between 2 spectra and find the RV.

    :spectrum1: The stellar spectrum
    :spectrum2: The model, sun or telluric
    :dv: The velocity step
    :returns: The RV shift
    """
    # Calculate the cross correlation
    s = False
    w, f = spectrum1
    tw, tf = spectrum2
    if not (len(w) and len(tw)):
        return 0, 0, 0, 0, 0
    c = 299792.458
    drvs = np.arange(rvmin, rvmax, drv)
    cc = np.zeros(len(drvs))
    for i, rv in enumerate(drvs):
        fi = sci.interp1d(tw * (1.0 + rv / c), tf)
        # Shifted template evaluated at location of spectrum
        try:
            fiw = fi(w)
            cc[i] = np.sum(f * fiw)
        except ValueError:
            s = True
            fiw = 0
            cc[i] = 0

    if s:
        print('Warning: Lower the bounds on RV')

    if not np.any(cc):
        return 0, 0, 0, 0, 0

    # Fit the CCF with a gaussian
    cc[cc == 0] = np.mean(cc)
    cc = (cc - min(cc)) / (max(cc) - min(cc))
    RV, g = _fit_ccf(drvs, cc)
    return int(RV), drvs, cc, drvs, g(drvs)


def _fit_ccf(rv, ccf):
    """Fit the CCF with a 1D gaussian.

    :rv: The RV vector
    :ccf: The CCF values
    :returns: The RV, and best fit gaussian

    """
    ampl = 1
    mean = rv[ccf == ampl]
    i_peak = np.where(ccf == ampl)[0]

    g_init = models.Gaussian1D(amplitude=ampl, mean=mean, stddev=5)
    fit_g = fitting.LevMarLSQFitter()

    try:
        g = fit_g(g_init, rv[i_peak - 10:i_peak + 10],
                  ccf[i_peak - 10:i_peak + 10])
    except TypeError:
        print('Warning: Not able to fit a gaussian to the CCF')
        return 0, g_init
    RV = g.mean.value
    return RV, g


def nrefrac(wavelength, density=1.0):
    """Calculate refractive index of air from Cauchy formula.

    Input:
    wavelength in Angstrom, density of air in amagat (relative to STP,
    e.g. ~10% decrease per 1000m above sea level). Returns N = (n-1) *
    1.e6.

    The IAU standard for conversion from air to vacuum wavelengths is given
    in Morton (1991, ApJS, 77, 119). For vacuum wavelengths (VAC) in
    Angstroms, convert to air wavelength (AIR) via:

    AIR = VAC / (1.0 + 2.735182E-4 + 131.4182 / VAC^2 + 2.76249E8 / VAC^4)
    """
    wl = np.array(wavelength)

    wl2inv = (1.e4 / wl) ** 2
    refracstp = 272.643 + 1.2288 * wl2inv + 3.555e-2 * wl2inv ** 2
    return density * refracstp


def dopplerShift(wvl, flux, v, edgeHandling='firstlast', fill_value=None):
    """Doppler shift a given spectrum.

    This code is taken from the PyAstronomy project:
    https://github.com/sczesla/PyAstronomy
    All credit to the author.

    A simple algorithm to apply a Doppler shift
    to a spectrum. This function, first, calculates
    the shifted wavelength axis and, second, obtains
    the new, shifted flux array at the old, unshifted
    wavelength points by linearly interpolating.

    Due to the shift, some bins at the edge of the
    spectrum cannot be interpolated, because they
    are outside the given input range. The default
    behavior of this function is to return numpy.NAN
    values at those points. One can, however, specify
    the `edgeHandling` parameter to choose a different
    handling of these points.

    If "firstlast" is specified for `edgeHandling`,
    the out-of-range points at the red or blue edge
    of the spectrum will be filled using the first
    (at the blue edge) and last (at the red edge) valid
    point in the shifted, i.e., the interpolated, spectrum.

    If "fillValue" is chosen for edge handling,
    the points under consideration will be filled with
    the value given through the `fillValue` keyword.

    .. warning:: Shifting a spectrum using linear
                interpolation has an effect on the
                noise of the spectrum. No treatment
                of such effects is implemented in this
                function.

    Parameters
    ----------
    wvl : array
        Input wavelengths in A.
    flux : array
        Input flux.
    v : float
        Doppler shift in km/s
    edgeHandling : string, {"fillValue", "firstlast"}, optional
        The method used to handle the edges of the
        output spectrum.
    fillValue : float, optional
        If the "fillValue" is specified as edge handling method,
        the value used to fill the edges of the output spectrum.

    Returns
    -------
    nflux : array
        The shifted flux array at the *old* input locations.
    wlprime : array
        The shifted wavelength axis.
    """
    # Shifted wavelength axis
    wlprime = wvl * (1.0 + v / 299792.458)

    f = sci.interp1d(wlprime, flux, bounds_error=False, fill_value=np.nan)
    nflux = f(wlprime)

    if edgeHandling == "firstlast":
        firsts = []
        # Search for first non-NaN value save indices of
        # leading NaN values
        for i, nfluxi in enumerate(nflux):
            if np.isnan(nfluxi):
                firsts.append(i)
            else:
                firstval = nfluxi
                break

        # Do the same for trailing NaNs
        lasts = []
        for i, nfluxi in enumerate(nflux[::-1]):
            if np.isnan(nfluxi):
                lasts.append(i)
            else:
                lastval = nfluxi
                break

        # Use first and last non-NaN value to
        # fill the nflux array
        if fill_value:
            nflux[firsts] = fill_value
            nflux[lasts] = fill_value
        else:
            nflux[firsts] = firstval
            nflux[lasts] = lastval
    return nflux, wlprime


def get_wavelength(hdr):
    """Return the wavelength vector calculated from the header of a FITS file.

    :hdr: Header from a FITS ('CRVAL1', 'CDELT1', and 'NAXIS1' is required as
            keywords)
    :returns: Equidistant wavelength vector

    """
    w0, dw, n = hdr['CRVAL1'], hdr['CDELT1'], hdr['NAXIS1']
    w1 = w0 + dw * n
    return np.linspace(w0, w1, n, endpoint=False)


@Gooey(program_name='Plot fits - Easy 1D fits plotting', default_size=(610, 730))
def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = GooeyParser(
        description='Plot 1D fits files with wavelength information in the header.')
    parser.add_argument('fname',
                        action='store',
                        widget='FileChooser',
                        help='Input fits file')
    parser.add_argument('fname2',
                        action='store',
                        widget='FileChooser',
                        help='Input fits file')
    parser.add_argument('-m', '--model',
                        default=False,
                        widget='FileChooser',
                        help='If not the Sun shoul be used as a model, put'
                        ' the model here (only support BT-Settl for the'
                        ' moment)')
    parser.add_argument('-s', '--sun',
                        help='Over plot solar spectrum',
                        action='store_true')
    parser.add_argument('-t', '--telluric',
                        help='Over plot telluric spectrum',
                        action='store_true')
    parser.add_argument('-r', '--rv',
                        help='RV shift to observed spectra in km/s',
                        default=False,
                        type=float)
    parser.add_argument('-r1', '--rv1',
                        help='RV shift to model/solar spectrum in km/s',
                        default=False,
                        type=float)
    parser.add_argument('-r2', '--rv2',
                        help='RV shift to telluric spectra in km/s',
                        default=False,
                        type=float)
    parser.add_argument('-l', '--lines',
                        help='Lines to plot on top (multiple lines is an'
                        ' option). If multiple lines needs to be plotted, then'
                        ' separate with a space',
                        default=False,
                        nargs='+',
                        type=float)
    parser.add_argument('-c', '--ccf',
                        default='none',
                        choices=['none', 'sun', 'model', 'telluric', 'both'],
                        help='Calculate the CCF for Sun/model or tellurics '
                        'or both.')
    parser.add_argument('--ftype', help='Select which type the fits file is',
                        choices=['ARES', 'CRIRES', 'DRACS'], default='CRIRES')
    parser.add_argument('--ftype2', help='Select which type the fits file is',
                        choices=['ARES', 'CRIRES', 'DRACS'], default='DRACS')
    parser.add_argument('--fitsext', help='Select fits extention, Default 0.',
                        choices=['0', '1', '2', '3', '4'], default='0')
    parser.add_argument('--fitsext2', help='Select fits extention, Default 0.',
                        choices=['0', '1', '2', '3', '4'], default='0')
    # parser.add_argument('--fitsext', default=0, type=int,
    #                    help='Select fits extention, 0 = Primary header')
    args = parser.parse_args()
    return args


def main(fname, fname2, lines=False, model=False, telluric=False, sun=False,
         rv=False, rv1=False, rv2=False, ccf='none', ftype='CRIRES', ftype2='DRACS',
         fitsext='0', fitsext2='0'):
    """Plot a fits file with extensive options.

    :fname: Input spectra
    :lines: Absorption lines
    :model: Model spectrum
    :telluric: Telluric spectrum
    :sun: Solar spectrum
    :rv: RV of input spectrum
    :rv1: RV of Solar/model spectrum
    :rv2: RV of telluric spectrum
    :ccf: Calculate CCF (sun, model, telluric, both)
    :ftype: Type of fits file (ARES or CRIRES)
    :fitsext: Slecet fits extention to use (0, 1, 2, 3, 4)
    :returns: RV if CCF have been calculated
    """
    print('\n-----------------------------------')
    path = os.path.expanduser('~/.plotfits/')
    pathsun = os.path.join(path, 'solarspectrum_01.fits')
    pathtel = os.path.join(path, 'telluric_NIR.fits')
    pathwave = os.path.join(path, 'WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
    if os.path.isdir(path):
        if sun and (not os.path.isfile(pathsun)):
            print('Downloading solar spectrum...')
            _download_spec(pathsun)
        if telluric and (not os.path.isfile(pathtel)):
            print('Downloading telluric spectrum...')
            _download_spec(pathtel)
        if model and (not os.path.isfile(pathwave)):
            print('Downloading wavelength vector for model...')
            url = 'ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS//WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'
            urllib.urlretrieve(url, pathwave)
    else:
        os.mkdir(path)
        print('{0!s} Created'.format(path))
        print('Downloading solar spectrum...')
        _download_spec(pathsun)
        print('Downloading telluric spectrum...')
        _download_spec(pathtel)

    fitsext = int(fitsext)
    fitsext2 = int(fitsext2)
    if ftype == 'ARES':
        i_flux = fits.getdata(fname, fitsext)
        hdr = fits.getheader(fname, fitsext)
        w = get_wavelength(hdr)
        label1 = "Star - ARES"
    elif ftype == 'CRIRES':
        d = fits.getdata(fname, fitsext)
        hdr = fits.getheader(fname, fitsext)
        i_flux = d['Extracted_OPT']
        w = d['Wavelength'] * 10
        label1 = "Star - CRIRES"
    elif ftype == 'DRACS':
        d = fits.getdata(fname)
        hdr = fits.getheader(fname)
        try:
            i_flux = d['Extracted_DRACS']
            w = d['Wavelength'] * 10
            label1 = "Star - CALIBRATED DRACS "
        except:
            try:
                i_flux = d["Combined"]
            except:
                i_flux = d
            w = np.linspace(hdr["HIERARCH ESO INS WLEN STRT"],
                            hdr["HIERARCH ESO INS WLEN END"], len(i_flux)) * 10
            label1 = "Star - DRACS (uncalib)"
    if ftype2 == 'ARES':
        i_flux2 = fits.getdata(fname2, fitsext2)
        hdr = fits.getheader(fname2, fitsext2)
        w = get_wavelength(hdr)
        label2 = "Star - ARES"
    elif ftype2 == 'CRIRES':
        d2 = fits.getdata(fname2, fitsext2)
        hdr2 = fits.getheader(fname2, fitsext2)
        i_flux2 = d2['Extracted_OPT']
        w2 = d2['Wavelength'] * 10
        label2 = "Star - CRIRES"
    elif ftype2 == 'DRACS':
        d2 = fits.getdata(fname2)
        hdr2 = fits.getheader(fname2)
        try:
            i_flux2 = d2['Extracted_DRACS']
            w2 = d2['Wavelength'] * 10
            label2 = "Star - CALIBRATED DRACS "
        except:
            try:
                i_flux2 = d2["Combined"]
            except:
                i_flux2 = d2
            w2 = np.linspace(hdr2['ESO INS WLEN MIN'],
                             hdr2['ESO INS WLEN MAX'], len(i_flux2)) * 10
            label2 = "Star - DRACS (uncalib)"
        # w = np.linspace(hdr['ESO INS WLEN MIN'], hdr['ESO INS WLEN MAX'], len(I)) * 10
    i_flux /= np.median(i_flux)
    i_flux2 /= np.median(i_flux2)
    # Normalization (use first 50 points below 1.2 as constant continuum)
    maxes = i_flux[(i_flux < 1.2)].argsort()[-50:][::-1]
    i_flux /= np.median(i_flux[maxes])
    maxes2 = i_flux2[(i_flux2 < 1.2)].argsort()[-50:][::-1]
    i_flux2 /= np.median(i_flux2[maxes2])
    # hdr = fits.getheader(fname)
    dw = 10  # Some extra coverage for RV shifts

    if rv:
        i_flux, w = dopplerShift(wvl=w, flux=i_flux, v=rv, fill_value=0.95)
        i_flux2, w2 = dopplerShift(wvl=w2, flux=i_flux2, v=rv, fill_value=0.95)
    w0, w1 = w[0] - dw, w[-1] + dw
    w02, w12 = w2[0] - dw, w2[-1] + dw

    if sun and not model:
        i_sun = fits.getdata(pathsun)
        hdr = fits.getheader(pathsun)
        w_sun = get_wavelength(hdr)
        i = (w_sun > w0) & (w_sun < w1)
        w_sun = w_sun[i]
        i_sun = i_sun[i]
        if len(w_sun) > 0:
            i_sun /= np.median(i_sun)
            if ccf in ['sun', 'both'] and rv1:
                print('Warning: RV set for Sun. Calculate RV with CCF')
            if rv1 and ccf not in ['sun', 'both']:
                i_sun, w_sun = dopplerShift(
                    wvl=w_sun, flux=i_sun, v=rv1, fill_value=0.95)
        else:
            print('Warning: Solar spectrum not available in wavelength range.')
            sun = False
    elif sun:
        print('Warning: Both solar spectrum and a model spectrum are selected. Using model spectrum.')
        sun = False

    if model:
        i_flux_mod = fits.getdata(model)
        hdr = fits.getheader(model)
        if 'WAVE' in hdr.keys():
            w_mod = fits.getdata(pathwave)
        else:
            w_mod = get_wavelength(hdr)
        nre = nrefrac(w_mod)  # Correction for vacuum to air (ground based)
        w_mod = w_mod / (1 + 1e-6 * nre)
        i = (w_mod > w0) & (w_mod < w1)
        w_mod = w_mod[i]
        i_flux_mod = i_flux_mod[i]
        if len(w_mod) > 0:
            # https://phoenix.ens-lyon.fr/Grids/FORMAT
            # i_flux_mod = 10 ** (i_flux_mod-8.0)
            i_flux_mod /= np.median(i_flux_mod)
            # Normalization (use first 50 points below 1.2 as continuum)
            maxes = i_flux_mod[(i_flux_mod < 1.2)].argsort()[-50:][::-1]
            i_flux_mod /= np.median(i_flux_mod[maxes])
            if ccf in ['model', 'both'] and rv1:
                print('Warning: RV set for model. Calculate RV with CCF')
            if rv1 and ccf not in ['model', 'both']:
                i_flux_mod, w_mod = dopplerShift(
                    wvl=w_mod, flux=i_flux_mod, v=rv1, fill_value=0.95)
        else:
            print('Warning: Model spectrum not available in wavelength range.')
            model = False

    if telluric:
        i_tel = fits.getdata(pathtel)
        hdr = fits.getheader(pathtel)
        w_tel = get_wavelength(hdr)
        i = (w_tel > w0) & (w_tel < w1)
        w_tel = w_tel[i]
        i_tel = i_tel[i]
        if len(w_tel) > 0:
            i_tel /= np.median(i_tel)
            if ccf in ['telluric', 'both'] and rv2:
                print('Warning: RV set for telluric, Calculate RV with CCF')
            if rv2 and ccf not in ['telluric', 'both']:
                i_tel, w_tel = dopplerShift(
                    wvl=w_tel, flux=i_tel, v=rv2, fill_value=0.95)
        else:
            print('Warning: Telluric spectrum not available in wavelength range.')
            telluric = False

    rvs = {}
    if ccf != 'none':
        if ccf in ['sun', 'both'] and sun:
            # remove tellurics from the Solar spectrum
            if telluric:
                print('Correcting solar spectrum for tellurics...')
                i_sun = i_sun / i_tel
            print('Calculating CCF for the Sun...')
            rv1, r_sun, c_sun, x_sun, y_sun = ccf_astro(
                (w, -i_flux + 1), (w_sun, -i_sun + 1))
            if rv1 != 0:
                print('Shifting solar spectrum...')
                i_sun, w_sun = dopplerShift(
                    w_sun, i_sun, v=rv1, fill_value=0.95)
                rvs['sun'] = rv1
                print('DONE')

        if ccf in ['model', 'both'] and model:
            print('Calculating CCF for the model...')
            rv1, r_mod, c_mod, x_mod, y_mod = ccf_astro(
                (w, -i_flux + 1), (w_mod, -i_flux_mod + 1))
            if rv1 != 0:
                print('Shifting model spectrum...')
                i_flux_mod, w_mod = dopplerShift(
                    w_mod, i_flux_mod, v=rv1, fill_value=0.95)
                rvs['model'] = rv1
                print('DONE')

        if ccf in ['telluric', 'both'] and telluric:
            print('Calculating CCF for the model...')
            rv2, r_tel, c_tel, x_tel, y_tel = ccf_astro(
                (w, -i_flux + 1), (w_tel, -i_tel + 1))
            if rv2 != 0:
                print('Shifting telluric spectrum...')
                i_tel, w_tel = dopplerShift(
                    w_tel, i_tel, v=rv2, fill_value=0.95)
                rvs['telluric'] = rv2
                print('DONE')

    if len(rvs) == 0:
        ccf = 'none'

    if ccf != 'none':
        from matplotlib.gridspec import GridSpec
        fig = plt.figure(figsize=(16, 5))
        gs = GridSpec(1, 5)
        if len(rvs) == 1:
            gs.update(wspace=0.25, hspace=0.35, left=0.05, right=0.99)
            ax1 = plt.subplot(gs[:, 0:-1])
            ax2 = plt.subplot(gs[:, -1])
            ax2.set_yticklabels([])
        elif len(rvs) == 2:
            gs.update(wspace=0.25, hspace=0.35, left=0.01, right=0.99)
            ax1 = plt.subplot(gs[:, 1:4])
            ax2 = plt.subplot(gs[:, 0])
            ax3 = plt.subplot(gs[:, -1])
            ax2.set_yticklabels([])
            ax3.set_yticklabels([])
    else:
        fig = plt.figure(figsize=(16, 5))
        ax1 = fig.add_subplot(111)

    # Start in pan mode with these two lines
    manager = plt.get_current_fig_manager()
    manager.toolbar.pan()

    # Use nice numbers on x axis (y axis is normalized)...
    x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax1.xaxis.set_major_formatter(x_formatter)

    if sun and not model:
        ax1.plot(w_sun, i_sun, '-g', lw=2, alpha=0.6, label='Sun')
    if telluric:
        ax1.plot(w_tel, i_tel, '-r', lw=2, alpha=0.5, label='Telluric')
    if model:
        ax1.plot(w_mod, i_flux_mod, '-g', lw=2, alpha=0.5, label='Model')
    ax1.plot(w, i_flux, '-k', lw=2, label=label1)
    ax1.plot(w2, i_flux2, '-m', lw=2, label=label2)

    # Add crosshair
    xlim = ax1.get_xlim()
    cursor = Cursor(ax1)
    plt.connect('motion_notify_event', cursor.mouse_move)
    ax1.set_xlim(xlim)

    if lines:
        y0, y1 = ax1.get_ylim()
        if rv1:
            shift = (1.0 + rv1 / 299792.458)
            for line in lines:
                ax1.vlines(line * shift, y0, y1,
                           linewidth=2, color='m', alpha=0.5)
                ax1.text(line * shift - 0.7, 1.2, str(line), rotation=90)
        else:
            for line in lines:
                ax1.vlines(line, y0, y1, linewidth=2, color='m', alpha=0.5)
                ax1.text(line - 0.7, 1.2, str(line), rotation=90)

    ax1.set_xlabel('Wavelength')
    ax1.set_ylabel('"Normalized" flux')

    if len(rvs) == 1:
        if 'sun' in rvs.keys():
            ax2.plot(r_sun, c_sun, '-k', lw=2)
            ax2.plot(x_sun, y_sun, '--r', lw=2)
            ax2.set_title('CCF (sun)')
        if 'model' in rvs.keys():
            ax2.plot(r_mod, c_mod, '-k', lw=2)
            ax2.plot(x_mod, y_mod, '--r', lw=2)
            ax2.set_title('CCF (mod)')
        if 'telluric' in rvs.keys():
            ax2.plot(r_tel, c_tel, '-k', lw=2)
            ax2.plot(x_tel, y_tel, '--r', lw=2)
            ax2.set_title('CCF (tel)')
        ax2.set_xlabel('RV [km/s]')

    elif len(rvs) == 2:
        if 'sun' in rvs.keys():
            ax2.plot(r_sun, c_sun, '-k', lw=2)
            ax2.plot(x_sun, y_sun, '--r', lw=2)
            ax2.set_title('CCF (sun)')
        if 'model' in rvs.keys():
            ax2.plot(r_mod, c_mod, '-k', lw=2)
            ax2.plot(x_mod, y_mod, '--r', lw=2)
            ax2.set_title('CCF (mod)')
        ax3.plot(r_tel, c_tel, '-k', lw=2)
        ax3.plot(x_tel, y_tel, '--r', lw=2)
        ax3.set_title('CCF (tel)')

        ax2.set_xlabel('RV [km/s]')
        ax3.set_xlabel('RV [km/s]')

    if rv:
        ax1.set_title('{0!s}\nRV correction: {1!s} km/s'.format(fname, rv))
    elif rv1 and rv2:
        ax1.set_title(
            '{0!s}\nSun/model: {1!s} km/s, telluric: {2!s} km/s'.format(fname, rv1, rv2))
    elif rv1:
        ax1.set_title('{0!s}\nSun/model: {1!s} km/s'.format(fname, rv1))
    elif rv2:
        ax1.set_title('{0!s}\nTelluric: {1!s} km/s'.format(fname, rv2))
    elif ccf == 'model':
        ax1.set_title('{0!s}\nModel(CCF): {1!s} km/s'.format(fname, rv1))
    elif ccf == 'sun':
        ax1.set_title('{0!s}\nSun(CCF): {1!s} km/s'.format(fname, rv1))
    elif ccf == 'telluric':
        ax1.set_title('{0!s}\nTelluric(CCF): {1!s} km/s'.format(fname, rv2))
    elif ccf == 'both':
        ax1.set_title(
            '{0!s}\nSun/model(CCF): {1!s} km/s, telluric(CCF): {2!s} km/s'.format(fname, rv1, rv2))
    else:
        ax1.set_title(fname)
    if sun or telluric or model:
        ax1.legend(loc=3, frameon=False)
    plt.show()

    return rvs


if __name__ == '__main__':
    args = vars(_parser())
    fname = args.pop('fname')
    opts = {k: args[k] for k in args}

    main(fname, **opts)
