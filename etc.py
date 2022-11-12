#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
TBA.

Note: The HAWKI instrument has no K filter, only Ks, so I assumed that the
      UX text was referring to the Ks filter and not K.

Created on Fri Nov 11 14:34:00 2022

@author: teutoburg
"""

__version__ = "0.1"

import numpy as np

from astropy import units as u

import hmbp

__author__ = "Fabian Haberhauer"
__copyright__ = "Copyright 2022"
__credits__ = []
__license__ = "GPL"
__maintainer__ = "Fabian Haberhauer"
__email__ = "fabian.haberhauer@univie.ac.at"
__status__ = "Prototype"


defaults = {"filter_name": "Ks",
            "instrument": "HAWKI",
            "observatory": "Paranal"}


def _noise(target_counts, sky_counts_pp, dit):
    # HACK: The following factors are directly taken from the ESO HAWK-I ETC
    drs_factor = 1.1
    read_out_noise = 5 * u.electron/u.pixel  # per DIT
    dark_current = .01 * u.electron/u.pixel/u.s

    # HACK: The pixel number below is most likely calculated from the PSF by
    #       official ESO HAWK-I ETC using a Image Quality FWHM of 2.43 arcsec.
    #       It has been fixed for the scope of this project.
    n_pixel = 1637 * u.pixel

    dark = dark_current * dit
    readout = read_out_noise.value**2 * u.electron/u.pixel
    background = drs_factor * n_pixel * (sky_counts_pp + dark + readout)

    return np.sqrt(target_counts + background)


def signal_to_noise(target_counts, sky_counts_pp, dit, n_dit):
    """Formula adapted from the official ESO HAWK-I ETC."""
    signal = np.sqrt(n_dit) * target_counts
    return (signal / _noise(target_counts, sky_counts_pp, dit)).value


def _to_electrons(flux, dit, aperture_area):
    # HACK: The following conversion factor has been derived to match the
    #       values to those produced by the official ESO HAWK-I ETC, it has no
    #       physical derivation! When used with AB magnitudes, the factor is
    #       closer to 0.4...
    eff = .3581 * u.electron / u.photon
    return flux * dit * aperture_area * eff


def create_sky(pixel_scale):
    therm_dict = {"incl_therm": "Y", "therm_t1": 285.0, "therm_e1": 0.20,
                  "therm_t2": 288.0, "therm_e2": 0.10,  "therm_t3": 33.0,
                  "therm_e3": 0.01}

    sky = hmbp.in_skycalc_background("Ks", instrument="HAWKI",
                                     observatory="Paranal",
                                     airmass=1.5, pwv=10.0,
                                     **therm_dict)  # * u.arcsec**(-2)
    return sky * dit / pixel_scale * u.electron/u.photon


dit = 60 * u.s
n_dit = 60

assert dit * n_dit == 1 * u.h  # total observing time should be 1 h

# conditions: https://www.eso.org/sci/facilities/paranal/astroclimate/site.html

aperture = 8.2 * u.m  # VLT aperture
aperture_area = np.pi/4 * (aperture**2)
pixel_scale = 4*2048*2048*u.pixel / aperture_area


flux = hmbp.for_flux_in_filter(flux=21*u.mag, **defaults)
sky = create_sky(pixel_scale)
target_counts = _to_electrons(flux, dit, aperture_area)
print(target_counts)
print(signal_to_noise(target_counts, sky, dit, n_dit))

