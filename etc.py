#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
TBA.

Note: The HAWKI instrument has no K filter, only Ks, so I assumed that the
      UX text was referring to the Ks filter and not K.

Created on Fri Nov 11 14:34:00 2022

@author: teutoburg
"""

__version__ = "0.2"

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



class HawkiEtc():
    default_kwargs = {"filter_name": "Ks",
                      "instrument": "HAWKI",
                      "observatory": "Paranal"}
    aperture = aperture = 8.2 * u.m  # VLT aperture

    def __init__(self):
        pass

    @property
    def aperture_area(self):
        """Aperture area of the telescope. Read-only property."""
        return 0.785398 * (self.aperture**2)  # pi/4 -> circle area

    @property
    def pixel_scale(self):
        """Pixel scale of the telescope. Read-only property."""
        return 4*2048*2048*u.pixel / self.aperture_area

    @staticmethod
    def _get_noise(target_counts, sky_counts_pp, dit):
        # HACK: The following factors are directly taken from ESO HAWK-I ETC
        drs_factor = 1.1
        read_out_noise = 5 * u.electron/u.pixel  # per DIT
        dark_current = .01 * u.electron/u.pixel/u.s

        # HACK: The pixel number below is most likely calculated from the PSF
        #       by official ESO HAWK-I ETC using a Image Quality FWHM of
        #       0.6 arcsec, or a seeing of 0.8 arcsec FWHM, which is described
        #       as the 50% value for observing conditions. It has been fixed
        #       for the scope of this project.
        n_pixel = 101 * u.pixel

        dark = dark_current * dit
        readout = read_out_noise.value**2 * u.electron/u.pixel
        background = drs_factor * n_pixel * (sky_counts_pp + dark + readout)

        return np.sqrt(target_counts + background)

    def signal_to_noise(self, target_counts, sky_counts_pp, dit, n_dit):
        """
        Calculate S/N ratio for given target and sky counts and exposure time.

        Formula adapted from the official ESO HAWK-I ETC.

        Parameters
        ----------
        target_counts : astropy.Quantity
            DESCRIPTION.
        sky_counts_pp : astropy.Quantity
            DESCRIPTION.
        dit : astropy.Quantity
            Detector integration time in seconds.
        n_dit : int
            Number of integrations.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        signal = np.sqrt(n_dit) * target_counts
        noise = self._get_noise(target_counts, sky_counts_pp, dit)
        return (signal / noise).value

    @staticmethod
    def _to_electrons(flux, dit, aperture_area):
        # HACK: The following conversion factor has been derived to match the
        #       values to those produced by the official ESO HAWK-I ETC, it
        #       has no physical derivation! When used with AB magnitudes,
        #       the factor is closer to 0.4...
        eff = .3581 * u.electron / u.photon
        return flux * dit * aperture_area * eff

    def create_sky(self, pixel_scale, dit, **kwargs):
        """
        Estimate number of electrons per pixel expected from sky background.

        Parameters
        ----------
        pixel_scale : astropy.Quantity
            Pixel scale of the optical system in pixel per square meter.
        dit : astropy.Quantity
            Detector integration time in seconds.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        sky_count : astropy.Quantity
            Sky level in electrons per pixel.

        """
        # TODO: investigate if this could be sped up by a decorator storing the
        #       output of last few calls for given parameters...
        defaults = {"filter_name": "Ks",
                    "instrument": "HAWKI",
                    "observatory": "Paranal",
                    "airmass": 1.5, "pwv": 10.0}
        therm_dict = {"incl_therm": "Y", "therm_t1": 285.0, "therm_e1": 0.20,
                      "therm_t2": 288.0, "therm_e2": 0.10,  "therm_t3": 33.0,
                      "therm_e3": 0.01}
        parameters = defaults | therm_dict | kwargs
        sky = hmbp.in_skycalc_background(**parameters)
        return sky * dit / pixel_scale * u.electron/u.photon


# dit = 60 * u.s
# n_dit = 60

# assert dit * n_dit == 1 * u.h  # total observing time should be 1 h

# flux = hmbp.for_flux_in_filter(flux=22*u.mag, **defaults)
# sky = create_sky(pixel_scale)
# target_counts = _to_electrons(flux, dit, aperture_area)
# print(target_counts)
# print(signal_to_noise(target_counts, sky, dit, n_dit))
