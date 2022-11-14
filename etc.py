#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
TBA.

Requires at least Python 3.9 for correct type hints with Astropy Quantities!

Note: The HAWKI instrument has no K filter, only Ks, so I assumed that the
      UX text was referring to the Ks filter and not K.

Created on Fri Nov 11 14:34:00 2022

@author: teutoburg
"""

__version__ = "0.2"

import numpy as np
from functools import lru_cache

from scipy.optimize import fsolve
from astropy import units as u
from astropy.units import Quantity

import hmbp

from utils import HiddenPrints

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
    def aperture_area(self) -> Quantity[u.m**2]:
        """Aperture area of the telescope. Read-only property."""
        return 0.785398 * (self.aperture**2)  # pi/4 -> circle area

    @property
    def pixel_scale(self) -> Quantity[u.pixel/u.m**2]:
        """Pixel scale of the telescope. Read-only property."""
        return 4*2048*2048*u.pixel / self.aperture_area

    @staticmethod
    def _get_noise(target_counts: Quantity[u.electron],
                   sky_counts_pp: Quantity[u.electron/u.pixel],
                   dit: Quantity[u.s]) -> Quantity[u.electron]:
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
        noise = target_counts + background
        return noise

    def signal_to_noise(self, target_counts: Quantity[u.electron],
                        sky_counts_pp: Quantity[u.electron/u.pixel],
                        dit: Quantity[u.s], n_dit: int) -> float:
        """
        Calculate S/N ratio for given target and sky counts and exposure time.

        Formula adapted from the official ESO HAWK-I ETC.

        Parameters
        ----------
        target_counts : astropy.Quantity
            Total electron count from source per DIT.
        sky_counts_pp : astropy.Quantity
            Total electron count per pixel from sky background per DIT.
        dit : astropy.Quantity
            Detector integration time in seconds, must be >= 10 s for non-
            destructive reading (NDR) as per the official ESO HAWK-I
            instrument description.
        n_dit : int
            Number of integrations, must be >= 1.

        Raises
        ------
        ValueError
            Raised if `dit` or `n_dit` are outside the allowed range.

        Returns
        -------
        s_n : float
            Signal/noise ratio (dimensionless).

        """
        if dit < 10 * u.s:
            raise ValueError("DIT must be at least 10 s for NDR.")
        if n_dit < 1:
            raise ValueError("NDIT must be at least 1.")

        signal = np.sqrt(n_dit) * target_counts
        noise = (self._get_noise(target_counts, sky_counts_pp, dit))**(1/2)
        s_n = float((signal / noise).value)
        return s_n

    def _to_electrons(self, flux: Quantity[u.photon/u.s/u.m**(-2)],
                      dit: Quantity[u.s]) -> Quantity[u.electron]:
        # HACK: The following conversion factor has been derived to match the
        #       values to those produced by the official ESO HAWK-I ETC, it
        #       has no physical derivation! When used with AB magnitudes,
        #       the factor is closer to 0.4...
        eff = .3581 * u.electron / u.photon
        counts = flux * dit * self.aperture_area * eff
        return counts

    @lru_cache
    def create_sky(self, dit: Quantity[u.s],
                   **kwargs) -> Quantity[u.electron/u.pixel]:
        """
        Estimate number of electrons per pixel expected from sky background.

        Parameters
        ----------
        pixel_scale : astropy.Quantity
            Pixel scale of the optical system in pixel per square meter.
        dit : astropy.Quantity
            Detector integration time in seconds.
        **kwargs : dict
            DESCRIPTION.

        Raises
        ------
        ValueError
            Raised if `dit` is outside the allowed range.

        Returns
        -------
        sky_count : astropy.Quantity
            Sky level in electrons per pixel.

        """
        # TODO: investigate if this could be sped up by a decorator storing the
        #       output of last few calls for given parameters...
        if dit < 10 * u.s:
            raise ValueError("DIT must be at least 10 s for NDR.")
        defaults = {"filter_name": "Ks",
                    "instrument": "HAWKI",
                    "observatory": "Paranal",
                    "airmass": 1.5, "pwv": 10.0}
        therm_dict = {"incl_therm": "Y", "therm_t1": 285.0, "therm_e1": 0.20,
                      "therm_t2": 288.0, "therm_e2": 0.10,  "therm_t3": 33.0,
                      "therm_e3": 0.01}
        parameters = defaults | therm_dict | kwargs
        sky = hmbp.in_skycalc_background(**parameters)
        sky_count = sky * dit / self.pixel_scale * u.electron/u.photon
        return sky_count

    def sn_for_mag(self, mag: Quantity[u.mag],
                   dit: Quantity[u.s], n_dit: int) -> float:
        mag = mag << u.mag
        sky_counts_pp = self.create_sky(dit)
        flux = hmbp.for_flux_in_filter(flux=mag, **calc.default_kwargs)
        target_counts = self._to_electrons(flux, dit)
        print(target_counts)
        s_n = self.signal_to_noise(target_counts, sky_counts_pp, dit, n_dit)
        return s_n

    def _sn_for_mag_for_root(self, mag: Quantity[u.mag],
                             dit: Quantity[u.s], n_dit: int,
                             limit_sn: float = 5.) -> float:
        mag = mag << u.mag
        # print(mag)
        mag = mag[0]
        with HiddenPrints():
            sky_counts_pp = self.create_sky(dit)
        flux = hmbp.for_flux_in_filter(flux=mag, **calc.default_kwargs)
        target_counts = self._to_electrons(flux, dit)
        s_n = self.signal_to_noise(target_counts, sky_counts_pp, dit, n_dit)
        s_n -= limit_sn
        return s_n

    def mag_for_sn(self, s_n: float,
                   dit: Quantity[u.s], n_dit: int) -> Quantity[u.mag]:
        mag = fsolve(self._sn_for_mag_for_root, 20, (dit, n_dit, s_n))[0]
        mag = mag << u.mag
        return mag


if __name__ == "__main__":
    # dit = 60 * u.s
    # n_dit = 60
    # assert dit * n_dit == 1 * u.h  # total observing time should be 1 h

    calc = HawkiEtc()
    sn = calc.sn_for_mag(19.5*u.mag, 60*u.s, 60)
    print(sn)
    limmag = calc.mag_for_sn(50., 60*u.s, 60)
    print(limmag)
