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

__version__ = "0.4"

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
    """
    Simplified Exposure Time Calculator class for the ESO HAWK-I instrument.

    Two main functions are currently implemented as methods of this class (see
    also below for example usage), ``sn_for_mag`` and ``mag_for_sn``, which
    will estimate signal-to-noise ratio for a given magnitude, or the limiting
    magnitude for a given signal-to-noise ratio respectively.

    The constructor for this class currently does not any parameters.

    All values which are not dimensionless make heavy use of ``astropy`` units.
    This is documented in all public methods.


    Limitations
    -----------
    Currently, only the K(s) band is tested, and used as a default value for
    ``filter_name``, along with the arguments ``"instrument" = "HAWKI"`` and
    ``"observatory" = "Paranal"``. These values can be accessed via the class
    attribute ``HawkiEtc.default_kwargs``, but behavior with other values is
    not tested. Similarly, the aperture diameter and detector pixel amount are
    fixed to the HAWK-I instrument and the VLT telescope.

    Other current limitations are: the conversion from photons to electrons
    is derived from the official ESO HAWK-I ETC, and lacks a more physical
    formula; the amount of pixels used for noise estimations is fixed to that
    corresponding to a seeing of 0.8 arcsec, which is described as the 50 %
    value of observing conditions for the site.

    Examples
    --------
    Setup instance and exposure time:
        >>> from astropy import units as u
        >>> from etc import HawkiEtc
        >>>
        >>> dit = 60 * u.s  # in seconds, needs Quantity
        >>> n_dit = 60  # simple integer
        >>> # total observing time in this example should be 1 h:
        >>> assert dit * n_dit == 1 * u.h
        >>>
        >>> calc = HawkiEtc()  # create ETC instance

    Calculate S/N ratio from given magnitude:
        >>> lim_mag = 21*u.mag  # only Vega mags currently fully supported!
        >>> sn_ratio = calc.sn_for_mag(lim_mag, dit, n_dit)
        >>> print(f"S/N: {sn_ratio:.2f}")  # formatted printing

    Calculate limiting magnitude from given S/N ratio:
        >>> sn_ratio = 5.0  # dimensionless
        >>> lim_mag = calc.mag_for_sn(sn_ratio, dit, n_dit)
        >>> print(f"{lim_mag:.2f}")  # formatted printing

    """

    default_kwargs = {"filter_name": "Ks",
                      "instrument": "HAWKI",
                      "observatory": "Paranal"}
    aperture = aperture = 8.2 * u.m  # VLT aperture
    pixels = 4*2048*2048 * u.pixel  # HAWK-I detector pixel amount (simplified)

    @property
    def aperture_area(self) -> Quantity[u.m**2]:
        """Aperture area of the telescope. Read-only property."""
        return 0.785398 * (self.aperture**2)  # pi/4 -> circle area

    @property
    def pixel_scale(self) -> Quantity[u.pixel/u.m**2]:
        """Pixel scale of the telescope. Read-only property."""
        return self.pixels / self.aperture_area

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

        signal = n_dit**(1/2) * target_counts
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
            Any keyword arguments taken by ``hmbp.in_skycalc_background``, e.g.
            `airmass` or `pwv`.

        Raises
        ------
        ValueError
            Raised if `dit` is outside the allowed range.

        Returns
        -------
        sky_count : astropy.Quantity
            Sky level in electrons per pixel.

        """
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
        """
        Calculate S/N ratio from given magnitude and integration time.

        Magnitude input `mag` is assumed to be Vega magnitudes. Float input
        will be converted to mag. AB magnitudes are currently not tested.

        Parameters
        ----------
        mag : astropy.Quantity
            Input (Vega) magnitude.
        dit : astropy.Quantity
            Detector integration time in seconds, must be >= 10 s for non-
            destructive reading (NDR) as per the official ESO HAWK-I
            instrument description.
        n_dit : int
            Number of integrations, must be >= 1.

        Returns
        -------
        s_n : float
            Signal-to-noise ratio (dimensionless).

        """
        mag = mag << u.mag
        sky_counts_pp = self.create_sky(dit)
        flux = hmbp.for_flux_in_filter(flux=mag, **self.default_kwargs)
        target_counts = self._to_electrons(flux, dit)
        s_n = self.signal_to_noise(target_counts, sky_counts_pp, dit, n_dit)
        return s_n

    def _sn_for_mag_for_root(self, mag: Quantity[u.mag],
                             dit: Quantity[u.s], n_dit: int,
                             limit_sn: float = 5.) -> float:
        mag = mag << u.mag
        mag = mag[0]
        with HiddenPrints():
            sky_counts_pp = self.create_sky(dit)
        flux = hmbp.for_flux_in_filter(flux=mag, **self.default_kwargs)
        target_counts = self._to_electrons(flux, dit)
        s_n = self.signal_to_noise(target_counts, sky_counts_pp, dit, n_dit)
        s_n -= limit_sn
        return s_n

    def mag_for_sn(self, s_n: float,
                   dit: Quantity[u.s], n_dit: int) -> Quantity[u.mag]:
        """
        Calculate limiting magnitude from given S/N ratio and integration time.

        Result is reached iteratively via solver. This may not always converge.
        Outout limiting magnitude is given in Vega magnitudes.

        Parameters
        ----------
        s_n : float
            Input signal-to-noise ratio (dimensionless), must be > 0.0.
        dit : astropy.Quantity
            Detector integration time in seconds, must be >= 10 s for non-
            destructive reading (NDR) as per the official ESO HAWK-I
            instrument description.
        n_dit : int
            Number of integrations, must be >= 1.

        Raises
        ------
        ValueError
            Raised if `s_n` is outside the allowed range.

        Returns
        -------
        mag : astropy.Quantity
            Output limiting (Vega) magnitude.

        """
        if s_n < 0:
            raise ValueError("S/N ratio must be greater than 0.0.")

        mag = fsolve(self._sn_for_mag_for_root, 20, (dit, n_dit, s_n))[0]
        mag = mag << u.mag
        return mag
