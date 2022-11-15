#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test suite for etc.

Created on Sat Nov 12 16:30:00 2022

@author: teutoburg
"""

__version__ = "0.3"

from unittest import main, TestCase

from astropy import units as u

import etc
from utils import HiddenPrints


class TestSignalToNoise(TestCase):
    def setUp(self):
        self.dit = 60 * u.s
        self.n_dit = 60
        self.target_counts = 2859.78 * u.electron
        self.sky_counts_pp = 165117.59 * u.electron/u.pixel

        self.etc = etc.HawkiEtc()

    def test_sn(self):
        desired = 5.1711
        actual = self.etc.signal_to_noise(self.target_counts,
                                          self.sky_counts_pp,
                                          self.dit, self.n_dit)
        self.assertAlmostEqual(desired, actual, 3)

    def test_invalid_dit(self):
        self.assertRaises(ValueError, self.etc.signal_to_noise,
                          self.target_counts, self.sky_counts_pp,
                          3 * u.s, self.n_dit)

    def test_invalid_n_dit(self):
        self.assertRaises(ValueError, self.etc.signal_to_noise,
                          self.target_counts, self.sky_counts_pp,
                          self.dit, -5)


class TestSky(TestCase):
    def setUp(self):
        self.dit = 60 * u.s
        self.etc = etc.HawkiEtc()

    def test_value(self):
        desired = 165117.562
        with HiddenPrints():
            actual = self.etc.create_sky(self.dit).value
        self.assertAlmostEqual(desired, actual, 3)

    def test_units(self):
        desired = u.Unit("electron / pix")
        with HiddenPrints():
            actual = self.etc.create_sky(self.dit).unit
        self.assertEqual(desired, actual)

    def test_invalid_dit(self):
        self.assertRaises(ValueError, self.etc.create_sky, 3 * u.s)


class TestElectrons(TestCase):
    def setUp(self):
        self.dit = 60 * u.s
        self.etc = etc.HawkiEtc()

    def test_value(self):
        fluxes = (2.52, 15.9, 100.3) * u.photon/u.s/u.m**2
        desireds = (2859., 18041., 113808.)
        for flux, desired in zip(fluxes, desireds):
            with self.subTest():
                actual = self.etc._to_electrons(flux, self.dit).value
                self.assertAlmostEqual(desired, actual, 0)

    def test_units(self):
        desired = u.Unit("electron")
        actual = self.etc._to_electrons(2.52 * u.photon/u.s/u.m**2,
                                        self.dit).unit
        self.assertEqual(desired, actual)


class TestSnMag(TestCase):
    def setUp(self):
        self.dit = 60 * u.s
        self.n_dit = 60
        self.etc = etc.HawkiEtc()

    def test_value(self):
        desired = 12.99
        lim_mag = 21 * u.mag
        with HiddenPrints():
            actual = self.etc.sn_for_mag(lim_mag, self.dit, self.n_dit)
        self.assertAlmostEqual(desired, actual, 2)


class TestMagSn(TestCase):
    def setUp(self):
        self.dit = 60 * u.s
        self.n_dit = 60
        self.etc = etc.HawkiEtc()

    def test_value(self):
        desired = 22.04
        s_n = 5.0
        with HiddenPrints():
            actual = self.etc.mag_for_sn(s_n, self.dit, self.n_dit).value
        self.assertAlmostEqual(desired, actual, 2)


# TODO: are we missing any other tests??

if __name__ == "__main__":
    main()
