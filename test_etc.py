#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
TBA.

Created on Sat Nov 12 16:30:00 2022

@author: teutoburg
"""

__version__ = "0.2"

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


# TODO: add tests for final sn from mag and mag from sn
# TODO: are we missing any other tests??

if __name__ == "__main__":
    main()
