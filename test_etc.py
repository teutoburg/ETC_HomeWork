#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
TBA.

Created on Sat Nov 12 16:30:00 2022

@author: teutoburg
"""

__version__ = "0.2"

import sys
import os
from unittest import main, TestCase

from astropy import units as u

import etc


class TestSignalToNoise(TestCase):
    def setUp(self):
        self.dit = 60 * u.s
        self.etc = etc.HawkiEtc()

    def test_sn(self):
        desired = 5.1711

        target_counts = 2859.78 * u.electron
        sky_counts_pp = 165117.59 * u.electron/u.pixel
        dit = 60 * u.s
        n_dit = 60

        actual = self.etc.signal_to_noise(target_counts, sky_counts_pp, dit, n_dit)
        self.assertAlmostEqual(desired, actual, 3)


class TestSky(TestCase):
    def setUp(self):
        self.dit = 60 * u.s
        self.etc = etc.HawkiEtc()

        # To silence hmbp internal printing
        sys.stdout = open(os.devnull, 'w')

    def tearDown(self):
        # Restore printing...
        sys.stdout = sys.__stdout__

    def test_value(self):
        desired = 165117.562
        actual = self.etc.create_sky(self.etc.pixel_scale, self.dit).value
        self.assertAlmostEqual(desired, actual, 3)

    def test_units(self):
        desired = u.Unit("electron / pix")
        actual = self.etc.create_sky(self.etc.pixel_scale, self.dit).unit
        self.assertEqual(desired, actual)


if __name__ == "__main__":
    main()
