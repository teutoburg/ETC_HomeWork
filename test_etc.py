#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
TBA.

Created on Sat Nov 12 16:30:00 2022

@author: teutoburg
"""

__version__ = "0.1"

from unittest import main, TestCase

from astropy import units as u

import etc

__author__ = "Fabian Haberhauer"
__copyright__ = "Copyright 2022"
__credits__ = []
__license__ = "GPL"
__maintainer__ = "Fabian Haberhauer"
__email__ = "fabian.haberhauer@univie.ac.at"
__status__ = "Prototype"


class TestSignalToNoise(TestCase):
    def test_something(self):
        desired = 5.1711

        target_counts = 2859.78 * u.electron
        sky_counts_pp = 165117.59 * u.electron/u.pixel
        dit = 60 * u.s
        n_dit = 60

        actual = etc.signal_to_noise(target_counts, sky_counts_pp, dit, n_dit)
        self.assertAlmostEqual(desired, actual, 3)


if __name__ == "__main__":
    main()
