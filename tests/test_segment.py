# This file is part of ts_simactuators.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import itertools
import math
import sys
import unittest

import numpy as np

from lsst.ts import simactuators


class TestSegment(unittest.TestCase):
    def check_segment(self, start_pos, start_vel, end_pos, end_vel, dt):
        """Check that a segment meets is constructor conditions
        """
        segment = simactuators.path.Segment(dt=dt, start_pos=start_pos, start_vel=start_vel,
                                            end_pos=end_pos, end_vel=end_vel, do_pos_lim=True)
        self.assertAlmostEqual(segment.dt, dt)
        self.assertAlmostEqual(segment.start_pos, start_pos)
        self.assertAlmostEqual(segment.end_pos, end_pos)
        self.assertAlmostEqual(segment.start_vel, start_vel)
        self.assertAlmostEqual(segment.end_vel, end_vel)
        start_accel = segment.start_accel
        jerk = segment.jerk
        desired_end_pos = start_pos + dt*(start_vel + dt*(0.5*start_accel + dt*(1/6)*jerk))
        desired_end_vel = start_vel + dt*(start_accel + dt*0.5*jerk)
        self.assertAlmostEqual(segment.end_pos, desired_end_pos)
        self.assertAlmostEqual(segment.end_vel, desired_end_vel)

        # estimate limits by computing at many points
        desired_pmin = min(start_pos, end_pos)
        desired_pmax = max(start_pos, end_pos)
        desired_vpeak = max(abs(start_vel), abs(end_vel))
        for t in np.linspace(start=0, stop=dt, num=100):
            pt = start_pos + t*(start_vel + t*(0.5*start_accel + t*jerk/6))
            vt = start_vel + t*(start_accel + t*0.5*jerk)
            desired_pmin = min(pt, desired_pmin)
            desired_pmax = max(pt, desired_pmax)
            desired_vpeak = max(abs(vt), desired_vpeak)
        self.assertAlmostEqual(segment.min_pos, desired_pmin, places=3)
        self.assertAlmostEqual(segment.max_pos, desired_pmax, places=3)
        self.assertAlmostEqual(segment.peak_vel, desired_vpeak, places=3)

        aB = start_accel + dt*jerk
        desired_apeak = max(abs(start_accel), abs(aB))
        self.assertAlmostEqual(segment.peak_accel, desired_apeak)

    def test_basics(self):
        # start_pos, start_vel, end_pos, end_vel, dt, do_pos_lim
        for start_pos, start_vel, end_pos, end_vel, dt in itertools.product(
            (0, -0.5, 0.2), (0, -0.2, 0.1), (0, 0.3, -0.6), (0, 0.3, -0.2), (1, 5),
        ):
            self.check_segment(dt=dt, start_pos=start_pos, start_vel=start_vel,
                               end_pos=end_pos, end_vel=end_vel)

    def test_invalid_inputs(self):
        for dt in (0, math.sqrt(sys.float_info.min)):
            with self.assertRaises(ValueError):
                simactuators.path.Segment(dt=dt, start_pos=1, start_vel=2,
                                          end_pos=1, end_vel=2, do_pos_lim=True)


if __name__ == '__main__':
    unittest.main()
