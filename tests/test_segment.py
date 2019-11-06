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


class TestPathSegment(unittest.TestCase):
    def check_from_end_conditions(self, start_pos, start_vel, end_pos, end_vel, dt):
        """Check PathSegment from_end_conditions and limits methods.
        """
        start_time = 45.3  # any value will do
        end_time = start_time + dt
        segment = simactuators.path.PathSegment.from_end_conditions(
            start_time=start_time, start_pos=start_pos, start_vel=start_vel,
            end_time=end_time, end_pos=end_pos, end_vel=end_vel)
        self.assertAlmostEqual(segment.start_time, start_time)
        self.assertAlmostEqual(segment.start_pos, start_pos)
        self.assertAlmostEqual(segment.start_vel, start_vel)

        pva_end = segment.pva(end_time)
        self.assertAlmostEqual(pva_end.pos, end_pos)
        self.assertAlmostEqual(pva_end.vel, end_vel)
        start_accel = segment.start_accel
        jerk = segment.jerk
        desired_end_pos = start_pos + dt*(start_vel + dt*(0.5*start_accel + dt*(1/6)*jerk))
        desired_end_vel = start_vel + dt*(start_accel + dt*0.5*jerk)
        self.assertAlmostEqual(pva_end.pos, desired_end_pos)
        self.assertAlmostEqual(pva_end.vel, desired_end_vel)

        # estimate position and velocity limits by computing at many points
        desired_min_pos = min(start_pos, end_pos)
        desired_max_pos = max(start_pos, end_pos)
        desired_max_vel = max(abs(start_vel), abs(end_vel))
        for t in np.linspace(start=0, stop=dt, num=100):
            pt = start_pos + t*(start_vel + t*(0.5*start_accel + t*jerk/6))
            vt = start_vel + t*(start_accel + t*0.5*jerk)
            desired_min_pos = min(pt, desired_min_pos)
            desired_max_pos = max(pt, desired_max_pos)
            desired_max_vel = max(abs(vt), desired_max_vel)
        aB = start_accel + dt*jerk
        desired_max_accel = max(abs(start_accel), abs(aB))

        limits = segment.limits(end_time)
        self.assertAlmostEqual(limits.min_pos, desired_min_pos, places=3)
        self.assertAlmostEqual(limits.max_pos, desired_max_pos, places=3)
        self.assertAlmostEqual(limits.max_vel, desired_max_vel, places=3)
        self.assertAlmostEqual(limits.max_accel, desired_max_accel)

    def test_basics(self):
        for start_pos, start_vel, end_pos, end_vel, dt in itertools.product(
            (0, -0.5, 0.2), (0, -0.2, 0.1), (0, 0.3, -0.6), (0, 0.3, -0.2), (1, 5),
        ):
            self.check_from_end_conditions(dt=dt, start_pos=start_pos, start_vel=start_vel,
                                           end_pos=end_pos, end_vel=end_vel)

    def test_invalid_inputs(self):
        """Test invalid inputs for PathSegment.from_end_conditions.

        The only thing that can go wrong is end_time - start_time <= 0.
        """
        for dt in (0, math.sqrt(sys.float_info.min)):
            start_time = 51.0  # arbitrary
            end_time = start_time + dt

            with self.assertRaises(ValueError):
                simactuators.path.PathSegment.from_end_conditions(
                    start_time, start_pos=1, start_vel=2,
                    end_time=end_time, end_pos=1, end_vel=2)


if __name__ == '__main__':
    unittest.main()
