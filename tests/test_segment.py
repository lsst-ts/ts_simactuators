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
        start_tai = 45.3  # any value will do
        end_tai = start_tai + dt
        segment = simactuators.path.PathSegment.from_end_conditions(
            start_tai=start_tai, start_pos=start_pos, start_vel=start_vel,
            end_tai=end_tai, end_pos=end_pos, end_vel=end_vel)
        self.assertAlmostEqual(segment.tai, start_tai)
        self.assertAlmostEqual(segment.pos, start_pos)
        self.assertAlmostEqual(segment.vel, start_vel)

        pva_end = segment.at(end_tai)
        self.assertAlmostEqual(pva_end.pos, end_pos)
        self.assertAlmostEqual(pva_end.vel, end_vel)
        start_accel = segment.accel
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

        limits = segment.limits(end_tai)
        self.assertAlmostEqual(limits.min_pos, desired_min_pos, places=3)
        self.assertAlmostEqual(limits.max_pos, desired_max_pos, places=3)
        self.assertAlmostEqual(limits.max_vel, desired_max_vel, places=3)
        self.assertAlmostEqual(limits.max_accel, desired_max_accel)

    def test_from_end_conditions(self):
        for start_pos, start_vel, end_pos, end_vel, dt in itertools.product(
            (0, -0.5, 0.2), (0, -0.2, 0.1), (0, 0.3, -0.6), (0, 0.3, -0.2), (1, 5),
        ):
            with self.subTest(start_pos=start_pos, start_vel=start_vel,
                              end_pos=end_pos, end_vel=end_vel, dt=dt):
                self.check_from_end_conditions(dt=dt, start_pos=start_pos, start_vel=start_vel,
                                               end_pos=end_pos, end_vel=end_vel)

    def test_basics(self):
        tai = 1573847242.4
        pos = 5.1
        vel = 4.3
        accel = 0.23
        jerk = 0.05
        segment = simactuators.path.PathSegment(tai=tai, pos=pos, vel=vel, accel=accel, jerk=jerk)
        self.assertEqual(segment.tai, tai)
        self.assertEqual(segment.pos, pos)
        self.assertEqual(segment.vel, vel)
        self.assertEqual(segment.accel, accel)
        self.assertEqual(segment.jerk, jerk)

        # The at command at the same time should return a copy
        copied_segment = segment.at(tai)
        self.assertEqual(copied_segment.tai, tai)
        self.assertAlmostEqual(copied_segment.pos, pos)
        self.assertAlmostEqual(copied_segment.vel, vel)
        self.assertAlmostEqual(copied_segment.accel, accel)
        self.assertEqual(copied_segment.jerk, jerk)

        for dt in (-5.1, -3.23, 1.23, 6.6):
            with self.subTest(dt=dt):
                new_tai = tai + dt
                new_segment = segment.at(new_tai)
                # Note: I am using slightly cruder math here than in
                # `PathSegment` (to avoid an exact copy),
                # so the computed values will be slightly different,
                # especially position and to a lessar extent velocity.
                expected_pos = pos + vel*dt + accel*dt*dt/2 + jerk*dt*dt*dt/6
                expected_vel = vel + accel*dt + jerk*dt*dt/2
                expected_accel = accel + jerk*dt
                self.assertEqual(new_segment.tai, new_tai)
                self.assertAlmostEqual(new_segment.pos, expected_pos, places=5)
                self.assertAlmostEqual(new_segment.vel, expected_vel, places=6)
                self.assertAlmostEqual(new_segment.accel, expected_accel)
                self.assertEqual(new_segment.jerk, jerk)

    def test_default_arguments(self):
        """Test constructing a PathSegment with default arguments.
        """
        full_kwargs = dict(
            tai=1573847242.4,
            pos=5.1,
            vel=4.3,
            accel=0.23,
            jerk=0.05,
        )
        for field_to_omit in full_kwargs:
            if field_to_omit == "tai":
                continue
            kwargs = full_kwargs.copy()
            del kwargs[field_to_omit]
            segment = simactuators.path.PathSegment(**kwargs)
            for field_name in full_kwargs:
                if field_name == field_to_omit:
                    self.assertEqual(getattr(segment, field_name), 0)
                else:
                    self.assertEqual(getattr(segment, field_name), full_kwargs[field_name])

    def test_invalid_inputs(self):
        """Test invalid inputs for PathSegment.from_end_conditions.

        The only thing that can go wrong is end_tai - start_tai <= 0.
        """
        for dt in (0, math.sqrt(sys.float_info.min)):
            start_tai = 51.0  # arbitrary
            end_tai = start_tai + dt

            with self.assertRaises(ValueError):
                simactuators.path.PathSegment.from_end_conditions(
                    start_tai, start_pos=1, start_vel=2,
                    end_tai=end_tai, end_pos=1, end_vel=2)


if __name__ == '__main__':
    unittest.main()
