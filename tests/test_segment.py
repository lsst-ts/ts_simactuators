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
    def check_from_end_conditions(self, start_position, start_velocity, end_position, end_velocity, dt):
        """Check PathSegment from_end_conditions and limits methods.
        """
        start_tai = 45.3  # any value will do
        end_tai = start_tai + dt
        segment = simactuators.path.PathSegment.from_end_conditions(
            start_tai=start_tai, start_position=start_position, start_velocity=start_velocity,
            end_tai=end_tai, end_position=end_position, end_velocity=end_velocity)
        self.assertAlmostEqual(segment.tai, start_tai)
        self.assertAlmostEqual(segment.pos, start_position)
        self.assertAlmostEqual(segment.vel, start_velocity)

        pva_end = segment.at(end_tai)
        self.assertAlmostEqual(pva_end.pos, end_position)
        self.assertAlmostEqual(pva_end.vel, end_velocity)
        start_accel = segment.accel
        jerk = segment.jerk
        desired_end_position = start_position + dt*(start_velocity + dt*(0.5*start_accel + dt*(1/6)*jerk))
        desired_end_velocity = start_velocity + dt*(start_accel + dt*0.5*jerk)
        self.assertAlmostEqual(pva_end.pos, desired_end_position)
        self.assertAlmostEqual(pva_end.vel, desired_end_velocity)

        # estimate position and velocity limits by computing at many points
        desired_min_position = min(start_position, end_position)
        desired_max_position = max(start_position, end_position)
        desired_max_velocity = max(abs(start_velocity), abs(end_velocity))
        for t in np.linspace(start=0, stop=dt, num=100):
            pt = start_position + t*(start_velocity + t*(0.5*start_accel + t*jerk/6))
            vt = start_velocity + t*(start_accel + t*0.5*jerk)
            desired_min_position = min(pt, desired_min_position)
            desired_max_position = max(pt, desired_max_position)
            desired_max_velocity = max(abs(vt), desired_max_velocity)
        aB = start_accel + dt*jerk
        desired_max_acceleration = max(abs(start_accel), abs(aB))

        limits = segment.limits(end_tai)
        self.assertAlmostEqual(limits.min_position, desired_min_position, places=3)
        self.assertAlmostEqual(limits.max_position, desired_max_position, places=3)
        self.assertAlmostEqual(limits.max_velocity, desired_max_velocity, places=3)
        self.assertAlmostEqual(limits.max_acceleration, desired_max_acceleration)

    def test_from_end_conditions(self):
        for start_position, start_velocity, end_position, end_velocity, dt in itertools.product(
            (0, -0.5, 0.2), (0, -0.2, 0.1), (0, 0.3, -0.6), (0, 0.3, -0.2), (1, 5),
        ):
            with self.subTest(start_position=start_position,
                              start_velocity=start_velocity,
                              end_position=end_position,
                              end_velocity=end_velocity,
                              dt=dt):
                self.check_from_end_conditions(dt=dt,
                                               start_position=start_position,
                                               start_velocity=start_velocity,
                                               end_position=end_position,
                                               end_velocity=end_velocity)

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
            with self.subTest(dt=dt):
                start_tai = 51.0  # arbitrary
                end_tai = start_tai + dt

                with self.assertRaises(ValueError):
                    simactuators.path.PathSegment.from_end_conditions(
                        start_tai, start_position=1, start_velocity=2,
                        end_tai=end_tai, end_position=1, end_velocity=2)


if __name__ == '__main__':
    unittest.main()
