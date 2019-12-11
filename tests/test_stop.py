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
import unittest

from lsst.ts import simactuators


class TestStop(unittest.TestCase):
    def check_path(self, path, tai, position, velocity, max_acceleration):
        """Check various aspects of a path

        Checks the following:

        - The initial time is correct.
        - Times increase monotonically.
        - The position and velocity at the end of each segment
          matches the start of the next.
        - The final position and velocity are correct.

        Parameters
        ----------
        path : `Path`
            Path to check
        tai : `float`
            TAI time (unix seconds, e.g. from time.time())
        position : `float` (optional)
            Position at ``tai`` (deg)
        velocity : `float` (optional)
            Velocity at ``tai`` (deg/sec)
        max_acceleration : `float` (optional)
            Maximum allowed acceleration (deg/sec^2)
        """
        self.assertAlmostEqual(path[0].tai, tai)
        self.assertAlmostEqual(path[0].position, position)
        self.assertAlmostEqual(path[0].velocity, velocity)

        self.assertIn(len(path), (1, 2))

        self.assertEqual(path[-1].velocity, 0)
        self.assertEqual(path[-1].acceleration, 0)

        if len(path) > 1:
            segment0 = path[0]
            segment1 = path[1]
            dt = segment1.tai - segment0.tai
            self.assertGreater(dt, 0)
            pred_p1 = segment0.position + dt*(segment0.velocity + dt*0.5*segment0.acceleration)
            pred_v1 = segment0.velocity + dt*segment0.acceleration
            self.assertAlmostEqual(segment1.position, pred_p1, places=4)
            self.assertAlmostEqual(segment1.velocity, pred_v1, places=4)

    def test_slew_to_stop(self):
        tai = 1550000000
        max_acceleration = 10

        for position, velocity in itertools.product(
            (-5, 0, 30),
            (-3, -1, 2, 4),
        ):
            with self.subTest(position=position, velocity=velocity):
                path = simactuators.path.stop(tai=tai,
                                              position=position,
                                              velocity=velocity,
                                              max_acceleration=max_acceleration)
                self.assertEqual(path.kind, simactuators.path.Kind.Stopping)
                self.assertEqual(len(path), 2)
                self.check_path(path,
                                tai=tai,
                                position=position,
                                velocity=velocity,
                                max_acceleration=max_acceleration)

    def test_already_stopped(self):
        """Test stop when already stopped."""
        # Arbitrary but reasonable values
        tai = 1550000000
        max_acceleration = 2

        for position in (-5, 0, 30):
            path = simactuators.path.stop(tai=tai,
                                          position=position,
                                          velocity=0,
                                          max_acceleration=max_acceleration)
            self.assertEqual(len(path), 1)
            self.assertEqual(path.kind, simactuators.path.Kind.Stopped)
            self.check_path(path,
                            tai=tai,
                            position=position,
                            velocity=0,
                            max_acceleration=max_acceleration)

    def test_invalid_inputs(self):
        # Arbitrary but reasonable values
        tai = 1530000000
        position = 1
        velocity = -2

        # max_acceleration must be >= 0
        with self.assertRaises(ValueError):
            simactuators.path.stop(tai=tai,
                                   position=position,
                                   velocity=velocity,
                                   max_acceleration=0)
        with self.assertRaises(ValueError):
            simactuators.path.stop(tai=tai,
                                   position=position,
                                   velocity=velocity,
                                   max_acceleration=-1)


if __name__ == '__main__':
    unittest.main()
