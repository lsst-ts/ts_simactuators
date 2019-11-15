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
    def check_path(self, path, tai, pos, vel, max_accel):
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
        pos : `float` (optional)
            Position at ``tai`` (deg)
        vel : `float` (optional)
            Velocity at ``tai`` (deg/sec)
        max_accel : `float` (optional)
            Maximum allowed acceleration (deg/sec^2)
        """
        self.assertAlmostEqual(path[0].tai, tai)
        self.assertAlmostEqual(path[0].pos, pos)
        self.assertAlmostEqual(path[0].vel, vel)

        self.assertIn(len(path), (1, 2))

        self.assertEqual(path[-1].vel, 0)
        self.assertEqual(path[-1].accel, 0)

        if len(path) > 1:
            pvat0 = path[0]
            pvat1 = path[1]
            dt = pvat1.tai - pvat0.tai
            self.assertGreater(dt, 0)
            pred_p1 = pvat0.pos + dt*(pvat0.vel + dt*0.5*pvat0.accel)
            pred_v1 = pvat0.vel + dt*pvat0.accel
            self.assertAlmostEqual(pvat1.pos, pred_p1, places=4)
            self.assertAlmostEqual(pvat1.vel, pred_v1, places=4)

    def test_slew_to_stop(self):
        tai = 1550000000
        max_accel = 10

        for pos, vel in itertools.product(
            (-5, 0, 30), (-3, -1, 2, 4),
        ):
            path = simactuators.path.stop(tai=tai, pos=pos,
                                          vel=vel, max_accel=max_accel)
            self.assertEqual(path.kind, simactuators.path.Kind.Stopping)
            self.assertEqual(len(path), 2)
            self.check_path(path, tai=tai,
                            pos=pos, vel=vel, max_accel=max_accel)

    def test_already_stopped(self):
        """Test stop when already stopped."""
        # Arbitrary but reasonable values
        tai = 1550000000
        max_accel = 2

        for pos in (-5, 0, 30):
            path = simactuators.path.stop(tai=tai, pos=pos,
                                          vel=0, max_accel=max_accel)
            self.assertEqual(len(path), 1)
            self.assertEqual(path.kind, simactuators.path.Kind.Stopped)
            self.check_path(path, tai=tai,
                            pos=pos, vel=0,
                            max_accel=max_accel)

    def test_invalid_inputs(self):
        # Arbitrary but reasonable values
        tai = 1530000000
        pos = 1
        vel = -2

        # max_accel must be >= 0
        with self.assertRaises(ValueError):
            simactuators.path.stop(tai=tai, pos=pos,
                                   vel=vel, max_accel=0)
        with self.assertRaises(ValueError):
            simactuators.path.stop(tai=tai, pos=pos,
                                   vel=vel, max_accel=-1)


if __name__ == '__main__':
    unittest.main()
