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
import unittest

from lsst.ts import simactuators


class TestSlew(unittest.TestCase):
    def check_path(self, path, tai, start_pos, start_vel, end_pos, end_vel, max_vel, max_accel):
        """Check various aspects of a path created by `slew`.

        Parameters
        ----------
        path : `Path`
            Path to check
        tai : `float`
            Start time of slew
        start_pos : `float`
            Initial position of path (deg)
        start_pos : `float`
            Position of A at time tai (deg)
        start_vel : `float`
            Velocity of A at time tai (deg/sec)
        end_pos : `float`
            Position of B at time tai (deg)
        end_vel : `float`
            Velocity of B at time tai (deg/sec)
        max_vel : `float`
            Maximum allowed velocity (deg/sec)
        max_accel : `float`
            Maximum allowed acceleration (deg/sec^2)

        Notes
        -----

        Checks the following:

        - The initial time is correct.
        - Times increase monotonically.
        - The position and velocity at the end of each segment
          matches the start of the next.
        - The final position and velocity are correct.
        """
        self.assertAlmostEqual(path[0].tai, tai)
        self.assertAlmostEqual(path[0].pos, start_pos)
        self.assertAlmostEqual(path[0].vel, start_vel)

        for i in range(len(path) - 1):
            segment0 = path[i]
            segment1 = path[i+1]
            dt = segment1.tai - segment0.tai
            self.assertGreater(dt, 0)
            pred_p1 = segment0.pos + dt*(segment0.vel + dt*0.5*segment0.accel)
            pred_v1 = segment0.vel + dt*segment0.accel
            self.assertAlmostEqual(segment1.pos, pred_p1, places=4)
            self.assertAlmostEqual(segment1.vel, pred_v1, places=4)

        for tpvaj in path:
            self.assertLessEqual(abs(tpvaj.vel), max_vel)
            self.assertLessEqual(abs(tpvaj.accel), max_accel)

    def test_no_slew(self):
        """Test moving from a point to itself (no slew needed)."""
        # Arbitrary but reasonable values
        tai = 1550000000
        max_vel = 3
        max_accel = 2

        for start_pos in (-5, 0, 30):
            for start_vel in (0, 1):
                path = simactuators.path.slew(tai=tai,
                                              start_pos=start_pos, start_vel=start_vel,
                                              end_pos=start_pos, end_vel=start_vel,
                                              max_vel=max_vel, max_accel=max_accel)
                self.assertEqual(path.kind, simactuators.path.Kind.Slewing)
                self.check_path(path, tai=tai,
                                start_pos=start_pos, start_vel=start_vel,
                                end_pos=start_pos, end_vel=start_vel,
                                max_vel=max_vel, max_accel=max_accel)
                self.assertEqual(len(path), 1)
                self.assertAlmostEqual(path[0].accel, 0)

    def test_long_fixed_points(self):
        """Test a fixed point to fixed point slew that is long enough
        to have a segment with constant velocity=+/-max_vel

        This case is trivial to guess the required answer.
        """
        # Arbitrary but reasonable values
        tai = 1540000000
        max_vel = 3.5
        max_accel = 2.1

        # compute expected delta time and distance covered
        # going from 0 velocity to full velocity at full acceleration
        dt_max_vel = max_vel/max_accel
        dp_max_vel = 0.5*max_accel*dt_max_vel**2

        for start_pos in (-5, 0):
            for dpos in (-2.001*dp_max_vel, 10*dp_max_vel):
                # time enough to ramp up to full speed
                # and stay there for at least a short time
                end_pos = start_pos + dpos
                path = simactuators.path.slew(tai=tai, start_pos=start_pos, start_vel=0,
                                              end_pos=end_pos, end_vel=0,
                                              max_vel=max_vel, max_accel=max_accel)
                self.check_path(path, tai=tai,
                                start_pos=start_pos, start_vel=0,
                                end_pos=end_pos, end_vel=0,
                                max_vel=max_vel, max_accel=max_accel)
                self.assertEqual(path.kind, simactuators.path.Kind.Slewing)
                self.assertEqual(len(path), 4)

                self.assertEqual(len(path), 4)
                self.assertAlmostEqual(path[0].tai, tai)
                self.assertAlmostEqual(path[0].pos, start_pos)
                self.assertAlmostEqual(path[0].vel, 0)
                self.assertAlmostEqual(path[0].accel, math.copysign(max_accel, dpos))

                predicted_dt1 = dt_max_vel
                predicted_t1 = tai + predicted_dt1
                predicted_dp1 = math.copysign(dp_max_vel, dpos)
                predicted_p1 = start_pos + predicted_dp1
                self.assertAlmostEqual(path[1].tai, predicted_t1, places=4)
                self.assertAlmostEqual(path[1].pos, predicted_p1)
                self.assertAlmostEqual(abs(path[1].vel), max_vel)
                self.assertAlmostEqual(path[1].accel, 0)

                predicted_abs_dp2 = abs(dpos) - 2*dp_max_vel
                predicted_dp2 = math.copysign(predicted_abs_dp2, dpos)
                predicted_p2 = path[1].pos + predicted_dp2
                predicted_dt2 = abs(predicted_dp2)/max_vel
                predicted_t2 = path[1].tai + predicted_dt2
                self.assertAlmostEqual(path[2].tai, predicted_t2, places=4)
                self.assertAlmostEqual(path[2].pos, predicted_p2)
                self.assertAlmostEqual(abs(path[2].vel), max_vel)
                self.assertAlmostEqual(path[2].accel, -path[0].accel)

                predicted_duration = 2*dt_max_vel + predicted_dt2
                predicted_t3 = tai + predicted_duration
                self.assertAlmostEqual(path[3].tai, predicted_t3, places=4)
                self.assertAlmostEqual(path[3].pos, end_pos)
                self.assertAlmostEqual(path[3].vel, 0)
                self.assertAlmostEqual(path[3].accel, 0)

    def test_short_fixed_points(self):
        """Test a fixed point to fixed point slew that is long enough
        to have a segment with constant velocity=+/-max_vel

        This case is trivial to guess the required answer.
        """
        # Arbitrary but reasonable values
        tai = 1560000000
        max_vel = 3.5
        max_accel = 2.1

        # compute expected delta time and distance covered
        # going from 0 velocity to full velocity at full acceleration
        dt_max_vel = max_vel/max_accel
        dp_max_vel = 0.5*max_accel*dt_max_vel**2

        for start_pos in (-5, 0):
            for dpos in (0.1*dp_max_vel, -0.9*dp_max_vel):
                # not enough time to ramp up to full speed
                end_pos = start_pos + dpos
                path = simactuators.path.slew(tai=tai,
                                              start_pos=start_pos, start_vel=0,
                                              end_pos=end_pos, end_vel=0,
                                              max_vel=max_vel, max_accel=max_accel)
                self.assertEqual(path.kind, simactuators.path.Kind.Slewing)
                self.check_path(path, tai=tai,
                                start_pos=start_pos, start_vel=0,
                                end_pos=end_pos, end_vel=0,
                                max_vel=max_vel, max_accel=max_accel)

                self.assertEqual(len(path), 3)
                self.assertAlmostEqual(path[0].tai, tai)
                self.assertAlmostEqual(path[0].pos, start_pos)
                self.assertAlmostEqual(path[0].vel, 0)
                self.assertAlmostEqual(path[0].accel, math.copysign(max_accel, dpos))

                predicted_dt1 = math.sqrt(abs(dpos)/max_accel)
                predicted_t1 = tai + predicted_dt1
                predicted_p1 = start_pos + dpos/2
                predicted_v1 = predicted_dt1*max_accel
                self.assertAlmostEqual(path[1].tai, predicted_t1, places=4)
                self.assertAlmostEqual(path[1].pos, predicted_p1)
                self.assertAlmostEqual(abs(path[1].vel), predicted_v1)
                self.assertAlmostEqual(path[1].accel, -path[0].accel)

                predicted_t2 = tai + 2*predicted_dt1
                self.assertAlmostEqual(path[2].tai, predicted_t2, places=4)
                self.assertAlmostEqual(path[2].pos, end_pos)
                self.assertAlmostEqual(path[2].vel, 0)
                self.assertAlmostEqual(path[2].accel, 0)

    def test_other_slews(self):
        # Arbitrary but reasonable values
        tai = 1560000000
        max_vel = 3.1
        max_accel = 1.76
        dt_max_vel = max_vel/max_accel

        for start_pos, dpos, start_vel, dvel in itertools.product(
            (-5, 0), (dt_max_vel*0.1, dt_max_vel*10), (-max_vel, -2, 0, 1, max_vel), (-1, 0, 2),
        ):
            end_pos = start_pos + dpos
            end_vel = start_vel + dvel
            if abs(end_vel) > max_vel/simactuators.path.SLEW_FUDGE:
                continue
            path = simactuators.path.slew(tai=tai,
                                          start_pos=start_pos, start_vel=start_vel,
                                          end_pos=end_pos, end_vel=end_vel,
                                          max_vel=max_vel, max_accel=max_accel)
            self.assertEqual(path.kind, simactuators.path.Kind.Slewing)
            self.check_path(path, tai=tai,
                            start_pos=start_pos, start_vel=start_vel,
                            end_pos=end_pos, end_vel=end_vel,
                            max_vel=max_vel, max_accel=max_accel)

    def test_invalid_inputs(self):
        # Arbitrary but reasonable values
        tai = 1530000000
        max_vel = 3.1
        max_accel = 1.76
        start_pos = 1
        end_pos = 2
        start_vel = -2
        end_vel = 3
        # Local version of SLEW_FUDGE with a bit of margin to avoid
        # test failure due to roundoff error
        fudge = simactuators.path.SLEW_FUDGE*1.000001

        # max_vel must be >= 0
        with self.assertRaises(ValueError):
            simactuators.path.slew(tai=tai,
                                   start_pos=start_pos, start_vel=start_vel,
                                   end_pos=end_pos, end_vel=end_vel,
                                   max_vel=0, max_accel=max_accel)
        with self.assertRaises(ValueError):
            simactuators.path.slew(tai=tai,
                                   start_pos=start_pos, start_vel=start_vel,
                                   end_pos=end_pos, end_vel=end_vel,
                                   max_vel=-1, max_accel=max_accel)

        # max_accel must be >= 0
        with self.assertRaises(ValueError):
            simactuators.path.slew(tai=tai,
                                   start_pos=start_pos, start_vel=start_vel,
                                   end_pos=end_pos, end_vel=end_vel,
                                   max_vel=max_vel, max_accel=0)
        with self.assertRaises(ValueError):
            simactuators.path.slew(tai=tai,
                                   start_pos=start_pos, start_vel=start_vel,
                                   end_pos=end_pos, end_vel=end_vel,
                                   max_vel=max_vel, max_accel=-1)

        for sign in (-1, 1):
            # |start_vel| must be < max_vel*FUDDGE
            with self.assertRaises(ValueError):
                simactuators.path.slew(tai=tai,
                                       start_pos=start_pos, start_vel=sign*max_vel*fudge,
                                       end_pos=end_pos, end_vel=end_vel,
                                       max_vel=max_vel, max_accel=max_accel)
            with self.assertRaises(ValueError):
                simactuators.path.slew(tai=tai,
                                       start_pos=start_pos, start_vel=sign*max_vel*2,
                                       end_pos=end_pos, end_vel=end_vel,
                                       max_vel=max_vel, max_accel=max_accel)
            # |end_vel| must be < max_vel/FUDDGE
            with self.assertRaises(ValueError):
                simactuators.path.slew(tai=tai,
                                       start_pos=start_pos, start_vel=start_vel,
                                       end_pos=end_pos, end_vel=sign*max_vel*fudge,
                                       max_vel=max_vel, max_accel=max_accel)


if __name__ == '__main__':
    unittest.main()
