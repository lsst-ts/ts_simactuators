# This file is part of ts_simactuators.
#
# Developed for the Rubin Observatory Telescope and Site System.
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
    def check_path(
        self,
        path,
        tai,
        start_position,
        start_velocity,
        end_position,
        end_velocity,
        max_velocity,
        max_acceleration,
    ):
        """Check various aspects of a path created by `slew`.

        Parameters
        ----------
        path : `Path`
            Path to check
        tai : `float`
            Start time of slew
        start_position : `float`
            Initial position of path (deg)
        start_position : `float`
            Position of A at time tai (deg)
        start_velocity : `float`
            Velocity of A at time tai (deg/sec)
        end_position : `float`
            Position of B at time tai (deg)
        end_velocity : `float`
            Velocity of B at time tai (deg/sec)
        max_velocity : `float`
            Maximum allowed velocity (deg/sec)
        max_acceleration : `float`
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
        self.assertAlmostEqual(path[0].position, start_position)
        self.assertAlmostEqual(path[0].velocity, start_velocity)

        for i in range(len(path) - 1):
            segment0 = path[i]
            segment1 = path[i + 1]
            dt = segment1.tai - segment0.tai
            self.assertGreater(dt, 0)
            pred_p1 = segment0.position + dt * (
                segment0.velocity + dt * 0.5 * segment0.acceleration
            )
            pred_v1 = segment0.velocity + dt * segment0.acceleration
            self.assertAlmostEqual(segment1.position, pred_p1, places=4)
            self.assertAlmostEqual(segment1.velocity, pred_v1, places=4)

        for segment in path:
            self.assertLessEqual(abs(segment.velocity), max_velocity)
            self.assertLessEqual(abs(segment.acceleration), max_acceleration)

    def test_no_slew(self):
        """Test moving from a point to itself (no slew needed)."""
        # Arbitrary but reasonable values
        tai = 1550000000
        max_velocity = 3
        max_acceleration = 2

        for start_position, start_velocity in itertools.product(
            (-5, 0, 30),
            (0, 1),
        ):
            with self.subTest(
                start_position=start_position, start_velocity=start_velocity
            ):
                path = simactuators.path.slew(
                    tai=tai,
                    start_position=start_position,
                    start_velocity=start_velocity,
                    end_position=start_position,
                    end_velocity=start_velocity,
                    max_velocity=max_velocity,
                    max_acceleration=max_acceleration,
                )
                self.assertEqual(path.kind, simactuators.path.Kind.Slewing)
                self.check_path(
                    path,
                    tai=tai,
                    start_position=start_position,
                    start_velocity=start_velocity,
                    end_position=start_position,
                    end_velocity=start_velocity,
                    max_velocity=max_velocity,
                    max_acceleration=max_acceleration,
                )
                self.assertEqual(len(path), 1)
                self.assertAlmostEqual(path[0].acceleration, 0)

    def test_long_fixed_points(self):
        """Test a fixed point to fixed point slew that is long enough
        to have a segment with constant velocity=+/-max_velocity

        This case is trivial to guess the required answer.
        """
        # Arbitrary but reasonable values
        tai = 1540000000
        max_velocity = 3.5
        max_acceleration = 2.1

        # compute expected delta time and distance covered
        # going from 0 velocity to full velocity at full acceleration
        dt_max_velocity = max_velocity / max_acceleration
        dp_max_velocity = 0.5 * max_acceleration * dt_max_velocity ** 2

        for start_position, dpos in itertools.product(
            (-5, 0),
            (-2.001 * dp_max_velocity, 10 * dp_max_velocity),
        ):
            with self.subTest(start_position=start_position, dpos=dpos):
                # time enough to ramp up to full speed
                # and stay there for at least a short time
                end_position = start_position + dpos
                path = simactuators.path.slew(
                    tai=tai,
                    start_position=start_position,
                    start_velocity=0,
                    end_position=end_position,
                    end_velocity=0,
                    max_velocity=max_velocity,
                    max_acceleration=max_acceleration,
                )
                self.check_path(
                    path,
                    tai=tai,
                    start_position=start_position,
                    start_velocity=0,
                    end_position=end_position,
                    end_velocity=0,
                    max_velocity=max_velocity,
                    max_acceleration=max_acceleration,
                )
                self.assertEqual(path.kind, simactuators.path.Kind.Slewing)
                self.assertEqual(len(path), 4)

                self.assertEqual(len(path), 4)
                self.assertAlmostEqual(path[0].tai, tai)
                self.assertAlmostEqual(path[0].position, start_position)
                self.assertAlmostEqual(path[0].velocity, 0)
                self.assertAlmostEqual(
                    path[0].acceleration, math.copysign(max_acceleration, dpos)
                )

                predicted_dt1 = dt_max_velocity
                predicted_t1 = tai + predicted_dt1
                predicted_dp1 = math.copysign(dp_max_velocity, dpos)
                predicted_p1 = start_position + predicted_dp1
                self.assertAlmostEqual(path[1].tai, predicted_t1, places=4)
                self.assertAlmostEqual(path[1].position, predicted_p1)
                self.assertAlmostEqual(abs(path[1].velocity), max_velocity)
                self.assertAlmostEqual(path[1].acceleration, 0)

                predicted_abs_dp2 = abs(dpos) - 2 * dp_max_velocity
                predicted_dp2 = math.copysign(predicted_abs_dp2, dpos)
                predicted_p2 = path[1].position + predicted_dp2
                predicted_dt2 = abs(predicted_dp2) / max_velocity
                predicted_t2 = path[1].tai + predicted_dt2
                self.assertAlmostEqual(path[2].tai, predicted_t2, places=4)
                self.assertAlmostEqual(path[2].position, predicted_p2)
                self.assertAlmostEqual(abs(path[2].velocity), max_velocity)
                self.assertAlmostEqual(path[2].acceleration, -path[0].acceleration)

                predicted_duration = 2 * dt_max_velocity + predicted_dt2
                predicted_t3 = tai + predicted_duration
                self.assertAlmostEqual(path[3].tai, predicted_t3, places=4)
                self.assertAlmostEqual(path[3].position, end_position)
                self.assertAlmostEqual(path[3].velocity, 0)
                self.assertAlmostEqual(path[3].acceleration, 0)

    def test_short_fixed_points(self):
        """Test a fixed point to fixed point slew that is long enough
        to have a segment with constant velocity=+/-max_velocity

        This case is trivial to guess the required answer.
        """
        # Arbitrary but reasonable values
        tai = 1560000000
        max_velocity = 3.5
        max_acceleration = 2.1

        # compute expected delta time and distance covered
        # going from 0 velocity to full velocity at full acceleration
        dt_max_velocity = max_velocity / max_acceleration
        dp_max_velocity = 0.5 * max_acceleration * dt_max_velocity ** 2

        for start_position, dpos in itertools.product(
            (-5, 0),
            (0.1 * dp_max_velocity, -0.9 * dp_max_velocity),
        ):
            with self.subTest(start_position=start_position, dpos=dpos):
                # not enough time to ramp up to full speed
                end_position = start_position + dpos
                path = simactuators.path.slew(
                    tai=tai,
                    start_position=start_position,
                    start_velocity=0,
                    end_position=end_position,
                    end_velocity=0,
                    max_velocity=max_velocity,
                    max_acceleration=max_acceleration,
                )
                self.assertEqual(path.kind, simactuators.path.Kind.Slewing)
                self.check_path(
                    path,
                    tai=tai,
                    start_position=start_position,
                    start_velocity=0,
                    end_position=end_position,
                    end_velocity=0,
                    max_velocity=max_velocity,
                    max_acceleration=max_acceleration,
                )

                self.assertEqual(len(path), 3)
                self.assertAlmostEqual(path[0].tai, tai)
                self.assertAlmostEqual(path[0].position, start_position)
                self.assertAlmostEqual(path[0].velocity, 0)
                self.assertAlmostEqual(
                    path[0].acceleration, math.copysign(max_acceleration, dpos)
                )

                predicted_dt1 = math.sqrt(abs(dpos) / max_acceleration)
                predicted_t1 = tai + predicted_dt1
                predicted_p1 = start_position + dpos / 2
                predicted_v1 = predicted_dt1 * max_acceleration
                self.assertAlmostEqual(path[1].tai, predicted_t1, places=4)
                self.assertAlmostEqual(path[1].position, predicted_p1)
                self.assertAlmostEqual(abs(path[1].velocity), predicted_v1)
                self.assertAlmostEqual(path[1].acceleration, -path[0].acceleration)

                predicted_t2 = tai + 2 * predicted_dt1
                self.assertAlmostEqual(path[2].tai, predicted_t2, places=4)
                self.assertAlmostEqual(path[2].position, end_position)
                self.assertAlmostEqual(path[2].velocity, 0)
                self.assertAlmostEqual(path[2].acceleration, 0)

    def test_other_slews(self):
        # Arbitrary but reasonable values
        tai = 1560000000
        max_velocity = 3.1
        max_acceleration = 1.76
        dt_max_velocity = max_velocity / max_acceleration

        for start_position, dpos, start_velocity, dvel in itertools.product(
            (-5, 0),
            (dt_max_velocity * 0.1, dt_max_velocity * 10),
            (-max_velocity, -2, 0, 1, max_velocity),
            (-1, 0, 2),
        ):
            with self.subTest(
                start_position=start_position,
                dpos=dpos,
                start_velocity=start_velocity,
                dvel=dvel,
            ):
                end_position = start_position + dpos
                end_velocity = start_velocity + dvel
                if abs(end_velocity) > max_velocity / simactuators.path.SLEW_FUDGE:
                    continue
                path = simactuators.path.slew(
                    tai=tai,
                    start_position=start_position,
                    start_velocity=start_velocity,
                    end_position=end_position,
                    end_velocity=end_velocity,
                    max_velocity=max_velocity,
                    max_acceleration=max_acceleration,
                )
                self.assertEqual(path.kind, simactuators.path.Kind.Slewing)
                self.check_path(
                    path,
                    tai=tai,
                    start_position=start_position,
                    start_velocity=start_velocity,
                    end_position=end_position,
                    end_velocity=end_velocity,
                    max_velocity=max_velocity,
                    max_acceleration=max_acceleration,
                )

    def test_invalid_inputs(self):
        # Arbitrary but reasonable values
        tai = 1530000000
        max_velocity = 3.1
        max_acceleration = 1.76
        start_position = 1
        end_position = 2
        start_velocity = -2
        end_velocity = 3
        # Local version of SLEW_FUDGE with a bit of margin to avoid
        # test failure due to roundoff error
        fudge = simactuators.path.SLEW_FUDGE * 1.000001

        # max_velocity must be >= 0
        with self.assertRaises(ValueError):
            simactuators.path.slew(
                tai=tai,
                start_position=start_position,
                start_velocity=start_velocity,
                end_position=end_position,
                end_velocity=end_velocity,
                max_velocity=0,
                max_acceleration=max_acceleration,
            )
        with self.assertRaises(ValueError):
            simactuators.path.slew(
                tai=tai,
                start_position=start_position,
                start_velocity=start_velocity,
                end_position=end_position,
                end_velocity=end_velocity,
                max_velocity=-1,
                max_acceleration=max_acceleration,
            )

        # max_acceleration must be >= 0
        with self.assertRaises(ValueError):
            simactuators.path.slew(
                tai=tai,
                start_position=start_position,
                start_velocity=start_velocity,
                end_position=end_position,
                end_velocity=end_velocity,
                max_velocity=max_velocity,
                max_acceleration=0,
            )
        with self.assertRaises(ValueError):
            simactuators.path.slew(
                tai=tai,
                start_position=start_position,
                start_velocity=start_velocity,
                end_position=end_position,
                end_velocity=end_velocity,
                max_velocity=max_velocity,
                max_acceleration=-1,
            )

        for sign in (-1, 1):
            # |start_velocity| must be < max_velocity*FUDDGE
            with self.assertRaises(ValueError):
                simactuators.path.slew(
                    tai=tai,
                    start_position=start_position,
                    start_velocity=sign * max_velocity * fudge,
                    end_position=end_position,
                    end_velocity=end_velocity,
                    max_velocity=max_velocity,
                    max_acceleration=max_acceleration,
                )
            with self.assertRaises(ValueError):
                simactuators.path.slew(
                    tai=tai,
                    start_position=start_position,
                    start_velocity=sign * max_velocity * 2,
                    end_position=end_position,
                    end_velocity=end_velocity,
                    max_velocity=max_velocity,
                    max_acceleration=max_acceleration,
                )
            # |end_velocity| must be < max_velocity/FUDDGE
            with self.assertRaises(ValueError):
                simactuators.path.slew(
                    tai=tai,
                    start_position=start_position,
                    start_velocity=start_velocity,
                    end_position=end_position,
                    end_velocity=sign * max_velocity * fudge,
                    max_velocity=max_velocity,
                    max_acceleration=max_acceleration,
                )


if __name__ == "__main__":
    unittest.main()
