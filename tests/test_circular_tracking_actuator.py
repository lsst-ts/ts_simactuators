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

from lsst.ts import salobj
from lsst.ts import simactuators


class TestCircularTrackingActuator(unittest.TestCase):
    def test_constructor(self):
        max_velocity = 3
        max_acceleration = 4
        dtmax_track = 0.5
        nsettle = 1
        tai = 0.5

        actuator = simactuators.CircularTrackingActuator(
            max_velocity=max_velocity,
            max_acceleration=max_acceleration,
            dtmax_track=dtmax_track,
            nsettle=nsettle,
            tai=tai,
        )
        self.assertEqual(actuator.max_velocity, max_velocity)
        self.assertEqual(actuator.max_acceleration, max_acceleration)
        self.assertEqual(actuator.dtmax_track, dtmax_track)
        self.assertEqual(actuator.nsettle, nsettle)
        self.assertEqual(actuator.target.tai, tai)
        self.assertEqual(actuator.target.position, 0)
        self.assertEqual(actuator.target.velocity, 0)
        self.assertEqual(actuator.target.acceleration, 0)
        self.assertEqual(actuator.target.jerk, 0)
        self.assertEqual(actuator.path[0].tai, tai)
        self.assertEqual(actuator.path[0].position, 0)
        self.assertEqual(actuator.path[0].velocity, 0)
        self.assertEqual(actuator.path[0].acceleration, 0)
        self.assertEqual(actuator.path[0].jerk, 0)
        self.assertEqual(len(actuator.path), 1)
        self.assertEqual(actuator.path.kind, actuator.Kind.Stopped)

        for start_position in (-10, -0.001, 1, 359, 360):
            expected_start_position = salobj.angle_wrap_nonnegative(start_position).deg
            actuator = simactuators.CircularTrackingActuator(
                max_velocity=max_velocity,
                max_acceleration=max_acceleration,
                dtmax_track=dtmax_track,
                nsettle=nsettle,
                tai=tai,
                start_position=start_position,
            )
            self.assertAlmostEqual(actuator.path[-1].position, expected_start_position)
            self.assertAlmostEqual(actuator.target.position, expected_start_position)

    def test_constructor_errors(self):
        max_velocity = 3
        max_acceleration = 4
        dtmax_track = 0.5

        # Check that max_velocity must be positive.
        with self.assertRaises(ValueError):
            simactuators.CircularTrackingActuator(
                max_velocity=0,
                max_acceleration=max_acceleration,
                dtmax_track=dtmax_track,
            )
        with self.assertRaises(ValueError):
            simactuators.CircularTrackingActuator(
                max_velocity=-1,
                max_acceleration=max_acceleration,
                dtmax_track=dtmax_track,
            )

        # Check that max_acceleration must be positive.
        with self.assertRaises(ValueError):
            simactuators.CircularTrackingActuator(
                max_velocity=max_velocity, max_acceleration=0, dtmax_track=dtmax_track,
            )
        with self.assertRaises(ValueError):
            simactuators.CircularTrackingActuator(
                max_velocity=max_velocity, max_acceleration=-1, dtmax_track=dtmax_track,
            )

    def test_set_target(self):
        for (
            start_position,
            end_position,
            start_velocity,
            end_velocity,
        ) in itertools.product(
            (0, -0.1, 0.1, -60, 359.9, 360),
            (0, -0.1, 0.1, 300, 359.9, 360),
            (0, -0.1, 0.1),
            (0, -0.1, 0.1),
        ):
            start_tai = salobj.current_tai()
            end_tai = start_tai + 1
            start_segment = simactuators.path.PathSegment(
                position=start_position, velocity=start_velocity, tai=start_tai
            )
            end_segment = simactuators.path.PathSegment(
                position=end_position, velocity=end_velocity, tai=end_tai
            )
            self.check_set_target(start_segment=start_segment, end_segment=end_segment)

    def check_set_target(self, start_segment, end_segment):
        self.assertGreater(end_segment.tai, start_segment.tai)
        max_velocity = 3
        max_acceleration = 4
        dtmax_track = 0.5
        actuator = simactuators.CircularTrackingActuator(
            max_velocity=max_velocity,
            max_acceleration=max_acceleration,
            dtmax_track=dtmax_track,
        )

        targets = dict()
        paths = dict()
        predicted_wrapped_target_position = salobj.angle_wrap_nonnegative(
            end_segment.position
        ).deg
        for direction in simactuators.Direction:
            actuator.target = start_segment
            actuator.path = simactuators.path.Path(
                start_segment, kind=actuator.Kind.Tracking
            )
            actuator.set_target(
                tai=end_segment.tai,
                position=end_segment.position,
                velocity=end_segment.velocity,
                direction=direction,
            )
            self.assertEqual(actuator.target.tai, end_segment.tai)
            self.assertEqual(actuator.target.velocity, end_segment.velocity)
            self.assertAlmostEqual(
                actuator.target.position, predicted_wrapped_target_position
            )
            targets[direction] = actuator.target
            paths[direction] = actuator.path

        # Check that NEAREST direction chooses the path of shortest duration
        end_tais = {direction: path[-1].tai for direction, path in paths.items()}
        min_end_tai = min(tai for tai in end_tais.values())
        self.assertEqual(end_tais[simactuators.Direction.NEAREST], min_end_tai)

        # Check that POSITIVE direction chooses positive direction at
        # end_segment.tai for the target vs the current position at end_tai,
        # and NEGATIVE the negative direction.
        # Use a bit of slop for roundoff error when target = current.
        for direction in (
            simactuators.Direction.NEGATIVE,
            simactuators.Direction.POSITIVE,
        ):
            target_position = targets[direction].at(end_segment.tai).position
            current_position = paths[direction].at(end_segment.tai).position
            slop = 1e-5
            if direction is simactuators.Direction.POSITIVE:
                self.assertGreaterEqual(target_position + slop, current_position)
            else:
                self.assertLessEqual(target_position - slop, current_position)


if __name__ == "__main__":
    unittest.main()
