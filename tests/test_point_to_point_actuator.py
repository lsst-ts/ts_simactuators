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

import asyncio
import time
import unittest

import asynctest

from lsst.ts import simactuators


class TestPointToPointActuator(asynctest.TestCase):
    def test_constructor(self):
        min_position = -1
        max_position = 5
        pos = 3
        speed = 1.5
        for good_pos in (pos, min_position, max_position):
            with self.subTest(good_pos=good_pos):
                actuator = simactuators.PointToPointActuator(min_position=min_position,
                                                             max_position=max_position,
                                                             pos=good_pos,
                                                             speed=speed)
                self.assertEqual(actuator.current_position, good_pos)
                self.assertFalse(actuator.moving)
                self.assertEqual(actuator.remaining_time, 0)

        for bad_min_position in (max_position, max_position+0.001):
            with self.subTest(bad_min_position=bad_min_position):
                with self.assertRaises(ValueError):
                    simactuators.PointToPointActuator(min_position=bad_min_position,
                                                      max_position=max_position,
                                                      pos=max_position,
                                                      speed=speed)

        for bad_max_position in (min_position, min_position-0.001):
            with self.subTest(bad_max_position=bad_max_position):
                with self.assertRaises(ValueError):
                    simactuators.PointToPointActuator(min_position=min_position,
                                                      max_position=bad_max_position,
                                                      pos=min_position,
                                                      speed=speed)

        for bad_pos in (min_position-0.001, max_position+0.001):
            with self.subTest(bad_pos=bad_pos):
                with self.assertRaises(ValueError):
                    simactuators.PointToPointActuator(min_position=min_position,
                                                      max_position=max_position,
                                                      pos=bad_pos,
                                                      speed=speed)

        for bad_speed in (0, -0.001):
            with self.subTest(bad_speed=bad_speed):
                with self.assertRaises(ValueError):
                    simactuators.PointToPointActuator(min_position=min_position,
                                                      max_position=max_position,
                                                      pos=pos,
                                                      speed=bad_speed)

    async def test_set_pos(self):
        min_position = -10
        max_position = 10
        pos = 3
        speed = 2
        actuator = simactuators.PointToPointActuator(min_position=min_position,
                                                     max_position=max_position,
                                                     pos=pos,
                                                     speed=speed)
        new_pos = 4
        # Keep track of how long it takes to call set_pos and remaining_time
        # so we know how picky to be when testing remaining_time
        call_start_time = time.monotonic()
        actuator.set_pos(new_pos)
        rem_time = actuator.remaining_time
        time_slop = time.monotonic() - call_start_time
        desired_rem_time = abs(new_pos - pos)/speed
        # It takes finite time to call set_pos and remaining_time;
        # that time slop determines how accurately we can predict
        # the value of remaining_time.
        time_error = abs(desired_rem_time - rem_time)
        self.assertLessEqual(time_error, time_slop)
        self.assertEqual(actuator.end_position, new_pos)
        self.assertTrue(actuator.moving)

        await asyncio.sleep(rem_time + 0.001)
        self.assertFalse(actuator.moving)
        self.assertEqual(actuator.remaining_time, 0)

    async def test_direction(self):
        min_position = -10
        max_position = 10
        pos = 3
        speed = 2
        actuator = simactuators.PointToPointActuator(min_position=min_position,
                                                     max_position=max_position,
                                                     pos=pos,
                                                     speed=speed)

        # If start_position == end_position direction is 1.
        self.assertEqual(actuator.start_position, actuator.end_position)
        self.assertEqual(actuator.direction, 1)

        # Move in the positive direction.
        actuator.set_pos(pos + 1)
        self.assertEqual(actuator.direction, 1)
        # Give the actuator time to move a bit.
        await asyncio.sleep(0.01)
        actuator.stop()
        self.assertEqual(actuator.direction, 1)

        # Move in the negative direction (back to pos).
        actuator.set_pos(pos)
        self.assertEqual(actuator.direction, -1)

    async def test_stop(self):
        min_position = -10
        max_position = 10
        pos = 3
        # Slow motion so plenty of time to stop
        speed = 0.1
        actuator = simactuators.PointToPointActuator(min_position=min_position,
                                                     max_position=max_position,
                                                     pos=pos,
                                                     speed=speed)
        new_pos = 4
        call_start_time = time.monotonic()
        actuator.set_pos(new_pos)
        self.assertEqual(actuator.end_position, new_pos)
        self.assertTrue(actuator.moving)

        move_time = 0.1
        predicted_end_position = pos + speed*move_time
        await asyncio.sleep(move_time)
        actuator.stop()
        call_duration = time.monotonic() - call_start_time
        # It takes a bit longer to call set_pos, sleep and stop
        # than the sleep time. That time slop determines how accurately
        # we can to predict the stopped position.
        time_slop = call_duration - move_time
        pos_slop = speed*time_slop
        self.assertFalse(actuator.moving)
        self.assertEqual(actuator.remaining_time, 0)
        # Check that stop sets the end position to the current position,
        # but does not change the start position.
        self.assertEqual(actuator.current_position, actuator.end_position)
        self.assertEqual(actuator.start_position, pos)
        pos_error = abs(actuator.end_position - predicted_end_position)
        self.assertLessEqual(pos_error, pos_slop)

        # Stopping a stopped actuator should have no effect.
        stopped_pos = actuator.current_position
        actuator.stop()
        self.assertEqual(actuator.start_position, pos)
        self.assertEqual(actuator.current_position, stopped_pos)
        self.assertEqual(actuator.end_position, stopped_pos)
        self.assertFalse(actuator.moving)
        self.assertEqual(actuator.remaining_time, 0)


if __name__ == '__main__':
    unittest.main()
