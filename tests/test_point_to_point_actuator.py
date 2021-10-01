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

import asyncio
import unittest

from lsst.ts import utils
from lsst.ts import simactuators


class TestPointToPointActuator(unittest.IsolatedAsyncioTestCase):
    def test_constructor(self) -> None:
        min_position = -1
        max_position = 5
        start_position = 3
        speed = 1.5
        for good_start_position in (start_position, min_position, max_position):
            with self.subTest(good_start_position=good_start_position):
                tai0 = utils.current_tai()
                actuator = simactuators.PointToPointActuator(
                    min_position=min_position,
                    max_position=max_position,
                    start_position=good_start_position,
                    speed=speed,
                )
                time_slop = utils.current_tai() - tai0
                self.assertAlmostEqual(actuator.start_tai, tai0, delta=time_slop)
                self.assertEqual(actuator.start_position, good_start_position)
                self.assertEqual(actuator.end_position, actuator.start_position)
                self.assertEqual(actuator.start_tai, actuator.end_tai)
                self.assertEqual(actuator.direction, simactuators.Direction.POSITIVE)
                for dt in (-1, 0, 1):
                    tai = actuator.start_tai + dt
                    self.assertEqual(actuator.position(tai), good_start_position)
                    self.assertEqual(actuator.velocity(tai), 0)
                    self.assertFalse(actuator.moving(tai))
                    predicted_remaining_time = 0 if dt >= 0 else -dt
                    self.assertAlmostEqual(
                        actuator.remaining_time(tai), predicted_remaining_time
                    )

                self.assertEqual(actuator.position(), good_start_position)
                self.assertEqual(actuator.velocity(), 0)
                self.assertFalse(actuator.moving())
                self.assertEqual(actuator.remaining_time(), 0)

        for bad_min_position in (max_position, max_position + 0.001):
            with self.subTest(bad_min_position=bad_min_position):
                with self.assertRaises(ValueError):
                    simactuators.PointToPointActuator(
                        min_position=bad_min_position,
                        max_position=max_position,
                        start_position=max_position,
                        speed=speed,
                    )

        for bad_max_position in (min_position, min_position - 0.001):
            with self.subTest(bad_max_position=bad_max_position):
                with self.assertRaises(ValueError):
                    simactuators.PointToPointActuator(
                        min_position=min_position,
                        max_position=bad_max_position,
                        start_position=min_position,
                        speed=speed,
                    )

        for bad_pos in (min_position - 0.001, max_position + 0.001):
            with self.subTest(bad_pos=bad_pos):
                with self.assertRaises(ValueError):
                    simactuators.PointToPointActuator(
                        min_position=min_position,
                        max_position=max_position,
                        start_position=bad_pos,
                        speed=speed,
                    )

        for bad_speed in (0, -0.001):
            with self.subTest(bad_speed=bad_speed):
                with self.assertRaises(ValueError):
                    simactuators.PointToPointActuator(
                        min_position=min_position,
                        max_position=max_position,
                        start_position=start_position,
                        speed=bad_speed,
                    )

    async def test_set_position(self) -> None:
        await self.check_set_position(start_position=3, end_position=4)
        await self.check_set_position(start_position=2, end_position=-5)

        # A move to the start position should have no effect.
        pos = 1
        actuator = simactuators.PointToPointActuator(
            min_position=pos - 1,
            max_position=pos + 1,
            start_position=pos,
            speed=2,
        )
        duration = actuator.set_position(pos)
        self.assertEqual(duration, 0)
        self.assertEqual(actuator.start_position, pos)
        self.assertEqual(actuator.end_position, pos)
        self.assertEqual(actuator.start_tai, actuator.end_tai)
        self.assertFalse(actuator.moving(actuator.start_tai))
        self.assertEqual(actuator.velocity(actuator.start_tai), 0)
        self.assertEqual(actuator.direction, simactuators.Direction.POSITIVE)

        # Check specifying an explicit start_tai;
        # pick a value different than the existing start_tai
        # so we can tell the difference.
        start_tai = actuator.start_tai + 5
        actuator.set_position(position=1, start_tai=start_tai)
        self.assertEqual(actuator.start_tai, start_tai)

    async def check_set_position(
        self, start_position: float, end_position: float
    ) -> None:
        if start_position == end_position:
            raise ValueError("start_position must not equal end_position")
        min_position = min(start_position, end_position) - 1
        max_position = max(start_position, end_position) + 1
        # Make the move take a reasonable amount of time
        speed = 2.0 / abs(end_position - start_position)
        actuator = simactuators.PointToPointActuator(
            min_position=min_position,
            max_position=max_position,
            start_position=start_position,
            speed=speed,
        )
        # Sleep a bit so actuator.start_tai will change by a noticeable amount
        # as a result of the set_position command.
        await asyncio.sleep(0.1)
        # Keep track of how long it takes to call `set_position`
        # so we know how picky to be when testing `start_tai`.
        tai0 = utils.current_tai()
        duration = actuator.set_position(end_position)
        remaining_time = actuator.remaining_time()
        time_slop = utils.current_tai() - tai0
        predicted_duration = abs(end_position - start_position) / speed
        self.assertAlmostEqual(duration, predicted_duration)
        self.assertAlmostEqual(actuator.start_tai, tai0, delta=time_slop)
        self.assertAlmostEqual(remaining_time, duration, delta=time_slop)
        self.assertAlmostEqual(
            actuator.start_tai + predicted_duration, actuator.end_tai, places=6
        )
        self.assertEqual(actuator.end_position, end_position)
        predicted_direction = (
            simactuators.Direction.POSITIVE
            if end_position >= start_position
            else simactuators.Direction.NEGATIVE
        )
        self.assertEqual(actuator.direction, predicted_direction)
        predicted_speed = speed * predicted_direction

        # The actuator should not be moving before or after the move
        # (but remaining_time is > 0 before the move).
        for tai in (
            actuator.start_tai - 1,
            actuator.start_tai,
            actuator.end_tai,
            actuator.end_tai + 1,
        ):
            self.assertFalse(actuator.moving(tai))
            self.assertEqual(actuator.velocity(tai), 0)
            predicted_remaining_time = (
                actuator.end_tai - tai if actuator.end_tai > tai else 0
            )
            self.assertAlmostEqual(
                actuator.remaining_time(tai), predicted_remaining_time
            )

        # The actuator should be moving during the move.
        for tai in (
            actuator.start_tai + 0.001,
            (actuator.start_tai + actuator.end_tai) / 2,
            actuator.end_tai - 0.001,
        ):
            self.assertTrue(actuator.moving(tai))
            self.assertEqual(actuator.velocity(tai), predicted_speed)
            predicted_remaining_time = actuator.end_tai - tai
            self.assertAlmostEqual(
                actuator.remaining_time(tai), predicted_remaining_time
            )

        # Try moving out of bounds.
        for bad_end_position in (
            min_position - 1,
            min_position - 0.0001,
            max_position + 0.0001,
            max_position + 1,
        ):
            with self.assertRaises(ValueError):
                actuator.set_position(bad_end_position)

    async def test_stop_default_tai(self) -> None:
        min_position = -10
        max_position = 10
        start_position = 3
        # Slow motion so plenty of time to stop
        speed = 0.1
        actuator = simactuators.PointToPointActuator(
            min_position=min_position,
            max_position=max_position,
            start_position=start_position,
            speed=speed,
        )
        end_position = 4
        duration = actuator.set_position(end_position)
        self.assertEqual(actuator.end_position, end_position)
        self.assertTrue(actuator.moving(utils.current_tai()))
        # Let the actuator move for some arbitrary time
        # that is less than the remaining time
        sleep_time = 0.21
        self.assertGreater(duration, sleep_time)
        await asyncio.sleep(sleep_time)
        tai0 = utils.current_tai()
        actuator.stop()
        time_slop = utils.current_tai() - tai0
        self.assertAlmostEqual(tai0, actuator.end_tai, delta=time_slop)
        position_slop = time_slop * speed
        self.assertAlmostEqual(
            actuator.end_position, actuator.position(tai0), delta=position_slop
        )

        move_time = actuator.end_tai - actuator.start_tai
        self.assertGreater(move_time, 0)
        predicted_end_position = start_position + speed * move_time
        self.assertAlmostEqual(predicted_end_position, actuator.end_position, places=6)
        self.assertFalse(actuator.moving(actuator.end_tai))
        self.assertEqual(actuator.position(actuator.end_tai), actuator.end_position)
        self.assertEqual(actuator.start_position, start_position)
        self.assertEqual(actuator.position(actuator.start_tai), actuator.start_position)

        # Stopping a stopped actuator should have no effect
        # except updating end_tai. Sleep long enough
        # that we can see the end_tai change.
        await asyncio.sleep(0.1)
        old_start_tai = actuator.start_tai
        old_pos = actuator.end_position
        tai0 = utils.current_tai()
        actuator.stop()
        time_slop = utils.current_tai()
        self.assertEqual(actuator.start_tai, old_start_tai)
        self.assertEqual(actuator.start_position, start_position)
        self.assertEqual(actuator.end_position, old_pos)
        self.assertAlmostEqual(actuator.end_tai, tai0, delta=time_slop)

    async def test_stop_specified_tai(self) -> None:
        min_position = -10
        max_position = 10
        start_position = 3
        # Slow motion so plenty of time to stop
        speed = 0.1
        actuator = simactuators.PointToPointActuator(
            min_position=min_position,
            max_position=max_position,
            start_position=start_position,
            speed=speed,
        )
        end_position = 4
        start_tai = 12.1  # arbitrary
        duration = actuator.set_position(end_position, start_tai=start_tai)
        self.assertEqual(actuator.end_position, end_position)
        predicted_duration = (end_position - start_position) / speed
        self.assertAlmostEqual(duration, predicted_duration)
        end_tai = start_tai + duration
        self.assertAlmostEqual(actuator.end_tai, end_tai)
        self.assertFalse(actuator.moving(tai=start_tai - 0.001))
        self.assertFalse(actuator.moving(tai=end_tai + 0.001))
        self.assertTrue(actuator.moving(tai=start_tai + 0.001))
        self.assertTrue(actuator.moving(tai=end_tai - 0.001))

        # Let the actuator move for some arbitrary time
        # that is less than the remaining time
        stop_dtime = 0.21
        stop_tai = stop_dtime + start_tai
        self.assertGreater(duration, stop_dtime)
        actuator.stop(tai=stop_tai)

        self.assertEqual(actuator.end_tai, stop_tai)
        self.assertFalse(actuator.moving(tai=start_tai - 0.001))
        self.assertFalse(actuator.moving(tai=stop_tai + 0.001))
        self.assertTrue(actuator.moving(tai=start_tai + 0.001))
        self.assertTrue(actuator.moving(tai=stop_tai - 0.001))
        self.assertEqual(actuator.end_tai, stop_tai)
        predicted_end_position = start_position + speed * stop_dtime
        self.assertAlmostEqual(actuator.end_position, predicted_end_position)
