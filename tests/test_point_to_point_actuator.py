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

import pytest
from lsst.ts import simactuators, utils


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
                assert actuator.start_tai == pytest.approx(tai0, abs=time_slop)
                assert actuator.start_position == good_start_position
                assert actuator.end_position == actuator.start_position
                assert actuator.start_tai == actuator.end_tai
                assert actuator.direction == simactuators.Direction.POSITIVE
                for dt in (-1, 0, 1):
                    tai = actuator.start_tai + dt
                    assert actuator.position(tai) == good_start_position
                    assert actuator.velocity(tai) == 0
                    assert not actuator.moving(tai)
                    predicted_remaining_time = 0 if dt >= 0 else -dt
                    assert actuator.remaining_time(tai) == pytest.approx(
                        predicted_remaining_time
                    )

                assert actuator.position() == good_start_position
                assert actuator.velocity() == 0
                assert not actuator.moving()
                assert actuator.remaining_time() == 0

        for bad_min_position in (max_position, max_position + 0.001):
            with self.subTest(bad_min_position=bad_min_position):
                with pytest.raises(ValueError):
                    simactuators.PointToPointActuator(
                        min_position=bad_min_position,
                        max_position=max_position,
                        start_position=max_position,
                        speed=speed,
                    )

        for bad_max_position in (min_position, min_position - 0.001):
            with self.subTest(bad_max_position=bad_max_position):
                with pytest.raises(ValueError):
                    simactuators.PointToPointActuator(
                        min_position=min_position,
                        max_position=bad_max_position,
                        start_position=min_position,
                        speed=speed,
                    )

        for bad_pos in (min_position - 0.001, max_position + 0.001):
            with self.subTest(bad_pos=bad_pos):
                with pytest.raises(ValueError):
                    simactuators.PointToPointActuator(
                        min_position=min_position,
                        max_position=max_position,
                        start_position=bad_pos,
                        speed=speed,
                    )

        for bad_speed in (0, -0.001):
            with self.subTest(bad_speed=bad_speed):
                with pytest.raises(ValueError):
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
        assert duration == 0
        assert actuator.start_position == pos
        assert actuator.end_position == pos
        assert actuator.start_tai == actuator.end_tai
        assert not actuator.moving(actuator.start_tai)
        assert actuator.velocity(actuator.start_tai) == 0
        assert actuator.direction == simactuators.Direction.POSITIVE

        # Check specifying an explicit start_tai;
        # pick a value different than the existing start_tai
        # so we can tell the difference.
        start_tai = actuator.start_tai + 5
        actuator.set_position(position=1, start_tai=start_tai)
        assert actuator.start_tai == start_tai

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
        assert duration == pytest.approx(predicted_duration)
        assert actuator.start_tai == pytest.approx(tai0, abs=time_slop)
        assert remaining_time == pytest.approx(duration, abs=time_slop)
        assert actuator.start_tai + predicted_duration == pytest.approx(
            actuator.end_tai, abs=0.000001
        )
        assert actuator.end_position == end_position
        predicted_direction = (
            simactuators.Direction.POSITIVE
            if end_position >= start_position
            else simactuators.Direction.NEGATIVE
        )
        assert actuator.direction == predicted_direction
        predicted_speed = speed * predicted_direction

        # The actuator should not be moving before or after the move
        # (but remaining_time is > 0 before the move).
        for tai in (
            actuator.start_tai - 1,
            actuator.start_tai,
            actuator.end_tai,
            actuator.end_tai + 1,
        ):
            assert not actuator.moving(tai)
            assert actuator.velocity(tai) == 0
            predicted_remaining_time = (
                actuator.end_tai - tai if actuator.end_tai > tai else 0
            )
            assert actuator.remaining_time(tai) == pytest.approx(
                predicted_remaining_time
            )

        # The actuator should be moving during the move.
        for tai in (
            actuator.start_tai + 0.001,
            (actuator.start_tai + actuator.end_tai) / 2,
            actuator.end_tai - 0.001,
        ):
            assert actuator.moving(tai)
            assert actuator.velocity(tai) == predicted_speed
            predicted_remaining_time = actuator.end_tai - tai
            assert actuator.remaining_time(tai) == pytest.approx(
                predicted_remaining_time
            )

        # Try moving out of bounds.
        for bad_end_position in (
            min_position - 1,
            min_position - 0.0001,
            max_position + 0.0001,
            max_position + 1,
        ):
            with pytest.raises(ValueError):
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
        assert actuator.end_position == end_position
        assert actuator.moving(utils.current_tai())
        # Let the actuator move for some arbitrary time
        # that is less than the remaining time
        sleep_time = 0.21
        assert duration > sleep_time
        await asyncio.sleep(sleep_time)
        tai0 = utils.current_tai()
        actuator.stop()
        time_slop = utils.current_tai() - tai0
        assert tai0 == pytest.approx(actuator.end_tai, abs=time_slop)
        position_slop = time_slop * speed
        assert actuator.end_position == pytest.approx(
            actuator.position(tai0), abs=position_slop
        )

        move_time = actuator.end_tai - actuator.start_tai
        assert move_time > 0
        predicted_end_position = start_position + speed * move_time
        assert predicted_end_position == pytest.approx(
            actuator.end_position, abs=0.0000001
        )
        assert not actuator.moving(actuator.end_tai)
        assert actuator.position(actuator.end_tai) == actuator.end_position
        assert actuator.start_position == start_position
        assert actuator.position(actuator.start_tai) == actuator.start_position

        # Stopping a stopped actuator should have no effect
        # except updating end_tai. Sleep long enough
        # that we can see the end_tai change.
        await asyncio.sleep(0.1)
        old_start_tai = actuator.start_tai
        old_pos = actuator.end_position
        tai0 = utils.current_tai()
        actuator.stop()
        time_slop = utils.current_tai()
        assert actuator.start_tai == old_start_tai
        assert actuator.start_position == start_position
        assert actuator.end_position == old_pos
        assert actuator.end_tai == pytest.approx(tai0, abs=time_slop)

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
        assert actuator.end_position == end_position
        predicted_duration = (end_position - start_position) / speed
        assert duration == pytest.approx(predicted_duration)
        end_tai = start_tai + duration
        assert actuator.end_tai == pytest.approx(end_tai)
        assert not actuator.moving(tai=start_tai - 0.001)
        assert not actuator.moving(tai=end_tai + 0.001)
        assert actuator.moving(tai=start_tai + 0.001)
        assert actuator.moving(tai=end_tai - 0.001)

        # Let the actuator move for some arbitrary time
        # that is less than the remaining time
        stop_dtime = 0.21
        stop_tai = stop_dtime + start_tai
        assert duration > stop_dtime
        actuator.stop(tai=stop_tai)

        assert actuator.end_tai == stop_tai
        assert not actuator.moving(tai=start_tai - 0.001)
        assert not actuator.moving(tai=stop_tai + 0.001)
        assert actuator.moving(tai=start_tai + 0.001)
        assert actuator.moving(tai=stop_tai - 0.001)
        assert actuator.end_tai == stop_tai
        predicted_end_position = start_position + speed * stop_dtime
        assert actuator.end_position == pytest.approx(predicted_end_position)
