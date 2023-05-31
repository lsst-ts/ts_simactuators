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
import typing
import unittest

import pytest
from lsst.ts import simactuators, utils


class TestCircularPointToPointActuator(unittest.IsolatedAsyncioTestCase):
    def test_constructor(self) -> None:
        speed = 1.5
        for start_position in (-360, -180, -1, 0, 359.99, 360):
            tai0 = utils.current_tai()
            actuator = simactuators.CircularPointToPointActuator(
                start_position=start_position,
                speed=speed,
            )
            time_slop = utils.current_tai() - tai0
            if 0 <= start_position < 360:
                assert actuator.start_position == start_position
            else:
                # The reported angle is wrapped
                utils.assert_angles_almost_equal(
                    actuator.start_position, start_position
                )
            assert actuator.end_position == actuator.start_position
            assert actuator.start_tai == pytest.approx(tai0, abs=time_slop)
            assert actuator.start_tai == actuator.end_tai
            assert actuator.direction == 1
            for dt in (-1, 0, 1):
                tai = actuator.start_tai + dt
                utils.assert_angles_almost_equal(actuator.position(tai), start_position)
                assert actuator.velocity(tai) == 0
                assert not actuator.moving(tai)

        start_position = 3
        for bad_speed in (0, -0.001):
            with self.subTest(bad_speed=bad_speed):
                with pytest.raises(ValueError):
                    simactuators.CircularPointToPointActuator(
                        start_position=start_position,
                        speed=bad_speed,
                    )

    async def test_set_position(self) -> None:
        directions = list(simactuators.Direction) + [None]
        tai = utils.current_tai()
        for direction in directions:
            kwargs = dict()
            if direction is not None:
                kwargs["direction"] = direction

            # Move in both directions, crossing the wrap at 0.
            # Test some with a specified start_tai and others without.
            await self.check_set_position(
                start_position=3, end_position=4, start_tai=tai, **kwargs
            )
            await self.check_set_position(
                start_position=4, end_position=3, start_tai=tai, **kwargs
            )
            await self.check_set_position(
                start_position=355, end_position=2, start_tai=tai, **kwargs
            )
            await self.check_set_position(start_position=2, end_position=355, **kwargs)
            await self.check_set_position(start_position=6, end_position=-5, **kwargs)
            await self.check_set_position(start_position=-5, end_position=6, **kwargs)

            # A move to the start position should have no effect.
            pos = 1
            actuator = simactuators.CircularPointToPointActuator(
                start_position=pos, speed=2
            )
            duration = actuator.set_position(pos, **kwargs)
            assert duration == 0
            utils.assert_angles_almost_equal(actuator.start_position, pos)
            utils.assert_angles_almost_equal(actuator.end_position, pos)
            assert actuator.start_tai == actuator.end_tai
            assert not actuator.moving(actuator.start_tai)
            assert actuator.velocity(actuator.start_tai) == 0
            assert actuator.direction == simactuators.Direction.POSITIVE

            # Check specifying an explicit start_tai;
            # pick a value different than the existing start_tai
            # so we can tell the difference.
            start_tai = actuator.start_tai + 5
            actuator.set_position(position=1, start_tai=start_tai, **kwargs)
            assert actuator.start_tai == start_tai

    async def check_set_position(
        self, start_position: float, end_position: float, **kwargs: typing.Any
    ) -> None:
        """Check the set_position command.

        Unlike `CircularPointToPointActuator`, direction=None is allowed
        and means that set_position will not be called with that argument.

        Parameters
        ----------
        start_position : `float`
            Initial position (degrees).
        end_position : `float'
            Target position (degrees).
        kwargs : `dict`
            Additional keyword arguments for set_position;
            may include "direction" and "start_tai".
        """
        if start_position == end_position:
            raise ValueError("start_position must not equal end_position")
        # Make the move a reasonable length
        speed = 2.0 / abs(end_position - start_position)
        actuator = simactuators.CircularPointToPointActuator(
            start_position=start_position,
            speed=speed,
        )
        direction = kwargs.get("direction", simactuators.Direction.NEAREST)
        # Sleep a bit so actuator.start_tai will change by a noticeable amount
        # as a result of the set_position command.
        await asyncio.sleep(0.1)
        # Keep track of how long it takes to call `set_position`
        # so we know how picky to be when testing `start_tai`.
        tai0 = utils.current_tai()
        duration = actuator.set_position(end_position, **kwargs)
        time_slop = utils.current_tai() - tai0
        min_delta_position = utils.angle_diff(end_position, start_position).deg
        delta_position = min_delta_position
        if direction == simactuators.Direction.POSITIVE and min_delta_position < 0:
            delta_position = 360 + min_delta_position
        elif direction == simactuators.Direction.NEGATIVE and min_delta_position > 0:
            delta_position = 360 - min_delta_position
        predicted_duration = abs(delta_position) / speed
        assert duration == pytest.approx(predicted_duration)
        if "start_tai" in kwargs:
            assert actuator.start_tai == kwargs["start_tai"]
        else:
            assert actuator.start_tai == pytest.approx(tai0, abs=time_slop)
        assert actuator.start_tai + predicted_duration == pytest.approx(
            actuator.end_tai, abs=0.000001
        )
        utils.assert_angles_almost_equal(actuator.end_position, end_position)
        if direction == simactuators.Direction.NEAREST:
            predicted_direction = (
                simactuators.Direction.POSITIVE
                if min_delta_position >= 0
                else simactuators.Direction.NEGATIVE
            )
        else:
            predicted_direction = direction
        assert actuator.direction == predicted_direction
        predicted_speed = speed * predicted_direction

        # The actuator should not be moving before or after the move.
        for tai in (
            actuator.start_tai - 1,
            actuator.start_tai,
            actuator.end_tai,
            actuator.end_tai + 1,
        ):
            assert not actuator.moving(tai)
            assert actuator.velocity(tai) == 0

        # The actuator should be moving during the move.
        for tai in (
            actuator.start_tai + 0.001,
            (actuator.start_tai + actuator.end_tai) / 2,
            actuator.end_tai - 0.001,
        ):
            assert actuator.moving(tai)
            assert actuator.velocity(tai) == predicted_speed

    async def test_stop(self) -> None:
        # Pick start and end positions that will not wrap
        # so we can use assertEqual where appropriate.
        start_position = 3
        target_position = 4
        # Slow motion so plenty of time to stop
        speed = 0.1
        actuator = simactuators.CircularPointToPointActuator(
            start_position=start_position,
            speed=speed,
        )
        duration = actuator.set_position(
            target_position, direction=simactuators.Direction.NEAREST
        )
        assert actuator.end_position == target_position
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
        utils.assert_angles_almost_equal(
            actuator.end_position, actuator.position(tai0), max_diff=position_slop
        )

        move_time = actuator.end_tai - actuator.start_tai
        assert move_time > 0
        predicted_end_position = start_position + speed * move_time
        assert predicted_end_position == pytest.approx(
            actuator.end_position, abs=0.000001
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


if __name__ == "__main__":
    unittest.main()
