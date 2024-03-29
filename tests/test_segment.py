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
import sys
import unittest

import numpy as np
import pytest
from lsst.ts import simactuators


class TestPathSegment(unittest.TestCase):
    def check_from_end_conditions(
        self,
        start_position: float,
        start_velocity: float,
        end_position: float,
        end_velocity: float,
        dt: float,
    ) -> None:
        """Check PathSegment from_end_conditions and limits methods."""
        start_tai = 45.3  # any value will do
        end_tai = start_tai + dt
        segment = simactuators.path.PathSegment.from_end_conditions(
            start_tai=start_tai,
            start_position=start_position,
            start_velocity=start_velocity,
            end_tai=end_tai,
            end_position=end_position,
            end_velocity=end_velocity,
        )
        assert segment.tai == pytest.approx(start_tai)
        assert segment.position == pytest.approx(start_position)
        assert segment.velocity == pytest.approx(start_velocity)

        end_segment = segment.at(end_tai)
        assert end_segment.position == pytest.approx(end_position)
        assert end_segment.velocity == pytest.approx(end_velocity)
        start_accel = segment.acceleration
        jerk = segment.jerk
        desired_end_position = start_position + dt * (
            start_velocity + dt * (0.5 * start_accel + dt * (1 / 6) * jerk)
        )
        desired_end_velocity = start_velocity + dt * (start_accel + dt * 0.5 * jerk)
        assert end_segment.position == pytest.approx(desired_end_position)
        assert end_segment.velocity == pytest.approx(desired_end_velocity)

        # estimate position and velocity limits by computing at many points
        desired_min_position = min(start_position, end_position)
        desired_max_position = max(start_position, end_position)
        desired_max_velocity = max(abs(start_velocity), abs(end_velocity))
        for t in np.linspace(start=0, stop=dt, num=100):
            pt = start_position + t * (
                start_velocity + t * (0.5 * start_accel + t * jerk / 6)
            )
            vt = start_velocity + t * (start_accel + t * 0.5 * jerk)
            desired_min_position = min(pt, desired_min_position)
            desired_max_position = max(pt, desired_max_position)
            desired_max_velocity = max(abs(vt), desired_max_velocity)
        aB = start_accel + dt * jerk
        desired_max_acceleration = max(abs(start_accel), abs(aB))

        limits = segment.limits(end_tai)
        assert limits.min_position == pytest.approx(desired_min_position, abs=0.001)
        assert limits.max_position == pytest.approx(desired_max_position, abs=0.001)
        assert limits.max_velocity == pytest.approx(desired_max_velocity, abs=0.001)
        assert limits.max_acceleration == pytest.approx(desired_max_acceleration)

    def test_from_end_conditions(self) -> None:
        for (
            start_position,
            start_velocity,
            end_position,
            end_velocity,
            dt,
        ) in itertools.product(
            (0, -0.5, 0.2),
            (0, -0.2, 0.1),
            (0, 0.3, -0.6),
            (0, 0.3, -0.2),
            (1, 5),
        ):
            with self.subTest(
                start_position=start_position,
                start_velocity=start_velocity,
                end_position=end_position,
                end_velocity=end_velocity,
                dt=dt,
            ):
                self.check_from_end_conditions(
                    dt=dt,
                    start_position=start_position,
                    start_velocity=start_velocity,
                    end_position=end_position,
                    end_velocity=end_velocity,
                )

    def test_basics(self) -> None:
        tai = 1573847242.4
        position = 5.1
        velocity = 4.3
        acceleration = 0.23
        jerk = 0.05
        segment = simactuators.path.PathSegment(
            tai=tai,
            position=position,
            velocity=velocity,
            acceleration=acceleration,
            jerk=jerk,
        )
        assert segment.tai == tai
        assert segment.position == position
        assert segment.velocity == velocity
        assert segment.acceleration == acceleration
        assert segment.jerk == jerk

        # The at command at the same time should return a copy
        copied_segment = segment.at(tai)
        assert copied_segment.tai == tai
        assert copied_segment.position == pytest.approx(position)
        assert copied_segment.velocity == pytest.approx(velocity)
        assert copied_segment.acceleration == pytest.approx(acceleration)
        assert copied_segment.jerk == jerk

        for dt in (-5.1, -3.23, 1.23, 6.6):
            with self.subTest(dt=dt):
                new_tai = tai + dt
                new_segment = segment.at(new_tai)
                # Note: I am using slightly cruder math here than in
                # `PathSegment` (to avoid an exact copy),
                # so the computed values will be slightly different,
                # especially position and to a lessar extent velocity.
                expected_pos = (
                    position
                    + velocity * dt
                    + acceleration * dt * dt / 2
                    + jerk * dt * dt * dt / 6
                )
                expected_vel = velocity + acceleration * dt + jerk * dt * dt / 2
                expected_accel = acceleration + jerk * dt
                assert new_segment.tai == new_tai
                assert new_segment.position == pytest.approx(expected_pos, abs=0.00001)
                assert new_segment.velocity == pytest.approx(
                    expected_vel, abs=0.0000001
                )
                assert new_segment.acceleration == pytest.approx(expected_accel)
                assert new_segment.jerk == jerk

    def test_default_arguments(self) -> None:
        """Test constructing a PathSegment with default arguments."""
        full_kwargs = dict(
            tai=1573847242.4,
            position=5.1,
            velocity=4.3,
            acceleration=0.23,
            jerk=0.05,
        )
        for field_to_omit in full_kwargs:
            if field_to_omit == "tai":
                continue
            with self.subTest(field_to_omit=field_to_omit):
                kwargs = full_kwargs.copy()
                del kwargs[field_to_omit]
                segment = simactuators.path.PathSegment(**kwargs)
                for field_name in full_kwargs:
                    if field_name == field_to_omit:
                        assert getattr(segment, field_name) == 0
                    else:
                        assert getattr(segment, field_name) == full_kwargs[field_name]

    def test_invalid_inputs(self) -> None:
        """Test invalid inputs for PathSegment.from_end_conditions.

        The only thing that can go wrong is end_tai - start_tai <= 0.
        """
        for dt in (0, math.sqrt(sys.float_info.min)):
            with self.subTest(dt=dt):
                start_tai = 51.0  # arbitrary
                end_tai = start_tai + dt

                with pytest.raises(ValueError):
                    simactuators.path.PathSegment.from_end_conditions(
                        start_tai,
                        start_position=1,
                        start_velocity=2,
                        end_tai=end_tai,
                        end_position=1,
                        end_velocity=2,
                    )
