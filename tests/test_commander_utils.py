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
import math
import unittest

import numpy as np
from lsst.ts import simactuators, utils


class TestCommanderUtils(unittest.IsolatedAsyncioTestCase):
    def test_cosine_generator_constructor_errors(self) -> None:
        center_positions = [1, -2, 3.3]
        amplitudes = [-2.1, 0, 1.2]
        max_speeds = [1, 0, -1]  # Check that sign is ignored
        simactuators.CosineGenerator(
            center_positions=center_positions,
            amplitudes=amplitudes,
            max_speeds=max_speeds,
        )

        with self.assertRaises(ValueError):
            simactuators.CosineGenerator(
                center_positions=center_positions[0:2],  # Too short
                amplitudes=amplitudes,
                max_speeds=max_speeds,
            )
        with self.assertRaises(ValueError):
            simactuators.CosineGenerator(
                center_positions=center_positions,
                amplitudes=amplitudes[0:2],  # Too short
                max_speeds=max_speeds,
            )
        with self.assertRaises(ValueError):
            simactuators.CosineGenerator(
                center_positions=center_positions,
                amplitudes=amplitudes,
                max_speeds=max_speeds[0:2],  # Too short
            )
        with self.assertRaises(ValueError):
            simactuators.CosineGenerator(
                center_positions=center_positions,
                amplitudes=amplitudes,
                max_speeds=[0, 0, 0],  # Zero speed for nonzero amplitude
            )

    def test_ramp_generator_constructor_errors(self) -> None:
        start_positions = [1, -2, 3.3]
        end_positions = [-3, -2, 4.1]
        speeds = [1, 0, -1]  # Check that sign is ignored
        simactuators.RampGenerator(
            start_positions=start_positions,
            end_positions=end_positions,
            speeds=speeds,
        )

        with self.assertRaises(ValueError):
            simactuators.RampGenerator(
                start_positions=start_positions[0:2],  # Too short
                end_positions=end_positions,
                speeds=speeds,
            )
        with self.assertRaises(ValueError):
            simactuators.RampGenerator(
                start_positions=start_positions,
                end_positions=end_positions[0:2],  # Too short
                speeds=speeds,
            )
        with self.assertRaises(ValueError):
            simactuators.RampGenerator(
                start_positions=start_positions,
                end_positions=end_positions,
                speeds=speeds[0:2],  # Too short
            )
        with self.assertRaises(ValueError):
            simactuators.RampGenerator(
                start_positions=start_positions,
                end_positions=end_positions,
                speeds=[0, 0, 0],  # Zero speed for nonzero delta
            )

    async def test_cosine_generator_path(self) -> None:
        center_positions = [1, -2, 3.3]
        amplitudes = [-1, 0, 1.1]
        max_speeds = [-2.2, 0, 2.1]  # Use the correct sign to simplify testing
        advance_time = 0.2
        nextra = 4

        start_and_end_position_arr = np.add(center_positions, amplitudes)

        predicted_durations = []
        for i in range(3):
            if amplitudes[i] == 0:
                duration = 0.0
            else:
                duration = abs(amplitudes[i] * 2 * math.pi / max_speeds[i])
            predicted_durations.append(duration)

        cosine_generator = simactuators.CosineGenerator(
            center_positions=center_positions,
            amplitudes=amplitudes,
            max_speeds=max_speeds,
            advance_time=advance_time,
            nextra=nextra,
        )

        predicted_duration = max(*predicted_durations)
        self.assertAlmostEqual(cosine_generator.duration, predicted_duration)

        positions_list = []
        velocities_list = []
        tais = []
        for positions, velocities, tai in cosine_generator():
            positions_list.append(positions)
            velocities_list.append(velocities)
            tais.append(tai)
            curr_tai = utils.current_tai()
            self.assertAlmostEqual(curr_tai + advance_time, tai, delta=0.1)
            await asyncio.sleep(0.02)

        start_tai = tais[0]
        np.testing.assert_allclose(positions_list[0], start_and_end_position_arr)
        previous_position_arr = None
        previous_velocity_arr = None
        previous_tai = None
        # Use simplistic extrapolation to predict the next position.
        # This catches gross scale errors.
        # Also record the extreme positions and velocities.
        min_position_arr = np.array(center_positions)
        max_position_arr = np.array(center_positions)
        min_velocity_arr = np.zeros(3)
        max_velocity_arr = np.zeros(3)
        for position_arr, velocity_arr, tai in zip(
            positions_list, velocities_list, tais
        ):
            min_position_arr = np.minimum(min_position_arr, position_arr)
            max_position_arr = np.maximum(max_position_arr, position_arr)
            min_velocity_arr = np.minimum(min_velocity_arr, velocity_arr)
            max_velocity_arr = np.maximum(max_velocity_arr, velocity_arr)
            dt_start = tai - start_tai
            if previous_tai is not None:
                dt_previous = tai - previous_tai
                predicted_position_arr = (
                    previous_position_arr + previous_velocity_arr * dt_previous
                )
            for i in range(len(center_positions)):
                if dt_start < predicted_durations[i]:
                    if previous_tai is not None:
                        self.assertAlmostEqual(
                            position_arr[i], predicted_position_arr[i], delta=0.01
                        )
                else:
                    self.assertAlmostEqual(
                        position_arr[i], start_and_end_position_arr[i]
                    )
                    self.assertEqual(velocity_arr[i], 0)
            previous_position_arr = position_arr
            previous_velocity_arr = velocity_arr
            previous_tai = tai

        np.testing.assert_allclose(
            max_position_arr, np.add(center_positions, np.abs(amplitudes)), atol=0.01
        )
        np.testing.assert_allclose(
            min_position_arr,
            np.subtract(center_positions, np.abs(amplitudes)),
            atol=0.01,
        )
        np.testing.assert_allclose(max_velocity_arr, np.abs(max_speeds), atol=0.01)
        np.testing.assert_allclose(min_velocity_arr, -np.abs(max_speeds), atol=0.01)

        # Check nextra
        for position_arr, velocity_arr, tai in zip(
            positions_list[-nextra:], velocities_list[-nextra:], tais[-nextra:]
        ):
            np.testing.assert_allclose(position_arr, start_and_end_position_arr)
            np.testing.assert_allclose(velocity_arr, np.zeros(3))

        self.assertFalse(np.all(positions_list[-nextra - 1] == 0))
        self.assertFalse(np.all(velocities_list[-nextra - 1] == 0))

    async def test_ramp_generator_path(self) -> None:
        start_positions = [1, -2, 3.3]
        end_positions = [-3, -2, 5.1]
        speeds = [-2.2, 0, 1]  # Use the correct sign to simplify testing
        advance_time = 0.2
        nextra = 4

        predicted_durations = []
        for i in range(3):
            dpos = end_positions[i] - start_positions[i]
            if dpos == 0:
                duration = 0.0
            else:
                duration = abs(dpos / speeds[i])
            predicted_durations.append(duration)

        ramp_generator = simactuators.RampGenerator(
            start_positions=start_positions,
            end_positions=end_positions,
            speeds=speeds,
            advance_time=advance_time,
            nextra=nextra,
        )

        predicted_duration = max(*predicted_durations)
        self.assertAlmostEqual(ramp_generator.duration, predicted_duration)

        positions_list = []
        velocities_list = []
        tais = []
        for positions, velocities, tai in ramp_generator():
            positions_list.append(positions)
            velocities_list.append(velocities)
            tais.append(tai)
            curr_tai = utils.current_tai()
            self.assertAlmostEqual(curr_tai + advance_time, tai, delta=0.1)
            await asyncio.sleep(0.1)

        np.testing.assert_allclose(positions_list[0], start_positions)
        start_tai = tais[0]
        for position_arr, velocity_arr, tai in zip(
            positions_list, velocities_list, tais
        ):
            dt = tai - start_tai
            for i in range(len(start_positions)):
                if dt < predicted_durations[i]:
                    self.assertAlmostEqual(
                        position_arr[i], start_positions[i] + speeds[i] * dt
                    )
                    self.assertAlmostEqual(velocity_arr[i], speeds[i])
                else:
                    self.assertAlmostEqual(position_arr[i], end_positions[i])
                    self.assertEqual(velocity_arr[i], 0)

        # Check nextra
        for position_arr, velocity_arr, tai in zip(
            positions_list[-nextra:], velocities_list[-nextra:], tais[-nextra:]
        ):
            np.testing.assert_allclose(position_arr, end_positions)
            np.testing.assert_allclose(velocity_arr, np.zeros(3))

        self.assertFalse(np.all(positions_list[-nextra - 1] == 0))
        self.assertFalse(np.all(velocities_list[-nextra - 1] == 0))
