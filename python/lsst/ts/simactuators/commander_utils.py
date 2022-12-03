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

__all__ = ["CosineGenerator", "RampGenerator"]

import typing

import numpy as np
from lsst.ts import utils


class CosineGenerator:
    def __init__(
        self,
        center_positions: typing.Sequence[float],
        amplitudes: typing.Sequence[float],
        max_speeds: typing.Sequence[float],
        advance_time: float = 0.1,
        nextra: int = 5,
    ) -> None:
        """Functor to generate position, velocity and time for one period
        of a cosine, in one or more axes.

        This is designed to be used in a CSC commander
        to provide a cosine tracking target path.

        Parameters
        ----------
        center_positions : `List` [`float`]
            Center position for each axis (deg).
        amplitudes : `List` [`float`]
            Peak to center amplitude of motion (deg). Sign matters.
            Motion starts and ends at amplitude + center_pos.
        max_speeds : `List` [`float`]
            Maximum speed for each axis (deg/sec). Sign is ignored.
            Values are ignored for any axis for which amplitude = 0
        advance_time : `float`
            Advance time for returned TAI timestamps (seconds).
            Each returned TAI = current TAI + advance_time.
        nextra : `int`
            Number of extra return values after all axes are done moving.
            These extra values have
            position = center_positions + amplitudes
            and velocity=[0, 0, ..., 0]

        Raises
        ------
        ValueError
            If len(center_positions) != len(amplitudes),
            len(center_positions) == len(max_speeds),
            amplitude is 0 for all axes,
            or max_speed is zero and amplitude is not for a given axis.
        """
        self.naxes = len(center_positions)
        if self.naxes != len(amplitudes):
            raise ValueError(
                f"center_positions={center_positions} and amplitudes={amplitudes} "
                "must have the same length"
            )
        if self.naxes != len(max_speeds):
            raise ValueError(
                f"center_positions={center_positions} and max_speeds={max_speeds} "
                "must have the same length"
            )
        for i in range(self.naxes):
            if max_speeds[i] == 0 and amplitudes[i] != 0:
                raise ValueError(
                    f"max_speeds={max_speeds} zero and amplitudes={amplitudes} nonzero "
                    "for one or more axes"
                )
        self.center_position_arr = np.array(center_positions, dtype=float)
        self.amplitude_arr = np.array(amplitudes, dtype=float)
        self.end_position_arr = self.center_position_arr + self.amplitude_arr
        self.max_velocity_arr = np.copysign(
            np.array(max_speeds, dtype=float), self.amplitude_arr
        )
        if np.all(self.amplitude_arr == 0):
            raise ValueError("amplitude must be nonzero for at least one axis")
        self.advance_time = float(advance_time)
        self.nextra = int(nextra)

        # Avoid needless divide-by-zero warnings.
        # Set safe_max_velocity_divisor_arr to an arbitrary nonzero value
        # for those axes with no motion (self.amplitude_arr == 0).
        safe_max_velocity_divisor_arr = np.where(
            self.amplitude_arr == 0, 1, self.max_velocity_arr
        )
        # Period; 0 for axes with no amplitude
        self.period_arr = (
            np.abs(self.amplitude_arr)
            * 2
            * np.pi
            / np.abs(safe_max_velocity_divisor_arr)
        )

        safe_period_divisor_arr = np.where(self.amplitude_arr == 0, 1, self.period_arr)
        # Inverse period; 0 for axes with no amplitude
        self.inverse_period_arr = np.where(
            self.amplitude_arr == 0, 1, 1 / safe_period_divisor_arr
        )

    @property
    def duration(self) -> float:
        return np.amax(self.period_arr)

    def __call__(
        self,
    ) -> typing.Generator[typing.Tuple[np.ndarray, np.ndarray, float], None, None]:
        """Generator of positions, velocities and time.

        Returns
        -------
        * positions : `numpy.ndarray`
            Position of each axis at tai_time.
        * velocities : `numpy.ndarray`
            Velocity of each axis at tai_time.
        * tai_time : `float`
            Current time + advance_time (TAI unix seconds).
        """
        tai = utils.current_tai() + self.advance_time
        start_tai = tai
        end_tai_arr = start_tai + self.period_arr
        max_tai = np.max(end_tai_arr)

        while tai < max_tai:
            frac_period_arr = (tai - start_tai) * self.inverse_period_arr
            angle_rad_arr = 2 * np.pi * frac_period_arr
            position_arr = (
                self.amplitude_arr * np.cos(angle_rad_arr) + self.center_position_arr
            )
            velocity_arr = -self.max_velocity_arr * np.sin(angle_rad_arr)
            for i in range(self.naxes):
                if self.amplitude_arr[i] == 0 or frac_period_arr[i] > 1:
                    position_arr[i] = self.end_position_arr[i]
                    velocity_arr[i] = 0
            yield (position_arr, velocity_arr, tai)
            tai = utils.current_tai() + self.advance_time

        zero_velocity_arr = np.zeros(self.naxes)
        for i in range(self.nextra):
            yield (self.end_position_arr, zero_velocity_arr, tai)
            tai = utils.current_tai() + self.advance_time


class RampGenerator:
    def __init__(
        self,
        start_positions: typing.Sequence[float],
        end_positions: typing.Sequence[float],
        speeds: typing.Sequence[float],
        advance_time: float = 0.1,
        nextra: int = 5,
    ) -> None:
        """Functor to generate position, velocity and TAI time
        for a constant-speed move in one or more axes.

        There is no jerk limiting at the beginning or end.

        This is designed to be used in a CSC commander
        to provide a linear tracking target path.

        Parameters
        ----------
        start_positions : `List` [`float`]
            Start position for each axis (deg).
        end_positions : `List` [`float`]
            End position for each axis (deg).
        speeds : `List` [`float`]
            Speed for each axis (deg/sec). Sign is ignored.
            Speed is ignored for any axis for which start_pos == end_pos.
        advance_time : `float`
            Advance time for returned TAI timestamps (seconds).
            Each returned TAI = current TAI + advance_time.
        nextra : `int`
            Number of extra return values after the moves end.
            These extra values all have position=end_positions
            and velocity=[0, 0, ..., 0]

        Raises
        ------
        ValueError
            If len(start_positions) != len(end_positions),
            len(start_positions) == len(speeds),
            start pos = end pos for all axes,
            or speed = 0 and start pos != end pos for any axis.
        """
        self.naxes = len(start_positions)
        if self.naxes != len(end_positions):
            raise ValueError(
                f"start_positions={start_positions} and end_positions={end_positions} "
                "must have the same length"
            )
        if self.naxes != len(speeds):
            raise ValueError(
                f"start_positions={start_positions} and self.speed_arr={speeds} "
                "must have the same length"
            )
        self.start_position_arr = np.array(start_positions, dtype=float)
        self.end_position_arr = np.array(end_positions, dtype=float)
        self.speed_arr = np.array(speeds, dtype=float)
        self.advance_time = float(advance_time)
        self.nextra = int(nextra)

        delta_position_arr = self.end_position_arr - self.start_position_arr
        # Avoid needless divide-by-zero warnings.
        # Set safe_speed_divisor_arr to an arbitrary nonzero value
        # for those axes with no motion (delta_position_arr == 0).
        safe_speed_divisor_arr = np.where(delta_position_arr == 0, 1, self.speed_arr)
        self.duration_arr = np.abs(delta_position_arr / safe_speed_divisor_arr)
        for i in range(self.naxes):
            if delta_position_arr[i] == 0:
                self.speed_arr[i] = 0
                self.duration_arr[i] = 0
            elif self.speed_arr[i] == 0:
                raise ValueError(f"Speed[{i}] is zero but end pos != start pos")
        if np.max(self.duration_arr) <= 0:
            raise ValueError("Nothing to move")

        # Speed with the correct sign
        self.velocity_arr = np.copysign(self.speed_arr, delta_position_arr)

    @property
    def duration(self) -> float:
        """Duration of the move (seconds), ignoring the extra elements."""
        return np.amax(self.duration_arr)

    def __call__(
        self,
    ) -> typing.Generator[typing.Tuple[np.ndarray, np.ndarray, float], None, None]:
        """Generator of positions, velocities and time.

        Returns
        -------
        * positions : `numpy.ndarray`
            Position of each axis at tai_time.
        * velocities : `numpy.ndarray`
            Velocity of each axis at tai_time.
        * tai_time : `float`
            Current time + advance_time (TAI unix seconds).
        """
        tai = utils.current_tai() + self.advance_time
        start_tai = tai
        end_tai_arr = start_tai + self.duration_arr
        max_tai = np.max(end_tai_arr)

        while tai < max_tai:
            position_arr = self.start_position_arr + self.velocity_arr * (
                tai - start_tai
            )
            velocity_arr = self.velocity_arr.copy()
            for i in range(self.naxes):
                if tai > end_tai_arr[i]:
                    position_arr[i] = self.end_position_arr[i]
                    velocity_arr[i] = 0
            yield (position_arr, velocity_arr, tai)
            tai = utils.current_tai() + self.advance_time

        zero_velocity_arr = np.zeros(self.naxes)
        for i in range(self.nextra):
            yield (self.end_position_arr, zero_velocity_arr, tai)
            tai = utils.current_tai() + self.advance_time
