# This file is part of ts_simactuators.
#
# Developed for the LSST Telescope and Site Systems.
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

__all__ = ["PointToPointActuator"]

import time


class PointToPointActuator:
    """Simulated actuator that moves to a specified position at constant
    velocity and halts.

    Parameters
    ----------
    min_position : `float`
        Minimum allowed position.
    max_position : `float`
        Maximum allowed position.
    pos : `float`
        Initial position.
    speed : `float`
        Speed of motion.

    Raises
    ------
    ValueError
        If ``speed <= 0``,
        ``min_position >= max_position``,
        or ``pos`` not in range ``[min_position, max_position]``.
    """
    def __init__(self, min_position, max_position, pos, speed):
        if speed <= 0:
            raise ValueError(f"speed={speed} must be positive")
        if min_position >= max_position:
            raise ValueError(f"min_position={min_position} must be < max_position={max_position}")
        if not min_position <= pos <= max_position:
            raise ValueError(f"pos={pos} must be in range "
                             f"[{min_position}, {max_position}]")

        self.min_position = min_position
        self.max_position = max_position
        self.speed = speed
        self._start_position = pos
        self._end_position = pos
        # End time of move, or 0 if not moving.
        self._end_time = 0

    @property
    def start_position(self):
        """Starting position of move.
        """
        return self._start_position

    @property
    def end_position(self):
        """Ending position of move.
        """
        return self._end_position

    def set_position(self, pos):
        """Set a new desired position.

        Raises
        ------
        ValueError
            If pos < self.min_position or > self.max_position.
        """
        if pos < self.min_position or pos > self.max_position:
            raise ValueError(f"pos={pos} not in range [{self.min_position}, {self.max_position}]")
        self._start_position = self.current_position
        self._end_position = pos
        dtime = abs(self.end_position - self.start_position) / self.speed
        self._end_time = time.monotonic() + dtime

    @property
    def current_position(self):
        """Current position.
        """
        rem_time = self.remaining_time
        if rem_time == 0:
            return self.end_position
        else:
            return self.end_position - self.direction*self.speed*rem_time

    @property
    def direction(self):
        """1 if moving or moved to greater position, -1 otherwise.
        """
        return 1 if self.end_position >= self.start_position else -1

    @property
    def moving(self):
        """Is the axis moving?
        """
        return self.remaining_time > 0

    def stop(self):
        """Stop motion instantly.

        Set end_position to the current position.
        """
        self._end_position = self.current_position
        self._end_time = 0

    @property
    def remaining_time(self):
        """Remaining time for this move (sec)."""
        if self._end_time == 0:
            return 0

        rem_time = self._end_time - time.monotonic()
        if rem_time <= 0:
            self._end_time = 0
            return 0

        return rem_time
