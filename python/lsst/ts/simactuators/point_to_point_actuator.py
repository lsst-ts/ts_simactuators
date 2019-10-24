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

    Poll to get the current state.

    Parameters
    ----------
    min_pos : `float`
        Minimum allowed position.
    max_pos : `float`
        Maximum allowed position.
    pos : `float`
        Initial position.
    speed : `float`
        Speed of motion.

    Raises
    ------
    ValueError
        If ``speed <= 0``,
        ``min_pos >= max_pos``,
        or ``pos`` not in range ``[min_pos, max_pos]``.
    """
    def __init__(self, min_pos, max_pos, pos, speed):
        if speed <= 0:
            raise ValueError(f"speed={speed} must be positive")
        if min_pos >= max_pos:
            raise ValueError(f"min_pos={min_pos} must be < max_pos={max_pos}")
        if not min_pos <= pos <= max_pos:
            raise ValueError(f"pos={pos} must be in range "
                             f"[{min_pos}, {max_pos}]")

        self.min_pos = min_pos
        self.max_pos = max_pos
        self.speed = speed
        self._start_pos = pos
        self._end_pos = pos
        # End time of move, or 0 if not moving.
        self._end_time = 0

    @property
    def start_pos(self):
        """Starting position of move.
        """
        return self._start_pos

    @property
    def end_pos(self):
        """Ending position of move.
        """
        return self._end_pos

    def set_pos(self, pos):
        """Set a new desired position.
        """
        if pos < self.min_pos or pos > self.max_pos:
            raise ValueError(f"pos={pos} not in range [{self.min_pos}, {self.max_pos}]")
        self._start_pos = self.curr_pos
        self._end_pos = pos
        dtime = abs(self.end_pos - self.start_pos) / self.speed
        self._end_time = time.monotonic() + dtime

    @property
    def curr_pos(self):
        """Current position.
        """
        rem_time = self.remaining_time
        if rem_time == 0:
            return self.end_pos
        else:
            return self.end_pos - self.direction*self.speed*rem_time

    @property
    def direction(self):
        """1 if moving or moved to greater position, -1 otherwise.
        """
        return 1 if self.end_pos >= self.start_pos else -1

    @property
    def moving(self):
        """Is the axis moving?
        """
        return self.remaining_time > 0

    def stop(self):
        """Stop motion instantly.

        Set end_pos to the current position.
        """
        self._end_pos = self.curr_pos
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
