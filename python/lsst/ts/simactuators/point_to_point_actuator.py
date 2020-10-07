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


from lsst.ts import salobj
from . import base_point_to_point_actuator


class PointToPointActuator(base_point_to_point_actuator.BasePointToPointActuator):
    """Simulated actuator that moves to a specified position at constant
    velocity and halts.

    Parameters
    ----------
    min_position : `float`
        Minimum allowed position.
    max_position : `float`
        Maximum allowed position.
    speed : `float`
        Speed of motion.
    start_position : `float` or `None`, optional
        Initial position. If `None` use 0 if 0 is in range
        ``[min_position, max_position]`` else use ``min_position``.

    Raises
    ------
    ValueError
        If speed <= 0,
        min_position >= max_position,
        start_position`` is not None and start_position < min_position
        or start_position > max_position.
    """

    def __init__(self, min_position, max_position, speed, start_position=None):
        if min_position >= max_position:
            raise ValueError(
                f"min_position={min_position} must be < max_position={max_position}"
            )
        if start_position is None:
            if min_position <= 0 <= max_position:
                start_position = 0
            else:
                start_position = min_position
        if not min_position <= start_position <= max_position:
            raise ValueError(
                f"start_position={start_position} must be in range "
                f"[{min_position}, {max_position}]"
            )

        self.min_position = min_position
        self.max_position = max_position
        super().__init__(start_position=start_position, speed=speed)

    def set_position(self, position, start_tai=None):
        """Set a new target position and return the move duration.

        Parameters
        ----------
        position : `float`
            Target position.
        start_tai : `float` or `None`, optional
            TAI date (unix seconds) of the start of the move.
            If `None` use the current time.

        Raises
        ------
        ValueError
            If position < self.min_position or > self.max_position.
        """
        if position < self.min_position or position > self.max_position:
            raise ValueError(
                f"position={position} not in range [{self.min_position}, {self.max_position}]"
            )
        if start_tai is None:
            start_tai = salobj.current_tai()
        start_position = self.position(start_tai)
        return self._set_position(
            start_position=start_position, start_tai=start_tai, end_position=position
        )
