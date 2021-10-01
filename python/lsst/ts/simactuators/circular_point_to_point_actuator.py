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

__all__ = ["CircularPointToPointActuator"]

from lsst.ts import utils
from . import base
from . import base_point_to_point_actuator


class CircularPointToPointActuator(
    base_point_to_point_actuator.BasePointToPointActuator
):
    """Simulated circular actuator that moves to a specified position at
    constant velocity and halts.

    Circular means motion around a circle, with no limits,
    such as the azimuth axis for some telescope enclosures.
    Positions are reported in the range [0, 360) degrees.

    Parameters
    ----------
    speed : `float`
        Speed of motion (degrees/second).
    start_position : `float` or `None`, optional
        Initial position (degrees).

    Raises
    ------
    ValueError
        If ``speed <= 0``
    """

    def __init__(self, speed, start_position=None):
        # The private fields _start_position and _end_position are not
        # necessarily in the range [0, 360). This is important for supporting
        # the ``direction`` parameter. The public accessors all wrap
        # position into the range [0, 360).
        if start_position is None:
            start_position = 0
        super().__init__(
            start_position=utils.angle_wrap_nonnegative(start_position).deg,
            speed=speed,
        )

    @property
    def start_position(self):
        """Starting position of move in the range [0, 360) degrees."""
        return utils.angle_wrap_nonnegative(self._start_position).deg

    @property
    def end_position(self):
        """Ending position of move, in the range [0, 360) degrees."""
        return utils.angle_wrap_nonnegative(self._end_position).deg

    def position(self, tai=None):
        """Current position.

        Parameters
        ----------
        tai : `float` or `None`, optional
            TAI date, unix seconds. Current time if `None`.
        """
        return utils.angle_wrap_nonnegative(super().position(tai)).deg

    def set_position(self, position, direction=base.Direction.NEAREST, start_tai=None):
        """Set a new target position and return the move duration.

        Parameters
        ----------
        position : `float`
            Commanded position, in degrees.
        direction : `Direction`, optional
            Direction of motion.
        start_tai : `float` or `None`, optional
            TAI date (unix seconds) of the start of the move.
            If `None` use the current time.
        """
        direction = base.Direction(direction)
        if start_tai is None:
            start_tai = utils.current_tai()
        start_position = utils.angle_wrap_nonnegative(self.position(start_tai)).deg
        delta = utils.angle_diff(position, start_position).deg
        if direction is base.Direction.POSITIVE:
            if delta < 0:
                delta += 360
        elif direction is base.Direction.NEGATIVE:
            if delta > 0:
                delta -= 360
        elif direction is not base.Direction.NEAREST:
            raise ValueError(f"Unsupported direction {direction!r}")
        return self._set_position(
            start_position=start_position,
            start_tai=start_tai,
            end_position=start_position + delta,
        )
