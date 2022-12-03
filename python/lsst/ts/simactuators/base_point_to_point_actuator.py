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

__all__ = ["BasePointToPointActuator"]

import typing

from lsst.ts import utils

from . import base


class BasePointToPointActuator:
    """Base class for point to point actuators that move to a specified
    position at constant velocity and halt.

    Parameters
    ----------
    start_position : `float`
        Initial position.
    speed : `float`
        Speed of motion.

    Raises
    ------
    ValueError
        If ``speed <= 0``

    Notes
    -----
    You may change the ``speed``, but we recommended only doing this
    when the actuator is stationary, because the actuator is not sophisticated
    enough to mimic actually changing actuator speed during a move.
    Changing speed does not change `start_position`, `end_position`,
    `start_tai`, or `end_tai`.

    The default value for all time arguments is current TAI, in unix seconds.
    However, you may use any zero point for time that you like, as long as
    you always specify the time argument. This can be handy in unit tests,
    as it can eliminate the need to wait for an actuator to move.
    """

    def __init__(self, start_position: float, speed: float) -> None:
        if speed <= 0:
            raise ValueError(f"speed={speed} must be positive")

        # Speed of a move
        self.speed = speed
        # Start position of a move.
        self._start_position = start_position
        # End position of a move; the current position if not moving.
        self._end_position = start_position
        # Start time of move.
        self._start_tai = utils.current_tai()
        # End time of move. Record this, rather than computing it on demand,
        # because it allows a user to change the speed without altering
        # the end time.
        self._end_tai = self._start_tai

    @property
    def direction(self) -> base.Direction:
        """Direction of current or most recent move, as a
        `lsst.ts.simactuators.Direction` enum value.

        `lsst.ts.simactuators.Direction.NEGATIVE` if moving or moved in the
        negative direction, `lsst.ts.simactuators.Direction.POSTIVE`
        otherwise, including for a null movement or if never moved.
        """
        return (
            base.Direction.POSITIVE
            if self._end_position >= self._start_position
            else base.Direction.NEGATIVE
        )

    @property
    def end_position(self) -> float:
        """Ending position of move."""
        return self._end_position

    @property
    def end_tai(self) -> float:
        """TAI date at end of move, unix seconds."""
        return self._end_tai

    @property
    def start_position(self) -> float:
        """Starting position of move."""
        return self._start_position

    @property
    def start_tai(self) -> float:
        """TAI date at start of move move recent move."""
        return self._start_tai

    def remaining_time(self, tai: typing.Optional[float] = None) -> float:
        """Remaining time for the move (seconds); 0 after the move.

        Parameters
        ----------
        tai : `float` or `None`, optional
            TAI date, unix seconds. Current time if `None`.
        """
        if tai is None:
            tai = utils.current_tai()
        return max(0, self.end_tai - tai)

    def moving(self, tai: typing.Optional[float] = None) -> bool:
        """Is the axis moving? False before and after the move.

        Parameters
        ----------
        tai : `float` or `None`, optional
            TAI date, unix seconds. Current time if `None`.
        """
        if tai is None:
            tai = utils.current_tai()
        return self.start_tai < tai < self.end_tai

    def position(self, tai: typing.Optional[float] = None) -> float:
        """Actual position.

        Parameters
        ----------
        tai : `float` or `None`, optional
            TAI date, unix seconds. Current time if `None`.
        """
        if tai is None:
            tai = utils.current_tai()
        if tai > self.end_tai:
            return self.end_position
        elif tai < self.start_tai:
            return self.start_position
        rem_time = self.end_tai - tai
        return self._end_position - self.direction * self.speed * rem_time

    def velocity(self, tai: typing.Optional[float] = None) -> float:
        """Actual velocity.

        Parameters
        ----------
        tai : `float` or `None`, optional
            TAI date, unix seconds. Current time if `None`.
        """
        if tai is None:
            tai = utils.current_tai()
        if not self.moving(tai):
            return 0
        return self.speed * self.direction

    def stop(self, tai: typing.Optional[float] = None) -> None:
        """Stop motion instantly.

        Set end_position to the current position,
        or to the position at the specified time.

        Parameters
        ----------
        tai : `float` or `None`, optional
            TAI date, unix seconds. Current time if `None`.
        """
        if tai is None:
            tai = utils.current_tai()
        self._end_position = self.position(tai)
        self._end_tai = tai

    def _set_position(
        self, start_position: float, start_tai: float, end_position: float
    ) -> float:
        """Set start and end positions and times and return move duration."""
        self._start_position = start_position
        self._start_tai = start_tai
        self._end_position = end_position
        duration = abs(self._end_position - self._start_position) / self.speed
        self._end_tai = self._start_tai + duration
        return duration
