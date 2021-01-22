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

__all__ = ["CircularTrackingActuator"]

import math

from lsst.ts import salobj

from . import base
from . import path
from . import tracking_actuator


class CircularTrackingActuator(tracking_actuator.TrackingActuator):
    """A version of a TrackingActuator that moves in a circle with no limits.

    Parameters
    ----------
    max_velocity : `float`
        Maximum allowed velocity (degree/second)
    max_acceleration : `float`
        Maximum allowed acceleration (degree/second^2)
    dtmax_track : `float`
        Maximum allowed time interval (tai - time of last segment)
        for `set_target` to compute a tracking path (sec); if this limit
        is not met then `set_target` computes a slewing path.
        This should be larger than the maximum expected time between calls to
        `set_target`, but not much more than that.
    nsettle : `int`, optional
        Number of calls to `set_target` after a slew finishes
        (meaning ``self.path.kind`` is tracking)
        before ``self.kind(tai)`` reports tracking instead of slewing.
    tai : `float`, optional
        TAI time for ``self.target`` and ``self.path``
        (unix seconds, e.g. from lsst.ts.salobj.current_tai()).
        If None then use current TAI.
        This is primarily for unit tests; None is usually what you want.
    start_position : `float` or `None`, optional
        Initial position. If `None` use 0.

    Raises
    ------
    ValueError
        If ``max_velocity <= 0`` or ``max_acceleration <= 0``.

    Notes
    -----
    Attributes:

    * ``target``: target set by `set_target` (a `path.PathSegment`).
    * ``path``: actual actuator path (a `path.Path`).

    `set_target` wraps ``target.position`` into the range [0, 360).
    That is the only guarantee about the wrap of of angles
    in ``target`` and ``path``. Thus the positions in ``path``
    can easily be out of that range, as can the computed ``path`` positions
    at later times.
    """

    def __init__(
        self,
        max_velocity,
        max_acceleration,
        dtmax_track,
        nsettle=2,
        tai=None,
        start_position=None,
    ):
        if start_position is None:
            wrapped_start_position = 0
        else:
            wrapped_start_position = salobj.angle_wrap_nonnegative(start_position).deg
        super().__init__(
            min_position=-math.inf,
            max_position=math.inf,
            max_velocity=max_velocity,
            max_acceleration=max_acceleration,
            dtmax_track=dtmax_track,
            nsettle=nsettle,
            tai=tai,
            start_position=wrapped_start_position,
        )

    def set_target(self, tai, position, velocity, direction=base.Direction.NEAREST):
        """Set the target position, velocity and time.

        The actuator will track, if possible, else slew to match the specified
        path.

        ``target.position`` is wrapped into the range [0, 360).

        Parameters
        ----------
        tai : `float`
            TAI time (unix seconds, e.g. from lsst.ts.salobj.current_tai()).
        position : `float`
            Position (deg)
        velocity : `float`
            Velocity (deg/sec)
        direction : `Direction`
            Desired direction for acquiring the target (which way to slew).
            `Direction.NEAREST` picks the slew with shortest duration.

        Raises
        ------
        ValueError
            If ``tai <= self.target.tai``,
            where ``self.target.tai`` is the time of
            the previous call to `set_target`.

        Notes
        -----
        The actuator will track if the following is true:

        * ``tai - self.target.tai < self.dtmax_track``
          where ``self.target.tai`` is the time of
          the previous call to `set_target`.
        * The tracking segment path obeys the velocity and acceleration limits.
        """
        wrapped_position = salobj.angle_wrap_nonnegative(position).deg
        new_path = self._compute_directed_path(
            tai=tai, position=wrapped_position, velocity=velocity, direction=direction
        )
        self.target = path.PathSegment(
            tai=tai, position=wrapped_position, velocity=velocity
        )
        self.path = new_path

    def _compute_directed_path(self, tai, position, velocity, direction):
        """Compute a path to a target.

        Parameters
        ----------
        tai : `float`
            TAI time (unix seconds, e.g. from lsst.ts.salobj.current_tai()).
        position : `float`
            Position (deg)
        velocity : `float`
            Velocity (deg/sec)
        direction : `Direction`
            Desired overall direction of motion.
            If `Direction.NEAREST` the direction is based on
            the path position at ``tai`` relative to the new position
            (velocit is ignored, to simplify the code).
        """
        current_position = self.path.at(tai).position
        delta_position = salobj.angle_diff(position, current_position).deg

        # wrapped_position is the target position wrapped to be
        # the appropriate delta from the current position
        # (which is computed using the current path,
        # regardless of how many times the actuator may have wrapped,
        # so as to avoid changing self.path in this call).
        if direction == base.Direction.NEAREST and delta_position != 0:
            wrapped_position1 = delta_position + current_position
            if delta_position < 0:
                wrapped_position2 = wrapped_position1 + 360
            else:
                wrapped_position2 = wrapped_position1 - 360
            paths = [
                self._compute_path(
                    tai=tai, position=wrapped_position, velocity=velocity
                )
                for wrapped_position in (wrapped_position1, wrapped_position2)
            ]
            paths.sort(key=lambda path: path[-1].tai)
            return paths[0]
        else:
            if direction == base.Direction.POSITIVE and delta_position < 0:
                delta_position += 360
            elif direction == base.Direction.NEGATIVE and delta_position > 0:
                delta_position -= 360

            wrapped_position = current_position + delta_position

            path = self._compute_path(
                tai=tai, position=wrapped_position, velocity=velocity
            )

            # Offset the computed path positions so they are in the
            # correct wrap relative to the ``position`` argument.
            offset = position - wrapped_position
            for segment in path.segments:
                segment.position += offset

            return path
