from __future__ import annotations

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

__all__ = ["TrackingActuator"]

import typing

from lsst.ts import utils

from . import path as path_m


class TrackingActuator:
    """Simulate an actuator that slews to and tracks a path defined by
    regular calls to `set_target`.

    Parameters
    ----------
    min_position : `float`
        Minimum allowed position
    max_position : `float`
        Maximum allowed position
    max_velocity : `float`
        Maximum allowed velocity (position units/second)
    max_acceleration : `float`
        Maximum allowed acceleration (position units/second^2)
    dtmax_track : `float`
        Maximum allowed time interval (tai - time of last segment)
        for `set_target` to compute a tracking path (second).
        If this limit is not met then `set_target` computes a slewing path.
        This should be larger than the maximum expected time between calls to
        `set_target`, but not much more than that.
    nsettle : `int`, optional
        Number of calls to `set_target` after a slew finishes
        (meaning ``self.path.kind`` is tracking)
        before ``self.kind(tai)`` reports tracking instead of slewing.
    tai : `float`, optional
        TAI time for ``self.target`` and ``self.path``
        (unix seconds, e.g. from lsst.ts.utils.current_tai()).
        If None then use current TAI.
        This is primarily for unit tests; None is usually what you want.
    start_position : `float` or `None`
        Initial position. If `None` use 0 if 0 is in range
        ``[min_position, max_position]`` else use ``min_position``.

    Raises
    ------
    ValueError
        If min_position >= max_position,
        max_velocity <= 0,
        max_acceleration <= 0,
        start_position is not None and  start_position < min_position,
        or start_position > max_position.

    Notes
    -----
    Attributes:

    * ``target``: target set by `set_target` (a `PathSegment`).
    * ``path``: actual actuator path (a `Path`).
    """

    Kind = path_m.Kind

    def __init__(
        self,
        min_position: float,
        max_position: float,
        max_velocity: float,
        max_acceleration: float,
        dtmax_track: float,
        nsettle: int = 2,
        tai: typing.Optional[float] = None,
        start_position: typing.Optional[float] = None,
    ) -> None:
        if min_position >= max_position:
            raise ValueError(
                f"min_position={min_position} must be < max_position={max_position}"
            )
        if max_velocity <= 0:
            raise ValueError(f"max_velocity={max_velocity} must be > 0")
        if max_acceleration <= 0:
            raise ValueError(f"max_velocity={max_velocity} must be > 0")
        self.min_position = min_position
        self.max_position = max_position
        self.max_velocity = max_velocity
        self.max_acceleration = max_acceleration
        self.dtmax_track = dtmax_track
        self.nsettle = nsettle
        if start_position is None:
            if min_position <= 0 <= max_position:
                start_position = 0
            else:
                start_position = min_position
        if start_position < min_position:
            raise ValueError(
                f"start_position={start_position} < min_position={min_position}"
            )
        if start_position > max_position:
            raise ValueError(
                f"start_position={start_position} > max_position={max_position}"
            )

        if tai is None:
            tai = utils.current_tai()
        self.target = path_m.PathSegment(tai=tai, position=start_position)
        self.path = path_m.Path(
            path_m.PathSegment(tai=tai, position=start_position), kind=self.Kind.Stopped
        )
        self._ntrack = 0

    def set_target(self, tai: float, position: float, velocity: float) -> None:
        """Set the target position, velocity and time.

        The actuator will track, if possible, else slew to match the specified
        path.

        Parameters
        ----------
        tai : `float`
            TAI time (unix seconds, e.g. from lsst.ts.utils.current_tai()).
        position : `float`
            Position (deg)
        velocity : `float`
            Velocity (deg/sec)

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
        * The tracking segment path obeys the position, velocity
          and acceleration limits.
        """
        new_path = self._compute_path(tai=tai, position=position, velocity=velocity)
        self.target = path_m.PathSegment(tai=tai, position=position, velocity=velocity)
        self.path = new_path

    @property
    def path(self) -> path_m.Path:
        """Get or set the actuator path, a `path.Path`."""
        return self._path

    @path.setter
    def path(self, path: path_m.Path) -> None:
        self._path = path
        if path.kind == self.Kind.Tracking:
            self._ntrack += 1
        else:
            self._ntrack = 0

    def stop(self, tai: typing.Optional[float] = None) -> None:
        """Stop the axis using maximum acceleration.

        Update the commanded position to match the end point of the stop.

        Parameters
        ----------
        tai : `float`, optional
            TAI time for ``self.target`` and ``self.path``
            (unix seconds, e.g. from lsst.ts.utils.current_tai()).
            If None then use current TAI.
            This is primarily for unit tests; None is usually what you want.
        """
        if tai is None:
            tai = utils.current_tai()
        curr_segment = self.path.at(tai)
        self.path = path_m.stop(
            tai=tai,
            position=curr_segment.position,
            velocity=curr_segment.velocity,
            max_acceleration=self.max_acceleration,
        )
        self.target = self.path[-1]

    def abort(
        self,
        tai: typing.Optional[float] = None,
        position: typing.Optional[float] = None,
    ) -> None:
        """Stop motion immediately, with infinite acceleration.

        Do not change the commanded position.

        Parameters
        ----------
        tai : `float`, optional
            TAI time for ``self.target`` and ``self.path``
            (unix seconds, e.g. from lsst.ts.utils.current_tai()).
            If None then use current TAI.
            This is primarily for unit tests; None is usually what you want.
        position : `float`, optional
            Position at which to stop (deg); if `None` then stop at position
            at time ``tai``.
        """
        if tai is None:
            tai = utils.current_tai()
        if position is None:
            position = self.path.at(tai).position
        self.path = path_m.Path(
            path_m.PathSegment(tai=tai, position=position), kind=self.Kind.Stopped
        )

    def kind(self, tai: typing.Optional[float] = None) -> path_m.Kind:
        """Kind of path at the specified time.

        Parameters
        ----------
        tai : `float`, optional
            TAI time at which to evaluate the kind of path
            (TAI unix seconds, e.g. from lsst.ts.utils.current_tai()).
            If None then use current TAI.
            Ignored unless stopping.

        The result will always match ``self.path.kind`` except as follows:

        - After a slew we report ``path.kind.Slewing`` until ``nsettle``
          consecutive calls to `set_target` result in a path that is tracking.
        - If self.path.kind is stopping and tai > start time of the
          last segment, then the kind is reported as stopped.
        """
        if tai is None:
            tai = utils.current_tai()
        if self.path.kind == self.Kind.Tracking:
            if self._ntrack > self.nsettle:
                return self.Kind.Tracking
            else:
                return self.Kind.Slewing
        elif self.path.kind == self.Kind.Stopping and tai > self.path[-1].tai:
            return self.Kind.Stopped
        return self.path.kind

    def _compute_path(
        self, tai: float, position: float, velocity: float
    ) -> path_m.Path:
        """Compute a trajectory path to the specified target position,
        velocity and time.

        The actuator will track, if possible, else slew to match the target
        path.

        Parameters
        ----------
        tai : `float`
            TAI time (unix seconds, e.g. from lsst.ts.utils.current_tai()).
        position : `float`
            Position (deg)
        velocity : `float`
            Velocity (deg/sec)

        Raises
        ------
        ValueError
            If ``tai <= self.target.tai``,
            where ``self.target.tai`` is the time of
            the previous call to `set_target`.
        """
        prev_tai = self.target.tai  # last commanded time
        dt = tai - prev_tai
        newcurr = None
        if dt <= 0:
            raise ValueError(f"New tai = {tai} <= previous target tai = {prev_tai}")
        if dt < self.dtmax_track:
            # Try tracking.
            prev_segment = self.path.at(prev_tai)
            tracking_segment = path_m.PathSegment.from_end_conditions(
                start_tai=prev_tai,
                start_position=prev_segment.position,
                start_velocity=prev_segment.velocity,
                end_tai=tai,
                end_position=position,
                end_velocity=velocity,
            )
            limits = tracking_segment.limits(tai)
            if (
                limits.max_velocity <= self.max_velocity
                and limits.max_acceleration <= self.max_acceleration
                and limits.min_position >= self.min_position
                and limits.max_position <= self.max_position
            ):
                # Tracking works.
                newcurr = path_m.Path(
                    tracking_segment,
                    path_m.PathSegment(tai=tai, position=position, velocity=velocity),
                    kind=self.Kind.Tracking,
                )

        if newcurr is None:
            # Tracking didn't work, so slew.
            curr_segment = self.path.at(tai)
            newcurr = path_m.slew(
                tai=tai,
                start_position=curr_segment.position,
                start_velocity=curr_segment.velocity,
                end_position=position,
                end_velocity=velocity,
                max_velocity=self.max_velocity,
                max_acceleration=self.max_acceleration,
            )

        return newcurr
