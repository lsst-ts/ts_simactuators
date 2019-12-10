# This file is part of ts_simactuators.
#
# Developed for the LSST Data Management System.
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

from lsst.ts import salobj

from . import path


class TrackingActuator:
    """Simulate an actuator that slews to and tracks a path defined by
    regular calls to `set_target`.

    Parameters
    ----------
    min_position : `float`
        Minimum allowed position (deg)
    max_position : `float`
        Maximum allowed position (deg)
    max_velocity : `float`
        Maximum allowed velocity (deg/sec)
    max_acceleration : `float`
        Maximum allowed acceleration (deg/sec^2)
    dtmax_track : `float`
        Maximum allowed time interval (tai - time of last segment)
        for `set_target` to compute a tracking path (sec); if this limit
        is not met then `set_target` computes a slewing path.
        This should be larger than the maximum expected time between calls to
        `set_target`, but not much more than that.
    nsettle : `int` (optional)
        Number of calls to `set_target` after a slew finishes
        (meaning ``self.current.kind`` is tracking)
        before ``self.kind(tai)`` reports tracking instead of slewing.
    tai : `float` (optional)
        TAI time for ``self.target`` and ``self.current``
        (unix seconds, e.g. from lsst.ts.salobj.curr_tai()).
        If None then use current TAI.
        This is primarily for unit tests; None is usually what you want.

    Raises
    ------
    ValueError
        If ``min_position >= max_position``,
        ``max_velocity <= 0``, or ``max_acceleration <= 0``.

    Notes
    -----
    Attributes:

    * ``target``: target set by `set_target` (a `path.PathSegment`).
    * ``current``: the current path (a `path.Path`).
    """
    Kind = path.Kind

    def __init__(self, min_position, max_position, max_velocity, max_acceleration,
                 dtmax_track, nsettle=2, tai=None):
        if min_position >= max_position:
            raise ValueError(f"min_position={min_position} must be < max_position={max_position}")
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

        if tai is None:
            tai = salobj.current_tai()
        if min_position <= 0 and 0 < max_position:
            pos = 0
        else:
            pos = min_position
        self.target = path.PathSegment(tai=tai, pos=pos)
        self.current = path.Path(path.PathSegment(tai=tai, pos=pos),
                                 kind=self.Kind.Stopped)
        self._ntrack = 0

    def set_target(self, tai, pos, vel):
        """Set the target position, velocity and time.

        The actuator will track, if possible, else slew to match the specified
        path.

        Parameters
        ----------
        tai : `float`
            TAI time (unix seconds, e.g. from lsst.ts.salobj.curr_tai()).
        pos : `float`
            Position (deg)
        vel : `float`
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
        prev_tai = self.target.tai  # last commanded time
        dt = tai - prev_tai
        newcurr = None
        if dt <= 0:
            raise ValueError(f"New tai = {tai} <= previous target tai = {prev_tai}")
        if dt < self.dtmax_track:
            # Try tracking.
            prev_segment = self.current.at(prev_tai)
            tracking_segment = path.PathSegment.from_end_conditions(
                start_tai=prev_tai,
                start_position=prev_segment.pos,
                start_velocity=prev_segment.vel,
                end_tai=tai,
                end_position=pos,
                end_velocity=vel)
            limits = tracking_segment.limits(tai)
            if limits.max_velocity <= self.max_velocity and limits.max_acceleration <= self.max_acceleration \
                    and limits.min_position >= self.min_position and limits.max_position <= self.max_position:
                # Tracking works.
                newcurr = path.Path(tracking_segment,
                                    path.PathSegment(tai=tai, pos=pos, vel=vel),
                                    kind=self.Kind.Tracking)

        if newcurr is None:
            # Tracking didn't work, so slew.
            curr_segment = self.current.at(tai)
            newcurr = path.slew(tai=tai, start_position=curr_segment.pos, start_velocity=curr_segment.vel,
                                end_position=pos, end_velocity=vel,
                                max_velocity=self.max_velocity, max_acceleration=self.max_acceleration)
        self.target = path.PathSegment(tai=tai, pos=pos, vel=vel)
        self.current = newcurr

    @property
    def current(self):
        """Get or set the current path, a `path.Path`."""
        return self._curr

    @current.setter
    def current(self, current):
        self._curr = current
        if current.kind == self.Kind.Tracking:
            self._ntrack += 1
        else:
            self._ntrack = 0

    def stop(self, tai=None):
        """Stop the axis using maximum acceleration.

        Update the commanded position to match the end point of the stop.

        Parameters
        ----------
        tai : `float` (optional)
            TAI time for ``self.target`` and ``self.current``
            (unix seconds, e.g. from lsst.ts.salobj.curr_tai()).
            If None then use current TAI.
            This is primarily for unit tests; None is usually what you want.
        """
        if tai is None:
            tai = salobj.current_tai()
        curr_segment = self.current.at(tai)
        self.current = path.stop(pos=curr_segment.pos, vel=curr_segment.vel, tai=tai,
                                 max_acceleration=self.max_acceleration)
        self.target = self.current[-1]

    def abort(self, tai=None, pos=None):
        """Stop motion immediately, with infinite acceleration.

        Do not change the commanded position.

        Parameters
        ----------
        tai : `float` (optional)
            TAI time for ``self.target`` and ``self.current``
            (unix seconds, e.g. from lsst.ts.salobj.curr_tai()).
            If None then use current TAI.
            This is primarily for unit tests; None is usually what you want.
        pos : `float` (optional)
            Position at which to stop (deg); if `None` then stop at position
            at time ``tai``.
        """
        if tai is None:
            tai = salobj.current_tai()
        if pos is None:
            pos = self.current.at(tai).pos
        self.current = path.Path(path.PathSegment(tai=tai, pos=pos), kind=self.Kind.Stopped)

    def kind(self, tai=None):
        """Kind of path at the specified time.

        Parameters
        ----------
        tai : `float` (optional)
            TAI time at which to evaluate the kind of path
            (TAI unix seconds, e.g. from lsst.ts.salobj.curr_tai()).
            If None then use current TAI.
            Ignored unless stopping.

        The result will always match ``self.current.kind`` except as follows:

        - After a slew we report ``path.kind.Slewing`` until ``nsettle``
          consecutive calls to `set_target` result in a path that is tracking.
        - If self.current.kind is stopping and tai > start time of the
          last segment, then the kind is reported as stopped.
        """
        if tai is None:
            tai = salobj.current_tai()
        if self.current.kind == self.Kind.Tracking:
            if self._ntrack > self.nsettle:
                return self.Kind.Tracking
            else:
                return self.Kind.Slewing
        elif self.current.kind == self.Kind.Stopping and tai > self.current[-1].tai:
            return self.Kind.Stopped
        return self.current.kind
