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
    regular calls to `set_cmd`, specifying position, velocity and time.

    Parameters
    ----------
    min_pos : `float`
        Minimum allowed position (deg)
    max_pos : `float`
        Maximum allowed position (deg)
    max_vel : `float`
        Maximum allowed velocity (deg/sec)
    max_accel : `float`
        Maximum allowed acceleration (deg/sec^2)
    dtmax_track : `float`
        Maximum allowed time interval (tai - time of last segment)
        for `set_cmd` to compute a tracking path (sec); if this limit
        is not met then `set_cmd` computes a slewing path.
        This should be larger than the maximum expected time between calls to
        `set_cmd`, but not much more than that.
    nsettle : `int` (optional)
        Number of calls to `set_cmd` after a slew finishes
        (meaning ``self.curr.kind`` is tracking)
        before ``self.kind(tai)`` reports tracking instead of slewing.
    tai : `float` (optional)
        Initial time for the ``self.cmd`` path segment and ``self.curr`` path
        (TAI unix seconds, e.g. from lsst.ts.salobj.curr_tai()).
        If None then use current TAI.
        This is primarily for unit tests; None is usually what you want.

    Raises
    ------
    ValueError
        If ``min_pos >= max_pos``, ``max_vel <= 0``, or ``max_accel <= 0``.

    Notes
    -----
    Attributes:

    * ``cmd``: commanded path set by `set_cmd` (a `path.PathSegment`).
    * ``curr``: the current path (a `path.Path`).
    """
    Kind = path.Kind

    def __init__(self, min_pos, max_pos, max_vel, max_accel, dtmax_track, nsettle=2, tai=None):
        if min_pos >= max_pos:
            raise ValueError(f"min_pos={min_pos} must be < max_pos={max_pos}")
        if max_vel <= 0:
            raise ValueError(f"max_vel={max_vel} must be > 0")
        if max_accel <= 0:
            raise ValueError(f"max_vel={max_vel} must be > 0")
        self.min_pos = min_pos
        self.max_pos = max_pos
        self.max_vel = max_vel
        self.max_accel = max_accel
        self.dtmax_track = dtmax_track
        self.nsettle = nsettle

        if tai is None:
            tai = salobj.current_tai()
        if min_pos <= 0 and 0 < max_pos:
            pos = 0
        else:
            pos = min_pos
        self.cmd = path.PathSegment(tai=tai, pos=pos)
        self.curr = path.Path(path.PathSegment(tai=tai, pos=pos),
                              kind=self.Kind.Stopped)
        self._ntrack = 0

    def set_cmd(self, tai, pos, vel):
        """Set a commanded position, velocity and time.

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
            If ``tai <= self.cmd.tai``,
            where ``self.cmd.tai`` is the time of
            the previous call to `set_cmd`.

        Notes
        -----
        The actuator will track if the following is true:

        * ``tai - self.cmd.tai < self.dtmax_track``
          where ``self.cmd.tai`` is the time of
          the previous call to `set_cmd`.
        * The tracking segment path obeys the position, velocity
          and acceleration limits.
        """
        prev_tai = self.cmd.tai  # last commanded time
        dt = tai - prev_tai
        newcurr = None
        if dt <= 0:
            raise ValueError(f"New tai = {tai} <= previous cmd tai = {prev_tai}")
        if dt < self.dtmax_track:
            # Try tracking.
            prev_segment = self.curr.at(prev_tai)
            tracking_segment = path.PathSegment.from_end_conditions(
                start_tai=prev_tai,
                start_pos=prev_segment.pos,
                start_vel=prev_segment.vel,
                end_tai=tai,
                end_pos=pos,
                end_vel=vel)
            limits = tracking_segment.limits(tai)
            if limits.max_vel <= self.max_vel and limits.max_accel <= self.max_accel \
                    and limits.min_pos >= self.min_pos and limits.max_pos <= self.max_pos:
                # Tracking works.
                newcurr = path.Path(tracking_segment,
                                    path.PathSegment(tai=tai, pos=pos, vel=vel),
                                    kind=self.Kind.Tracking)

        if newcurr is None:
            # Tracking didn't work, so slew.
            curr_segment = self.curr.at(tai)
            newcurr = path.slew(tai=tai, start_pos=curr_segment.pos, start_vel=curr_segment.vel,
                                end_pos=pos, end_vel=vel,
                                max_vel=self.max_vel, max_accel=self.max_accel)
        self.cmd = path.PathSegment(tai=tai, pos=pos, vel=vel)
        self.curr = newcurr

    @property
    def curr(self):
        """Get or set the current path, a `path.Path`."""
        return self._curr

    @curr.setter
    def curr(self, curr):
        self._curr = curr
        if curr.kind == self.Kind.Tracking:
            self._ntrack += 1
        else:
            self._ntrack = 0

    def stop(self, tai=None):
        """Stop the axis using maximum acceleration.

        Update the commanded position to match the end point of the stop.

        Parameters
        ----------
        tai : `float` (optional)
            Initial time for ``self.cmd`` path segment and ``self.curr`` path
            (TAI unix seconds, e.g. from lsst.ts.salobj.curr_tai()).
            If None then use current TAI.
            This is primarily for unit tests; None is usually what you want.
        """
        if tai is None:
            tai = salobj.current_tai()
        curr_segment = self.curr.at(tai)
        self.curr = path.stop(pos=curr_segment.pos, vel=curr_segment.vel, tai=tai,
                              max_accel=self.max_accel)
        self.cmd = self.curr[-1]

    def abort(self, tai=None, pos=None):
        """Stop motion immediately, with infinite acceleration.

        Do not change the commanded position.

        Parameters
        ----------
        tai : `float` (optional)
            Initial time for ``self.cmd`` path segment and ``self.curr`` path
            (TAI unix seconds, e.g. from lsst.ts.salobj.curr_tai()).
            If None then use current TAI.
            This is primarily for unit tests; None is usually what you want.
        pos : `float` (optional)
            Position at which to stop (deg); if `None` then stop at position
            at time ``tai``.
        """
        if tai is None:
            tai = salobj.current_tai()
        if pos is None:
            pos = self.curr.at(tai).pos
        self.curr = path.Path(path.PathSegment(tai=tai, pos=pos), kind=self.Kind.Stopped)

    def kind(self, tai=None):
        """Kind of path at the specified time.

        Parameters
        ----------
        tai : `float` (optional)
            Time at which to evaluate the kind of path
            (TAI unix seconds, e.g. from lsst.ts.salobj.curr_tai()).
            If None then use current TAI.
            Ignored unless stopping.

        The result will always match ``self.curr.kind`` except as follows:

        - After a slew we report ``path.kind.Slewing`` until ``nsettle``
          consecutive calls to `set_cmd` result in a path that is tracking.
        - if self.curr.kind is stopping and tai > start time of last segment,
          the kind is reported as stopped.
        """
        if tai is None:
            tai = salobj.current_tai()
        if self.curr.kind == self.Kind.Tracking:
            if self._ntrack > self.nsettle:
                return self.Kind.Tracking
            else:
                return self.Kind.Slewing
        elif self.curr.kind == self.Kind.Stopping and tai > self.curr[-1].tai:
            return self.Kind.Stopped
        return self.curr.kind
