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

from . import path

from lsst.ts import salobj


class TrackingActuator:
    """Simulate an actuator that slews to and tracks a path defined
    by regular position, velocity, time updates.

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
        Maximum allowed tB-tA for `set_cmd` to compute a tracking path (sec);
        if this limit is not met then `set_cmd` computes a slewing path.
        This should be larger than the maximum expected time between calls to
        `set_cmd`, but not much more than that.
    nsettle : `int` (optional)
        Number of calls to `set_cmd` after a slew finishes
        (meaning ``self.curr.kind`` is tracking)
        before ``self.kind(t)`` reports tracking instead of slewing.
    t : `float` (optional)
        TAITime for initial `cmd` `TPVAJ` and `curr` `Path` (TAI seconds);
        if None then use current TAI.
        This is primarily for unit tests; None is usually what you want.

    Raises
    ------
    ValueError
        If ``min_pos >= max_pos``, ``max_vel <= 0``, or ``max_accel <= 0``.
    """
    Kind = path.Kind

    def __init__(self, min_pos, max_pos, max_vel, max_accel, dtmax_track, nsettle=2, t0=None):
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

        if t0 is None:
            t0 = salobj.current_tai()
        if min_pos <= 0 and 0 < max_pos:
            pos0 = 0
        else:
            pos0 = min_pos
        self.cmd = path.TPVAJ(t0=t0, pos0=pos0)
        self.curr = path.Path(path.TPVAJ(t0=t0, pos0=pos0), kind=self.Kind.Stopped)
        self._ntrack = 0

    def set_cmd(self, pos, vel, t):
        """Set commanded position and velocity.

        Parameters
        ----------
        pos : `float`
            Position (deg)
        vel : `float`
            Velocity (deg/sec)
        t : `float`
            Time as a unix time, e.g. from `lsst.ts.salobj.current_tai` (sec)
        """
        tA = self.cmd.t0  # last commanded time
        dt = t - tA
        newcurr = None
        if dt <= 0:
            raise RuntimeError(f"New t = {t} <= previous cmd t = {tA}")
        if dt < self.dtmax_track:
            # try tracking
            pcurr_tA, vcurr_tA = self.curr.pva(tA)[0:2]
            segment = path.Segment(dt=dt, start_pos=pcurr_tA, end_pos=pos, start_vel=vcurr_tA,
                                   end_vel=vel, do_pos_lim=True)
            if segment.peak_vel <= self.max_vel and segment.peak_accel <= self.max_accel \
                    and segment.min_pos >= self.min_pos and segment.max_pos <= self.max_pos:
                # tracking works
                newcurr = path.Path(
                    path.TPVAJ(t0=tA, pos0=pcurr_tA, vel0=vcurr_tA,
                               accel0=segment.start_accel, jerk=segment.jerk),
                    path.TPVAJ(t0=t, pos0=pos, vel0=vel),
                    kind=self.Kind.Tracking)

        if newcurr is None:
            # tracking didn't work, so slew
            pcurr_t, vcurr_t = self.curr.pva(t)[0:2]
            newcurr = path.slew(t0=t, start_pos=pcurr_t, start_vel=vcurr_t, end_pos=pos, end_vel=vel,
                                max_vel=self.max_vel, max_accel=self.max_accel)
        self.cmd = path.TPVAJ(t0=t, pos0=pos, vel0=vel)
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

    def stop(self, t=None):
        """Stop the axis using maximum acceleration.

        Update the commanded position to match the end point of the stop.

        Parameters
        ----------
        t : `float` (optional)
            Time for initial `cmd` `TPVAJ` and `curr` `Path`;
            if None then use current TAI.
            This is primarily for unit tests; None is usually what you want.
        """
        if t is None:
            t = salobj.current_tai()
        pos0, vel0 = self.curr.pva(t)[0:2]
        self.curr = path.stop(pos0=pos0, vel0=vel0, t0=t, max_accel=self.max_accel)
        self.cmd = self.curr[-1]

    def abort(self, t=None, pos=None):
        """Stop motion immediately.

        Do not change the commanded position.

        Parameters
        ----------
        t : `float` (optional)
            Time for initial `cmd` `TPVAJ` and `curr` `Path`;
            if None then use current TAI.
            This is primarily for unit tests; None is usually what you want.
        pos : `float` (optional)
            Position at which to stop (deg); if `None` then stop at position
            at time ``t``.
        """
        if t is None:
            t = salobj.current_tai()
        if pos is None:
            pos = self.curr.pva(t)[0]
        self.curr = path.Path(path.TPVAJ(t0=t, pos0=pos), kind=self.Kind.Stopped)

    def kind(self, t=None):
        """Kind of path we are currently following.

        The answer will match ``self.curr.kind`` except as follows:

        - After a slew we report ``path.kind.Slewing`` until ``nsettle``
          consecutive calls to `set_cmd` result in a path that is tracking.
        - if self.curr.kind is stopping and t > start time of last segment,
          the kind is reported as stopped.
        """
        if t is None:
            t = salobj.current_tai()
        if self.curr.kind == self.Kind.Tracking:
            if self._ntrack > self.nsettle:
                return self.Kind.Tracking
            else:
                return self.Kind.Slewing
        elif self.curr.kind == self.Kind.Stopping and t > self.curr[-1].t0:
            return self.Kind.Stopped
        return self.curr.kind
