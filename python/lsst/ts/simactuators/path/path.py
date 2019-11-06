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

__all__ = ["Kind", "Path"]

import bisect
import enum


class Kind(enum.Enum):
    """Kind of path.
    """
    Stopped = enum.auto()
    Tracking = enum.auto()
    Slewing = enum.auto()
    Stopping = enum.auto()


class Path:
    """A path defined by a sequence of one or more `PathSegment`s.

    Parameters
    ----------
    segments : ``iterable`` of `PathSegment`
        Segments in the path. For the `pva` method to work correctly,
        times must be in increasing order, but this is not checked.
    kind : `Kind`
        Kind of path
    """
    def __init__(self, *segments, kind):
        if len(segments) < 1:
            raise RuntimeError(f"segments={segments} needs at least one element")
        self.segments = segments
        self.kind = Kind(kind)
        self.ts = [segment.start_time for segment in segments]

    def pva(self, t):
        """Compute position, velocity and acceleration at a given time.

        Parameters
        ----------
        t : `float`
            Time (TAI unix seconds, e.g. from lsst.ts.salobj.curr_tai()).

        Returns
        -------
        pva : `PosVelAccel`
            Position, velocity and acceleration at time ``t``,
            extrapolated if necessary.
        """
        ind = bisect.bisect(self.ts, t)
        if ind > 0:
            ind -= 1
        return self.segments[ind].pva(t)

    def __len__(self):
        return len(self.segments)

    def __getitem__(self, ind):
        """Indexed access to the PVATs that make up the path."""
        return self.segments[ind]

    def __repr__(self):
        segments_str = ", ".join(repr(segment) for segment in self.segments[:-1])
        return f"Path({segments_str}, kind={self.kind})"
