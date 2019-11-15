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
        Segments in the path. Times must be in increasing order.
    kind : `Kind`
        Kind of path

    Raises
    ------
    ValueError
        If no segments are supplied
        or the segments do not have increasing ``tai``.
    """
    def __init__(self, *segments, kind):
        if len(segments) < 1:
            raise ValueError(f"segments={segments} needs at least one element")
        prev_tai = None
        for segment in segments:
            if prev_tai is not None:
                if segment.tai <= prev_tai:
                    raise ValueError(f"segment start times not in increasing order")
                prev_tai = segment.tai

        self.segments = segments
        self.kind = Kind(kind)
        self.tais = [segment.tai for segment in segments]

    def at(self, tai):
        """Compute a path segment at the specified time.

        Parameters
        ----------
        tai : `float`
            TAI time (unix seconds, e.g. from lsst.ts.salobj.curr_tai()).

        Returns
        -------
        segment : `PathSegment`
            Path segement at the specified time, extrapolated if necessary.
        """
        ind = bisect.bisect(self.tais, tai)
        if ind > 0:
            ind -= 1
        return self.segments[ind].at(tai)

    def __len__(self):
        return len(self.segments)

    def __getitem__(self, ind):
        """Indexed access to the PVATs that make up the path."""
        return self.segments[ind]

    def __repr__(self):
        segments_str = ", ".join(repr(segment) for segment in self.segments[:-1])
        return f"Path({segments_str}, kind={self.kind})"
