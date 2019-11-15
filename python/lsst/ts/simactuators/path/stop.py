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

__all__ = ["stop"]

import math

from . import path_segment
from . import path


def stop(tai, pos, vel, max_acceleration):
    """Compute a path to stop as quickly as possible,
    starting from a path of constant velocity.

    Parameters
    ----------
    tai : `float`
        TAI time (unix seconds, e.g. from lsst.ts.salobj.curr_tai()).
    pos : `float`
        Position at time ``tai`` (deg)
    vel : `float`
        Velocity at time ``tai`` (deg/sec)
    max_acceleration : `float`
        Maximum allowed acceleration (deg/sec^2)

    Returns
    -------
    path : `Path`
        A path with 1-2 segments, ending with a segment of zero acceleration
        and zero velocity.

    Raises
    ------
    ValueError
        If |max_acceleration| <= 0
    """
    if max_acceleration <= 0.0:
        raise ValueError(f"max_acceleration={max_acceleration} < 0")

    if vel == 0:
        return path.Path(path_segment.PathSegment(tai=tai, pos=pos),
                         kind=path.Kind.Stopped)

    segments = []
    # tai is large enough that small time changes
    # have poor accuracy, so to improve accuracy
    # compute the path segments using tai = 0
    # then offset all the times before returning the path
    dt = abs(vel)/max_acceleration
    accel = -math.copysign(max_acceleration, vel)
    p1 = pos + dt*(vel + dt*0.5*accel)
    segments.append(path_segment.PathSegment(tai=0, pos=pos,
                                             vel=vel, accel=accel))
    segments.append(path_segment.PathSegment(tai=dt, pos=p1, vel=0))

    for tpvaj in segments:
        tpvaj.tai += tai

    return path.Path(*segments, kind=path.Kind.Stopping)
