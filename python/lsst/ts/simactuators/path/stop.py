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


def stop(start_time, start_pos, start_vel, max_accel):
    """Compute a path to stop as quickly as possible.

    Parameters
    ----------
    start_time : `float`
        Start time of path.
    start_pos : `float`
        Starting position (deg)
    start_vel : `float`
        Starting velocity (deg/sec)
    max_accel : `float`
        Maximum allowed acceleration (deg/sec^2)

    Returns
    -------
    path : `Path`
        A path with 1-2 segments, ending with a segment of zero acceleration
        and zero velocity.

    Raises
    ------
    ValueError if |max_accel| <= 0
    """
    if max_accel <= 0.0:
        raise ValueError(f"max_accel={max_accel} < 0")

    if start_vel == 0:
        return path.Path(path_segment.PathSegment(start_time=start_time, start_pos=start_pos),
                         kind=path.Kind.Stopped)

    segments = []
    # start_time is typically large enough that small time changes
    # have poor accuracy, so to improve accuracy
    # compute the path segments using start_time = 0
    # then offset all the times before returning the path
    dt = abs(start_vel)/max_accel
    accel = -math.copysign(max_accel, start_vel)
    p1 = start_pos + dt*(start_vel + dt*0.5*accel)
    segments.append(path_segment.PathSegment(start_time=0, start_pos=start_pos,
                                             start_vel=start_vel, start_accel=accel))
    segments.append(path_segment.PathSegment(start_time=dt, start_pos=p1, start_vel=0))

    for tpvaj in segments:
        tpvaj.start_time += start_time

    return path.Path(*segments, kind=path.Kind.Stopping)
