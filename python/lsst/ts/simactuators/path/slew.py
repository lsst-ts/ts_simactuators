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

__all__ = ["slew", "SLEW_FUDGE"]

import math
import sys

from . import path
from . import path_segment


SLEW_FUDGE = 1.05
"""Fudge factor to avoid borderline cases; should be a bit larger than 1"""


def slew(tai, start_pos, start_vel, end_pos, end_vel, max_vel, max_accel):
    """Compute a trapezoidal slew from a start path of constant velocity
    to and end path of constant velocity.

    Parameters
    ----------
    tai : `float`
        Start time of slew.
    start_pos : `float`
        Position of A at time ``tai`` (deg)
    start_vel : `float`
        Velocity of A at time ``tai`` (deg/sec)
    end_pos : `float`
        Position of B at time ``tai`` (deg)
    end_vel : `float`
        Velocity of B at time ``tai`` (deg/sec)
    max_vel : `float`
        Maximum allowed velocity (deg/sec)
    max_accel : `float`
        Maximum allowed acceleration (deg/sec^2)

    Returns
    -------
    path : `Path`
        A path with 1-4 segments, ending with a segment of zero acceleration
        that matches B.

    Raises
    ------
    ValueError if any of the following are true:

    * max_vel <= 0
    * max_accel <= 0
    * |start_vel| > max_vel*SLEW_FUDGE
    * |end_vel| > max_vel/SLEW_FUDGE

    RuntimeError if internal consistency checks fail.

    Notes
    -----
    Dies if max_vel or max_accel are so small as to cause under- or over-flow.

    * Details *

    The magic number `SLEW_FUDGE` is used to avoid borderline cases.
    It should be a number a bit bigger than one.

    How it works:
    The slew begins by tracking object A, and ends by tracking object B.
    Both objects are assumed to be moving at constant velocity.

    A trapezoidal slew has three segments, two of constant acceleration
    separated by a constant velocity segment. It is called "trapezoidal"
    because that is the shape of the velocity vs. time curve.

    Here are the initial velocity and constant acceleration for each segment,
    and the duration of that segment, in the notation of this subroutine
    (each of the first 3 segments is present only if its duration is nonzero,
    and the last is absent if the final segment already has zero acceleration):

        segment   vel        accel          duration
           1      start_vel  start_accel    dt1
           2      vel_mid    0              dt2
           3      vel_mid    end_accel      dt3
           4      end_vel    0              unlimited
    """
    if max_vel <= 0.0:
        raise ValueError(f"max_vel={max_vel} < 0")
    if max_accel <= 0.0:
        raise ValueError(f"max_accel={max_accel} < 0")

    # Check velocities; errors are: |end_vel| SLEW_FUDGE > max_vel,
    # |start_vel| SLEW_FUDGE > max_vel (can lead to dt1 < 0 for type 2 slews).
    if abs(start_vel) > max_vel*SLEW_FUDGE:
        raise ValueError(f"Telescope is moving too fast (|{start_vel:0.4f}| > {SLEW_FUDGE} * {max_vel}).")
    if abs(end_vel)*SLEW_FUDGE > max_vel:
        raise ValueError(f"Target is moving too fast (|{end_vel:0.4f}| * {SLEW_FUDGE} > {max_vel}; "
                         "telescope cannot acquire it.")

    # Compute end_velA, half_end_velAsq, sign_dpBAi and sign_end_velA
    # and handle null slews (dpBAi and end_velA both zero).
    # "A" refers to the initial path (start_pos, etc.) and "B" to the end path
    # and "dp" is a change in position.
    dpBAi = end_pos - start_pos
    end_velA = end_vel - start_vel
    half_end_velAsq = 0.5*end_velA*end_velA
    if dpBAi != 0.0 and end_velA != 0.0:
        sign_dpBAi = math.copysign(1.0, dpBAi)
        sign_end_velA = math.copysign(1.0, end_velA)
    elif dpBAi != 0.0:
        sign_dpBAi = math.copysign(1.0, dpBAi)
        sign_end_velA = sign_dpBAi
    elif end_velA != 0.0:
        sign_end_velA = math.copysign(1.0, end_velA)
        sign_dpBAi = sign_end_velA
    else:
        return path.Path(path_segment.PathSegment(tai=tai, pos=start_pos, vel=start_vel),
                         kind=path.Kind.Slewing)

    # Compute start_accel and accel3
    # if sign(dpBAi) = sign(end_velA), slew is type 1
    # a solution is sure because dt3 has no upper limit over range of soln
    if sign_dpBAi == sign_end_velA:
        start_accel = sign_dpBAi * max_accel

    # Else sign(dpBAi) = -sign(end_velA) so we use type 2, 3 or 4 slew...

    # A type 2 slew has a maximum dt2 dependent on initial conditions;
    # the biggest dt2 occurs at largest |a|, |a| = max_accel,
    # and smallest |vPB|, |vPB| = |end_velA|
    # so test at that point to see if solutions exist with dt2 > 0
    elif abs(start_vel) * SLEW_FUDGE < max_vel and 0 <= (max_accel*abs(dpBAi) - half_end_velAsq):
        start_accel = sign_dpBAi * max_accel

    # A type 3 slew only exists if max_accel is small enough.
    elif max_accel*abs(dpBAi)*SLEW_FUDGE <= half_end_velAsq:
        start_accel = -sign_dpBAi * max_accel

    # A type 4 slew requires reducing accel. to obtain a solution.
    else:
        # The equation for start_accel is sure to give
        # |start_accel| < max_accel(otherwise slew would have been type 3),
        # so this equation is guranteed to not overflow
        start_accel = -half_end_velAsq / (SLEW_FUDGE*dpBAi)
    accel3 = -start_accel

    # Make sure velocity / acceleration divisions will not overflow;
    # this is especially important for slew type 4 because acceleration
    # gets reduced, but could also catch stupid max_accel or max_vel inputs.
    # Note that sys.float_info.max is intended to be a tiny delta-time.
    max_vdiff = max_vel + abs(start_vel) + abs(end_vel)
    if max_vdiff >= min(abs(start_accel), 1.0)*sys.float_info.max:
        raise RuntimeError("Computed slew time is ridiculous.")

    # Compute dt2 and vel_mid
    # first assume that dt2 = 0 and compute vel_mid;
    # if resulting vel_mid is too big, reduce it to maximum allowed
    # and compute corresponding increased dt2.
    dt2 = 0
    vPB_temp = (0.5*start_accel*dt2)**2 + half_end_velAsq + start_accel*dpBAi
    if vPB_temp < 0.0:
        raise RuntimeError("Bug! Tried to compute square root of negative value.")
    vPB = math.copysign(math.sqrt(vPB_temp), start_accel) - 0.5*start_accel*dt2
    vel_mid = vPB + end_vel
    if abs(vel_mid) > max_vel:
        # |vel_mid| is too big, and so must be reduced.
        # Note that |end_vel| < max_vel / SLEW_FUDGE (as tested far above),
        # so |vel_mid| is guaranteed to be reducible to max_vel without causing
        # vPB to approach zero (much less change sign).
        # The division velocity / acceleration was proved safe above.
        # Thus dt2 may be computed without overflow.
        vel_mid = math.copysign(max_vel, vel_mid)
        vPB = vel_mid - end_vel
        dt2 = (dpBAi + ((half_end_velAsq - vPB*vPB)/start_accel)) / vPB

    # Compute dt1 and dt3. Note that the following divisions
    # were proved safe from overflow above.
    vPA = vPB + end_velA
    dt1 = vPA/start_accel
    dt3 = -vPB/accel3

    # Perform sanity checks.
    if dt1 < 0.0 or dt2 < 0.0 or dt3 < 0.0:
        raise RuntimeError("Bug! Computed negative duration for one or more segments.")
    if abs(vel_mid) > max_vel:
        raise RuntimeError("Bug! Computed velocity greater than max velocity.")
    if abs(start_accel) > max_accel:
        raise RuntimeError("Bug! Computed acceleration greater than max acceleration.")

    # Start_time is typically large enough that small time changes
    # have poor accuracy, so to improve accuracy
    # compute the path segments using tai = 0
    # then offset all the times before returning the path.
    segments = []
    if dt1 > 0:
        # There is a segment 1 with acceleration start_accel
        segments.append(path_segment.PathSegment(tai=0, pos=start_pos,
                                                 vel=start_vel, accel=start_accel))
        segment2 = segments[-1].at(dt1)
    else:
        segment2 = path_segment.PathSegment(tai=segment2.tai, pos=start_pos, vel=start_vel)
    if dt2 > 0:
        # There is a segment 2 with no acceleration
        segments.append(path_segment.PathSegment(tai=segment2.tai, pos=segment2.pos, vel=segment2.vel))
        segment3 = segments[-1].at(segment2.tai + dt2)
    else:
        segment3 = segment2
    if dt3 > 0:
        # There is a segment 3 with acceleration accel3
        segments.append(path_segment.PathSegment(tai=segment3.tai, pos=segment3.pos,
                                                 vel=segment3.vel, accel=accel3))
        segment4 = segments[-1].at(segment3.tai + dt3)
    else:
        segment4 = segment3
    if segments[-1].accel != 0:
        # If the last PathSegment has non-zero acceleraton then append
        # a segment with zero acceleration.
        segments.append(path_segment.PathSegment(tai=segment4.tai, pos=segment4.pos, vel=segment4.vel))
    for segment in segments:
        segment.tai += tai

    return path.Path(*segments, kind=path.Kind.Slewing)
