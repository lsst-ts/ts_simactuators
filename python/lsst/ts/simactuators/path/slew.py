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


SLEW_FUDGE = 1.05
"""Fudge factor to avoid borderline cases; should be a bit larger than 1"""


def slew(t0, start_pos, start_vel, end_pos, end_vel, max_vel, max_accel):
    """Compute a trapezoidal slew from A to B.

    A and B are paths of constant velocity.

    Parameters
    ----------
    t0 : `float`
        Start time of slew.
    start_pos : `float`
        Position of A at time t0 (deg)
    start_vel : `float`
        Velocity of A at time t0 (deg/sec)
    end_pos : `float`
        Position of B at time t0 (deg)
    end_vel : `float`
        Velocity of B at time t0 (deg/sec)
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
    because that is the shape of the v vs. t curve.

    Here are the initial velocity and constant acceleration for each segment,
    and the duration of that segment, in the notation of this subroutine
    (each of the first 3 segments is present only if its duration is nonzero,
    and the last is absent if the final segment already has zero acceleration):

        segment   v   a   duration
           1      start_vel  a1  dt1
           2      vP  0   dt2
           3      vP  a2  dt3
           4      end_vel  0   unlimited
    """
    if max_vel <= 0.0:
        raise ValueError(f"max_vel={max_vel} < 0")
    if max_accel <= 0.0:
        raise ValueError(f"max_accel={max_accel} < 0")

    # check velocities; errors are: |end_vel| SLEW_FUDGE > max_vel,
    # |start_vel| SLEW_FUDGE > max_vel (can lead to dt1 < 0 for type 2 slews)
    if abs(start_vel) > max_vel*SLEW_FUDGE:
        raise ValueError(f"Telescope is moving too fast (|{start_vel:0.4f}| > {SLEW_FUDGE} * {max_vel}).")
    if abs(end_vel)*SLEW_FUDGE > max_vel:
        raise ValueError(f"Target is moving too fast (|{end_vel:0.4f}| * {SLEW_FUDGE} > {max_vel}; "
                         "telescope cannot acquire it.")

    # compute end_velA, half_end_velAsq, sign_rBAi and sign_end_velA
    # and handle null slews (rBAi and end_velA both zero)
    rBAi = end_pos - start_pos
    end_velA = end_vel - start_vel
    half_end_velAsq = 0.5*end_velA*end_velA
    if rBAi != 0.0 and end_velA != 0.0:
        sign_rBAi = math.copysign(1.0, rBAi)
        sign_end_velA = math.copysign(1.0, end_velA)
    elif rBAi != 0.0:
        sign_rBAi = math.copysign(1.0, rBAi)
        sign_end_velA = sign_rBAi
    elif end_velA != 0.0:
        sign_end_velA = math.copysign(1.0, end_velA)
        sign_rBAi = sign_end_velA
    else:
        return path.Path(path.TPVAJ(t0=t0, pos0=start_pos, vel0=start_vel), kind=path.Kind.Slewing)

    # compute a1 and a3
    # if sign(rBAi) = sign(end_velA), slew is type 1
    # a solution is sure because dt3 has no upper limit over range of soln
    if sign_rBAi == sign_end_velA:
        a1 = sign_rBAi * max_accel

    # else sign(rBAi) = -sign(end_velA) so we use type 2, 3 or 4 slew...

    # a type 2 slew has a maximum dt2 dependent on initial conditions;
    # the biggest dt2 occurs at largest |a|, |a| = max_accel,
    # and smallest |vPB|, |vPB| = |end_velA|
    # so test at that point to see if solutions exist with dt2 > 0
    elif abs(start_vel) * SLEW_FUDGE < max_vel and 0 <= (max_accel*abs(rBAi) - half_end_velAsq):
        a1 = sign_rBAi * max_accel

    # a type 3 slew only exists if max_accel is small enough
    elif max_accel*abs(rBAi)*SLEW_FUDGE <= half_end_velAsq:
        a1 = -sign_rBAi * max_accel

    # a type 4 slew requires reducing accel. to obtain a solution
    else:
        # since the equation for a1 is sure to give |a1| < max_accel
        # (otherwise slew would have been type 3)
        # the equation is guranteed to not overflow
        a1 = -half_end_velAsq / (SLEW_FUDGE*rBAi)
    a3 = -a1

    # make sure velocity / acceleration divisions will not overflow;
    # this is especially important for slew type 4 because acceleration
    # gets reduced, but could also catch stupid max_accel or max_vel inputs
    max_vdiff = max_vel + abs(start_vel) + abs(end_vel)
    if max_vdiff >= min(abs(a1), 1.0)*sys.float_info.max:
        raise RuntimeError('Computed slew time is ridiculous')

    # compute dt2 and vP
    # first assume that dt2 = 0 and compute vP;
    # if resulting vP is too big, reduce it to maximum allowed
    # and compute corresponding increased dt2
    dt2 = 0
    vPB_temp = (0.5*a1*dt2)**2 + half_end_velAsq + a1*rBAi
    if vPB_temp < 0.0:
        raise RuntimeError('Bug! Tried to compute square root of negative value')
    vPB = math.copysign(math.sqrt(vPB_temp), a1) - 0.5*a1*dt2
    vP = vPB + end_vel
    if abs(vP) > max_vel:
        # |vP| is too big, and so must be reduced.
        # Note that |end_vel| < max_vel / SLEW_FUDGE (as tested far above),
        # so |vP| is guaranteed to be reducible to max_vel without causing
        # vPB to approach zero (much less change sign).
        # The division velocity / acceleration was proved safe above.
        # Thus dt2 may be computed without overflow.
        vP = math.copysign(max_vel, vP)
        vPB = vP - end_vel
        dt2 = (rBAi + ((half_end_velAsq - vPB*vPB)/a1)) / vPB

    # compute dt1 and dt3
    # the following divisions were proved safe from overflow above,
    # just after computing a1 and a3
    vPA = vPB + end_velA
    dt1 = vPA/a1
    dt3 = -vPB/a3

    # sanity checks
    if dt1 < 0.0 or dt2 < 0.0 or dt3 < 0.0:
        raise RuntimeError('Bug! Computed negative duration for one or more segments')
    if abs(vP) > max_vel:
        raise RuntimeError('Bug! Computed velocity greater than max velocity')
    if abs(a1) > max_accel:
        raise RuntimeError('Bug! Computed acceleration greater than max acceleration')

    # t0 is typically large enough that small time changes
    # have poor accuracy, so to improve accuracy
    # compute the path segments using t0 = 0
    # then offset all the times before returning the path
    tpvajs = []
    if dt1 > 0:
        tpvajs.append(path.TPVAJ(t0=0, pos0=start_pos, vel0=start_vel, accel0=a1))
        t20 = dt1
        p20, v20 = tpvajs[-1].pva(t20)[0:2]
    else:
        t20 = 0
        p20 = start_pos
        v20 = start_vel
    if dt2 > 0:
        tpvajs.append(path.TPVAJ(t0=t20, pos0=p20, vel0=v20))
        t30 = t20 + dt2
        p30, v30 = tpvajs[-1].pva(t30)[0:2]
    else:
        t30 = t20
        p30 = p20
        v30 = v20
    if dt3 > 0:
        tpvajs.append(path.TPVAJ(t0=t30, pos0=p30, vel0=v30, accel0=a3))
        t40 = t30 + dt3
        p40, v40 = tpvajs[-1].pva(t40)[0:2]
    else:
        t40 = t30
        p40 = p30
        v40 = v30
    if tpvajs[-1].accel0 != 0:
        # if the last TPVAJ has non-zero acceleraton then append
        # a segment with zero acceleration
        tpvajs.append(path.TPVAJ(t0=t40, pos0=p40, vel0=v40))
    for tpvaj in tpvajs:
        tpvaj.t0 += t0

    return path.Path(*tpvajs, kind=path.Kind.Slewing)
