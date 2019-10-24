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

import math
import sys

__all__ = ["Segment"]


class Segment(object):
    """A path segment of constant jerk; typically used for tracking.

    Also computes maximum ``|velocity|`` and ``|acceleration|`` and,
    optionally, minimum and maximum position during the segment.

    Parameters
    ----------
    dt : `float`
       Duration of motion (sec)
    start_pos : `float`
       Position of starting point at time t = 0 (deg)
    start_vel : `float`
       Velocity of starting point at time t = 0 (deg/sec)
    end_pos : `float`
       Position of ending point at time t = dt (deg)
    end_vel : `float`
       Velocity of ending point at time t = dt (deg/sec)
    do_pos_lim : `bool`
       Compute minimum and maximum position?

    Raises
    ------
    ValueError
        If ``dt < sys.float_info.min``
    """
    def __init__(self, dt, start_pos, start_vel, end_pos, end_vel, do_pos_lim):
        self.dt = dt
        self.start_pos = start_pos
        self.start_vel = start_vel
        self.end_pos = end_pos
        self.end_vel = end_vel

        dt_sq = dt * dt

        # crude test for overflow of vx, start_accel, and jerk;
        # assumes |numerator| < sqrt(bignum),
        # tests |denominator| > sqrt(smallnum)
        # use <= to simplify unit tests
        if dt <= math.sqrt(sys.float_info.min):
            raise ValueError(f"dt={dt} <= math.sqrt(sys.float_info.min))={math.sqrt(sys.float_info.min)}")

        # compute start_accel, jerk and aB
        vx = (end_pos - start_pos) / dt
        start_accel = (3.0 * vx - (2.0 * start_vel + end_vel)) * 2.0 / dt
        jerk = (((start_vel + end_vel) / 2.0) - vx) * (12.0 / dt_sq)
        aB = start_accel + (jerk * dt)

        # Compute maximum |velocity| (peak_vel); this may occur
        # at the endpoints or at time t_vex = -start_accel/jerk.
        # Compute t_vex and vex = v(t_vex); if t_vex is not in range [0, dt),
        # set t_vex = 0, so that vex = start_vel
        if abs(start_accel) < abs(jerk * dt):
            t_vex = max(-start_accel/jerk, 0.0)
        else:
            t_vex = 0.0
        vex = start_vel + t_vex * (start_accel + (t_vex / 2.0) * jerk)

        self.start_accel = start_accel
        self.jerk = jerk
        self.min_pos = None
        self.max_pos = None
        self.peak_vel = max(abs(start_vel), abs(end_vel), abs(vex))
        self.peak_accel = max(abs(start_accel), abs(aB))

        # If desired, compute minimum and maximum position (min_pos & max_pos)
        # min_pos and max_pos may occur at the endpoints or at times t_pex1
        # or t_pex2 (the two solutions to the quadratic equation v(t) = 0).
        # Note that lim (jerk->0) t_pexArr = -start_vel/start_accel,
        # yet the equation used below is ill-behaved at small jerk.
        # Also, it is hard to distinguish t_pexArr not in range [0,dt]
        # and t_pexArr not computable. Rather than try, I simply use
        #  -start_vel/start_accel whenever the standard equation
        # would not give me a reasonable answer.
        if do_pos_lim:
            t_pexArr = [0]*2
            numArr = [0]*2
            pexArr = [0]*2

            # compute the two times t_pexArr,
            # and positions pexArr = p(t_pexArr);
            # if a t_pexArr is out of range [0, dt), set it to 0
            # (so its pexArr = start_pos)
            if abs(start_vel) < abs(start_accel * dt):
                t_pex_zeroj = max(-start_vel / start_accel, 0.0)
            else:
                t_pex_zeroj = 0.0
            sqrt_arg = (start_accel * start_accel) - (2.0 * start_vel * jerk)
            if sqrt_arg < 0.0:
                t_pexArr[0] = 0.0
                t_pexArr[1] = 0.0
            else:
                sqrt_val = math.sqrt(sqrt_arg)
                numArr[0] = -start_accel - sqrt_val
                numArr[1] = -start_accel + sqrt_val
                for branch in range(2):
                    if abs(numArr[branch]) < abs(jerk * dt):
                        t_pexArr[branch] = max(0.0, numArr[branch] / jerk)
                    else:
                        t_pexArr[branch] = t_pex_zeroj
            for branch in range(2):
                pexArr[branch] = start_pos + t_pexArr[branch] \
                    * (start_vel + (t_pexArr[branch] / 2.0) * (start_accel + (t_pexArr[branch] / 3.0) * jerk))

            self.min_pos = min(start_pos, end_pos, pexArr[0], pexArr[1])
            self.max_pos = max(start_pos, end_pos, pexArr[0], pexArr[1])

    def __repr__(self):
        return f"Segment(dt={self.dt}, start_pos={self.start_pos}, start_vel={self.start_vel}, " \
               f"end_pos={self.end_pos}, end_vel={self.end_vel}, do_pos_lim={self.do_pos_lim})"
