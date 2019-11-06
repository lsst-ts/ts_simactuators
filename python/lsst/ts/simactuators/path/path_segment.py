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

__all__ = ["PosVelAccel", "MotionLimits", "PathSegment"]


class PosVelAccel:
    """Position, velocity and acceleration.

    Parameters
    ----------
    pos : `float` (optional)
        Position (deg)
    vel : `float` (optional)
        Velocity (deg/sec)
    accel : `float` (optional)
        Acceleration (deg/sec^2)
    """
    def __init__(self, pos=0, vel=0, accel=0):
        self.pos = pos
        self.vel = vel
        self.accel = accel


class MotionLimits:
    """Limits of motion.

    Parameters
    ----------
    min_pos : `float`
        Minimum position (deg)
    max_pos : `float`
        Maximum position (deg)
    max_vel : `float`
        Maximum absolute value of velocity (deg/sec)
    max_accel : `float`
        Maximum absolute value of acceleration (deg/sec/sec)
    """
    def __init__(self, min_pos, max_pos, max_vel, max_accel):
        self.min_pos = min_pos
        self.max_pos = max_pos
        self.max_vel = max_vel
        self.max_accel = max_accel


class PathSegment:
    """A segment of a path, motion with constant jerk.

    Parameters
    ----------
    start_time : `float`
        Initial time (TAI unix seconds, e.g. from lsst.ts.salobj.curr_tai()).
    start_pos : `float` (optional)
        Initial position (deg)
    start_vel : `float` (optional)
        Initial velocity (deg/sec)
    start_accel : `float` (optional)
        Initial acceleration (deg/sec^2)
    jerk : `float` (optional)
        Jerk (deg/sec^3)
    """
    def __init__(self, start_time, start_pos=0, start_vel=0, start_accel=0, jerk=0):
        self.start_time = float(start_time)
        self.start_pos = float(start_pos)
        self.start_vel = float(start_vel)
        self.start_accel = float(start_accel)
        self.jerk = float(jerk)

    @classmethod
    def from_end_conditions(cls, start_time, start_pos, start_vel, end_time, end_pos, end_vel):
        """Create a path segment from end conditions.

        Parameters
        ----------
        start_time : `float`
            Initial time (TAI unix seconds, e.g. from
            lsst.ts.salobj.curr_tai()).
        start_pos : `float` (optional)
            Initial position (deg)
        start_vel : `float` (optional)
            Initial velocity (deg/sec)
        end_time : `float`
            Final time (TAI unix seconds, e.g. from lsst.ts.salobj.curr_tai()).
        end_pos : `float` (optional)
            Final position (deg)
        end_vel : `float` (optional)
            Final velocity (deg/sec)

        Raises
        ------
        ValueError
            If end_time - start_time <= 0.
        """
        dt = end_time - start_time
        # Avoid overflow
        if dt <= math.sqrt(sys.float_info.min):
            raise ValueError(f"dt={dt} <= math.sqrt(sys.float_info.min))={math.sqrt(sys.float_info.min)}")

        mean_vel = (end_pos - start_pos) / dt
        start_accel = (3*mean_vel - (2*start_vel + end_vel))*2/dt
        jerk = (((start_vel + end_vel)/2) - mean_vel)*12/(dt*dt)

        return cls(start_time=start_time,
                   start_pos=start_pos,
                   start_vel=start_vel,
                   start_accel=start_accel,
                   jerk=jerk)

    def limits(self, end_time):
        """Compute limits of motion given an end time.

        Parameters
        ----------
        end_time : `float`
            End time (TAI unix seconds, e.g. from lsst.ts.salobj.curr_tai()).

        Returns
        -------
        limits : `Limits`
            Motion limits between start and end time, inclusive.

        Raises
        ------
        ValueError
            If ``end_time - self.start_time <= math.sqrt(sys.float_info.min)``
            (to avoid overflow).
        """
        dt = end_time - self.start_time
        # Avoid overflow
        if dt <= math.sqrt(sys.float_info.min):
            raise ValueError(f"dt={dt} <= math.sqrt(sys.float_info.min))={math.sqrt(sys.float_info.min)}")

        # Compute maximum |velocity| (max_vel); this may occur
        # at the endpoints or at time t_vex = -start_accel/jerk.
        # Compute t_vex and vex = v(t_vex); if t_vex is not in range [0, dt),
        # set t_vex = 0, so that vex = start_vel
        end_pva = self.pva(end_time)

        if abs(self.start_accel) < abs(self.jerk * dt):
            t_vex = max(-self.start_accel/self.jerk, 0.0)
        else:
            t_vex = 0.0
        vex = self.start_vel + t_vex*(self.start_accel + (t_vex/2)*self.jerk)
        max_vel = max(abs(self.start_vel), abs(end_pva.vel), abs(vex))
        max_accel = max(abs(self.start_accel), abs(end_pva.accel))

        t_pexArr = [0]*2
        numArr = [0]*2
        pexArr = [0]*2

        # Compute the two times t_pexArr,
        # and positions pexArr = p(t_pexArr).
        # If a t_pexArr is out of range [0, dt), set it to 0
        # (so its pexArr = self.start_pos).
        if abs(self.start_vel) < abs(self.start_accel * dt):
            t_pex_zeroj = max(-self.start_vel / self.start_accel, 0.0)
        else:
            t_pex_zeroj = 0.0
        sqrt_arg = (self.start_accel * self.start_accel) - (2.0 * self.start_vel * self.jerk)
        if sqrt_arg < 0.0:
            t_pexArr[0] = 0.0
            t_pexArr[1] = 0.0
        else:
            sqrt_val = math.sqrt(sqrt_arg)
            numArr[0] = -self.start_accel - sqrt_val
            numArr[1] = -self.start_accel + sqrt_val
            for branch in range(2):
                if abs(numArr[branch]) < abs(self.jerk * dt):
                    t_pexArr[branch] = max(0.0, numArr[branch] / self.jerk)
                else:
                    t_pexArr[branch] = t_pex_zeroj
        for branch in range(2):
            t_branch = t_pexArr[branch]
            vel_branch = self.start_vel + (t_branch/2)*(self.start_accel + (t_branch/3)*self.jerk)
            pexArr[branch] = self.start_pos + t_branch*vel_branch

        min_pos = min(self.start_pos, end_pva.pos, pexArr[0], pexArr[1])
        max_pos = max(self.start_pos, end_pva.pos, pexArr[0], pexArr[1])

        return MotionLimits(min_pos=min_pos,
                            max_pos=max_pos,
                            max_vel=max_vel,
                            max_accel=max_accel)

    def pva(self, t):
        """Compute position, velocity and acceleration at a given time.

        Parameters
        ----------
        t : `float`
            Time (TAI unix seconds, e.g. from lsst.ts.salobj.curr_tai()).

        Returns
        -------
        pva : `PosVelAccel`
            Position, velocity and acceleration at the specified time.
        """
        dt = t - self.start_time
        return PosVelAccel(
            pos=self.start_pos + dt*(self.start_vel + dt*(0.5*self.start_accel + dt*self.jerk/6.0)),
            vel=self.start_vel + dt*(self.start_accel + dt*(0.5*self.jerk)),
            accel=self.start_accel + dt*self.jerk,
        )

    def __repr__(self):
        fields = [f"start_time={self.start_time}"]
        for name in ("start_pos", "start_vel", "start_accel", "jerk"):
            val = getattr(self, name)
            if val != 0:
                fields.append(f"{name}={val}")
        return f"PathSegment({', '.join(fields)})"
