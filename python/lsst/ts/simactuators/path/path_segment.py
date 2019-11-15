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

__all__ = ["MotionLimits", "PathSegment"]


class MotionLimits:
    """Limits of motion.

    Parameters
    ----------
    min_position : `float`
        Minimum position (deg)
    max_position : `float`
        Maximum position (deg)
    max_velocity : `float`
        Maximum absolute value of velocity (deg/sec)
    max_acceleration : `float`
        Maximum absolute value of acceleration (deg/sec/sec)
    """
    def __init__(self, min_position, max_position, max_velocity, max_acceleration):
        self.min_position = min_position
        self.max_position = max_position
        self.max_velocity = max_velocity
        self.max_acceleration = max_acceleration


class PathSegment:
    """A segment of a `Path`, a path of constant jerk.

    Parameters
    ----------
    tai : `float`
        TAI time (unix seconds, e.g. from lsst.ts.salobj.curr_tai()).
    pos : `float` at time ``tai`` (optional)
        Position (deg)
    vel : `float` at time ``tai`` (optional)
        Velocity at time ``tai`` (deg/sec)
    accel : `float` (optional)
        Acceleration at time ``tai`` (deg/sec^2)
    jerk : `float` (optional)
        Jerk (deg/sec^3)
    """
    def __init__(self, tai, pos=0, vel=0, accel=0, jerk=0):
        self.tai = float(tai)
        self.pos = float(pos)
        self.vel = float(vel)
        self.accel = float(accel)
        self.jerk = float(jerk)

    @classmethod
    def from_end_conditions(cls, start_tai, start_position, start_velocity,
                            end_tai, end_position, end_velocity):
        """Create a path segment from end conditions.

        Parameters
        ----------
        start_tai : `float`
            Initial time (TAI unix seconds, e.g. from
            lsst.ts.salobj.curr_tai()).
        start_position : `float` (optional)
            Position at ``start_tai`` (deg)
        start_velocity : `float` (optional)
            Velocity at ``start_tai`` (deg/sec)
        end_tai : `float`
            Final time (TAI unix seconds, e.g. from lsst.ts.salobj.curr_tai()).
        end_position : `float` (optional)
            Position at ``end_tai`` (deg)
        end_velocity : `float` (optional)
            Velocity  at ``end_tai`` (deg/sec)

        Raises
        ------
        ValueError
            If end_tai - start_tai <= 0.
        """
        dt = end_tai - start_tai
        # Avoid overflow
        if dt <= math.sqrt(sys.float_info.min):
            raise ValueError(f"dt={dt} <= math.sqrt(sys.float_info.min))={math.sqrt(sys.float_info.min)}")

        mean_vel = (end_position - start_position) / dt
        start_accel = (3*mean_vel - (2*start_velocity + end_velocity))*2/dt
        jerk = (((start_velocity + end_velocity)/2) - mean_vel)*12/(dt*dt)

        return cls(tai=start_tai,
                   pos=start_position,
                   vel=start_velocity,
                   accel=start_accel,
                   jerk=jerk)

    def limits(self, end_tai):
        """Compute limits of motion given an end time.

        Parameters
        ----------
        end_tai : `float`
            End time (TAI unix seconds, e.g. from lsst.ts.salobj.curr_tai()).

        Returns
        -------
        limits : `Limits`
            Motion limits between start and end time, inclusive.

        Raises
        ------
        ValueError
            If ``end_tai - self.tai <= math.sqrt(sys.float_info.min)``
            (to avoid overflow).
        """
        dt = end_tai - self.tai
        # Avoid overflow
        if dt <= math.sqrt(sys.float_info.min):
            raise ValueError(f"dt={dt} <= math.sqrt(sys.float_info.min))={math.sqrt(sys.float_info.min)}")

        # Compute maximum |velocity| (max_velocity); this may occur
        # at the endpoints or at time t_vex = -start_accel/jerk.
        # Compute t_vex and vex = v(t_vex); if t_vex is not in range [0, dt),
        # set t_vex = 0, so that vex = start_velocity
        end_pva = self.at(end_tai)

        if abs(self.accel) < abs(self.jerk * dt):
            t_vex = max(-self.accel/self.jerk, 0.0)
        else:
            t_vex = 0.0
        vex = self.vel + t_vex*(self.accel + (t_vex/2)*self.jerk)
        max_velocity = max(abs(self.vel), abs(end_pva.vel), abs(vex))
        max_acceleration = max(abs(self.accel), abs(end_pva.accel))

        t_pexArr = [0]*2
        numArr = [0]*2
        pexArr = [0]*2

        # Compute the two times t_pexArr,
        # and positions pexArr = p(t_pexArr).
        # If a t_pexArr is out of range [0, dt), set it to 0
        # (so its pexArr = self.pos).
        if abs(self.vel) < abs(self.accel * dt):
            t_pex_zeroj = max(-self.vel / self.accel, 0.0)
        else:
            t_pex_zeroj = 0.0
        sqrt_arg = (self.accel * self.accel) - (2.0 * self.vel * self.jerk)
        if sqrt_arg < 0.0:
            t_pexArr[0] = 0.0
            t_pexArr[1] = 0.0
        else:
            sqrt_val = math.sqrt(sqrt_arg)
            numArr[0] = -self.accel - sqrt_val
            numArr[1] = -self.accel + sqrt_val
            for branch in range(2):
                if abs(numArr[branch]) < abs(self.jerk * dt):
                    t_pexArr[branch] = max(0.0, numArr[branch] / self.jerk)
                else:
                    t_pexArr[branch] = t_pex_zeroj
        for branch in range(2):
            t_branch = t_pexArr[branch]
            vel_branch = self.vel + (t_branch/2)*(self.accel + (t_branch/3)*self.jerk)
            pexArr[branch] = self.pos + t_branch*vel_branch

        min_position = min(self.pos, end_pva.pos, pexArr[0], pexArr[1])
        max_position = max(self.pos, end_pva.pos, pexArr[0], pexArr[1])

        return MotionLimits(min_position=min_position,
                            max_position=max_position,
                            max_velocity=max_velocity,
                            max_acceleration=max_acceleration)

    def at(self, tai):
        """Return a copy with the specified time.

        Parameters
        ----------
        tai : `float`
            Time (TAI unix seconds, e.g. from lsst.ts.salobj.curr_tai()).

        Returns
        -------
        copy : `PathSegment`
            Copy at the specified time.
        """
        dt = tai - self.tai
        return PathSegment(
            tai=tai,
            pos=self.pos + dt*(self.vel + dt*(0.5*self.accel + dt*self.jerk/6.0)),
            vel=self.vel + dt*(self.accel + dt*(0.5*self.jerk)),
            accel=self.accel + dt*self.jerk,
            jerk=self.jerk,
        )

    def __repr__(self):
        fields = [f"tai={self.tai}"]
        for name in ("pos", "vel", "accel", "jerk"):
            val = getattr(self, name)
            if val != 0:
                fields.append(f"{name}={val}")
        return f"PathSegment({', '.join(fields)})"
