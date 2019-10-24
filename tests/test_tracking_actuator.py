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

import itertools
import math
import unittest

import numpy as np

from lsst.ts import simactuators


class SinFunctor:
    """Functor to compute sine wave position and associated velocity

    pos = pos_off + vel_off*(t-t0) + pos_ampl*sin(theta)
    vel = vel_off + vel_amp*sin(theta)
    where theta = 2 pi (t-t0) / period
    and vel_ampl = 2 pi ampl / period
    """
    def __init__(self, pos_off, pos_ampl, vel_off, period, t0):
        assert period > 0
        self.pos_off = pos_off
        self.pos_ampl = pos_ampl
        self.vel_off = vel_off
        self.period = period
        self.t0 = t0
        self.vel_ampl = self.pos_ampl * 2 * math.pi / self.period

    def __call__(self, t):
        """Compute pos, vel at the specified time."""
        dt = t - self.t0
        theta = 2*math.pi*dt/self.period
        return (
            self.pos_off + self.vel_off*dt + self.pos_ampl*math.sin(theta),
            self.vel_off + self.vel_ampl*math.cos(theta),
        )


class TestTrackingActuator(unittest.TestCase):

    def test_constructor(self):
        min_pos = -1
        max_pos = 2
        max_vel = 3
        max_accel = 4
        dtmax_track = 0.5
        nsettle = 1
        t0 = 0.5

        actuator = simactuators.TrackingActuator(
            min_pos=min_pos, max_pos=max_pos, max_vel=max_vel, max_accel=max_accel,
            dtmax_track=dtmax_track, nsettle=nsettle, t0=t0)
        self.assertEqual(actuator.min_pos, min_pos)
        self.assertEqual(actuator.max_pos, max_pos)
        self.assertEqual(actuator.max_vel, max_vel)
        self.assertEqual(actuator.max_accel, max_accel)
        self.assertEqual(actuator.dtmax_track, dtmax_track)
        self.assertEqual(actuator.nsettle, nsettle)
        self.assertEqual(actuator.cmd.t0, t0)
        self.assertEqual(actuator.cmd.pos0, 0)
        self.assertEqual(actuator.cmd.vel0, 0)
        self.assertEqual(actuator.cmd.accel0, 0)
        self.assertEqual(actuator.cmd.jerk, 0)
        self.assertEqual(actuator.curr[0].t0, t0)
        self.assertEqual(actuator.curr[0].pos0, 0)
        self.assertEqual(actuator.curr[0].vel0, 0)
        self.assertEqual(actuator.curr[0].accel0, 0)
        self.assertEqual(actuator.curr[0].jerk, 0)
        self.assertEqual(len(actuator.curr), 1)
        self.assertEqual(actuator.curr.kind, actuator.Kind.Stopped)

        # pos0 is 0 if in range [min_pos, max_pos) else min_pos
        for min_pos, max_pos in itertools.product((-10, 0, 10), (-9, -1, 9)):
            if max_pos <= min_pos:
                continue
            if min_pos <= 0 <= max_pos:
                expected_p0 = 0
            else:
                expected_p0 = min_pos
            actuator = simactuators.TrackingActuator(min_pos=min_pos, max_pos=max_pos,
                                                     max_vel=max_vel, max_accel=max_accel,
                                                     dtmax_track=dtmax_track, nsettle=nsettle, t0=t0)
            self.assertAlmostEqual(actuator.curr[-1].pos0, expected_p0)
            self.assertAlmostEqual(actuator.cmd.pos0, expected_p0)

    def test_constructor_errors(self):
        min_pos = -1
        max_pos = 2
        max_vel = 3
        max_accel = 4
        dtmax_track = 0.5
        nsettle = 1
        t0 = 0.5

        # require min_pos < max_pos
        with self.assertRaises(ValueError):
            simactuators.TrackingActuator(
                min_pos=max_pos, max_pos=max_pos, max_vel=max_vel, max_accel=max_accel,
                dtmax_track=dtmax_track, nsettle=nsettle, t0=t0)
        with self.assertRaises(ValueError):
            simactuators.TrackingActuator(
                min_pos=max_pos + 1, max_pos=max_pos, max_vel=max_vel, max_accel=max_accel,
                dtmax_track=dtmax_track, nsettle=nsettle, t0=t0)

        # require max_vel > 0
        with self.assertRaises(ValueError):
            simactuators.TrackingActuator(
                min_pos=min_pos, max_pos=max_pos, max_vel=0, max_accel=max_accel,
                dtmax_track=dtmax_track, nsettle=nsettle, t0=t0)
        with self.assertRaises(ValueError):
            simactuators.TrackingActuator(
                min_pos=min_pos, max_pos=max_pos, max_vel=-1, max_accel=max_accel,
                dtmax_track=dtmax_track, nsettle=nsettle, t0=t0)

        # require max_accel > 0
        with self.assertRaises(ValueError):
            simactuators.TrackingActuator(
                min_pos=min_pos, max_pos=max_pos, max_vel=max_vel, max_accel=0,
                dtmax_track=dtmax_track, nsettle=nsettle, t0=t0)
        with self.assertRaises(ValueError):
            simactuators.TrackingActuator(
                min_pos=min_pos, max_pos=max_pos, max_vel=max_vel, max_accel=-1,
                dtmax_track=dtmax_track, nsettle=nsettle, t0=t0)

    def test_matched_start(self):
        """Test slew_cmd for a sine path where initial position and velocity
        for curr and cmd match.

        Follow a sine wave path for one period. There should be only tracking.
        """
        print("\ntest_matched_start")
        period = 40  # period of sine wave (sec)
        min_pos = -100  # large enough to not be relevant
        max_pos = 100  # ditto
        max_vel = 7  # maximum allowed velocity (deg/sec)
        max_accel = 10  # maximum allowed acceleration (deg/sec^2)
        cmd_interval = 0.05  # interval between commanded points (sec)
        for pos_ampl, vel_off, frac_phase in itertools.product(
            (0.1, 0.3), (0, 0.1, -0.2), (0, 0.2, 0.7),
        ):
            with self.subTest(pos_ampl=pos_ampl, vel_off=vel_off, frac_phase=frac_phase):
                self.check_sin_path(
                    pos_off=0, pos_ampl=pos_ampl, vel_off=0,
                    period=period, frac_phase=frac_phase,
                    cmd_interval=cmd_interval,
                    min_pos=min_pos, max_pos=max_pos, max_vel=max_vel, max_accel=max_accel,
                    nsettle=1,
                    max_nslew=0,
                    max_pos_err=0.00001, max_vel_err=0.001,
                )

    def test_mismatched_start(self):
        """Test slew_cmd for a sine path where initial position and velocity
        for curr and cmd do not match.

        Follow a sine wave path for one period. There should be an initial
        slew followed by only tracking.
        """
        print("\ntest_mismatched_start")
        period = 40
        min_pos = -100  # large enough to not be relevant
        max_pos = 100  # ditto
        max_vel = 5
        max_accel = 5
        cmd_interval = 0.05
        for pos_off, pos_ampl, vel_off, frac_phase in itertools.product(
            (1, -30, 30), (0.1, 0.3), (0, -0.1, 0.2), (0, 0.2, 0.7),
        ):
            with self.subTest(pos_off=pos_off, pos_ampl=pos_ampl, vel_off=vel_off, frac_phase=frac_phase):
                self.check_sin_path(
                    pos_off=pos_off, pos_ampl=pos_ampl, vel_off=vel_off,
                    period=period, frac_phase=frac_phase,
                    cmd_interval=cmd_interval,
                    min_pos=min_pos, max_pos=max_pos, max_vel=max_vel, max_accel=max_accel,
                    nsettle=1,
                    max_nslew=200,
                    max_pos_err=0.00001, max_vel_err=0.001,
                )

    def test_stop(self):
        for pos_off, vel_off in itertools.product(
            (30, -30), (0, 1, -1),
        ):
            cmd_interval = 0.05
            with self.subTest(pos_off=pos_off, vel_off=vel_off):
                actuator = simactuators.TrackingActuator(
                    min_pos=-100, max_pos=100,
                    max_vel=5, max_accel=10,
                    dtmax_track=0.05, nsettle=1, t0=0)
                self.assertEqual(actuator.kind(0), actuator.Kind.Stopped)
                for i in range(5):
                    t = i*cmd_interval + 1  # +1 to avoid 0
                    cmd_pos = pos_off + t*vel_off
                    actuator.set_cmd(pos=cmd_pos, vel=vel_off, t=t)
                self.assertEqual(actuator.kind(t), actuator.Kind.Slewing)

                # expected starting position and velocity for the halt
                t_start_halt = t + cmd_interval
                p_start_halt, v_start_halt = actuator.curr.pva(t_start_halt)[0:2]
                actuator.stop(t=t_start_halt)
                self.assertEqual(actuator.kind(t_start_halt), actuator.Kind.Stopping)
                pcurr_t, vcurr_t = actuator.curr.pva(t_start_halt)[0:2]
                self.assertAlmostEqual(pcurr_t, p_start_halt)
                self.assertAlmostEqual(vcurr_t, v_start_halt)

                # the ending velocity and acceleration are, of course, 0
                # and the curr and cmd positions should match at the end
                t_end_halt = actuator.curr[-1].t0
                p_end_halt, v_end_halt, a_end_halt = actuator.curr.pva(t_end_halt)
                self.assertEqual(v_end_halt, 0)
                self.assertEqual(a_end_halt, 0)
                pcmd, vcmd, acmd = actuator.cmd.pva(t)
                self.assertAlmostEqual(pcmd, p_end_halt)
                self.assertEqual(vcmd, 0)
                self.assertEqual(acmd, 0)
                self.assertGreater(t_end_halt, t_start_halt)
                # the end kind should be Stopped;
                # add a margin to the time  to avoid roundoff error
                self.assertEqual(actuator.kind(t_end_halt + 0.0001), actuator.Kind.Stopped)

    def test_abort(self):
        for pos_off, vel_off in itertools.product(
            (30, -30), (0, 1, -1),
        ):
            cmd_interval = 0.05
            with self.subTest(pos_off=pos_off, vel_off=vel_off):
                actuator = simactuators.TrackingActuator(
                    min_pos=-100, max_pos=100,
                    max_vel=5, max_accel=10,
                    dtmax_track=0.05, nsettle=1, t0=0)
                self.assertEqual(actuator.kind(0), actuator.Kind.Stopped)
                for i in range(5):
                    t = i*cmd_interval + 1  # +1 to avoid 0
                    cmd_pos = pos_off + t*vel_off
                    actuator.set_cmd(pos=cmd_pos, vel=vel_off, t=t)
                self.assertEqual(actuator.kind(t), actuator.Kind.Slewing)

                # expected starting position and velocity for the halt
                t_abort = t + cmd_interval
                p_start_halt, v_start_halt = actuator.curr.pva(t_abort)[0:2]
                actuator.abort(t=t_abort)
                self.assertEqual(actuator.kind(t_abort), actuator.Kind.Stopped)
                pcurr_t, vcurr_t = actuator.curr.pva(t_abort)[0:2]
                self.assertAlmostEqual(pcurr_t, p_start_halt)
                self.assertAlmostEqual(vcurr_t, 0)
                self.assertEqual(len(actuator.curr), 1)

                # abort does not change actuator.cmd
                cmd_pos, cmd_vel = actuator.cmd.pva(t_abort)[0:2]
                desired_cmd_pos = pos_off + t_abort*vel_off
                self.assertAlmostEqual(cmd_pos, desired_cmd_pos)
                self.assertAlmostEqual(cmd_vel, vel_off)

    def check_sin_path(self, pos_off, pos_ampl, vel_off, period, frac_phase, cmd_interval,
                       min_pos, max_pos, max_vel, max_accel, nsettle, max_nslew,
                       max_pos_err, max_vel_err):
        """Check slewing and tracking to a sinusoidal path.

        Parameters
        ----------
        pos_off : `float`
            Position offset between start and center of sine wave (deg).
            The difference is split evenly between current and commanded
            (start and end).
        pos_ampl : `float`
            Amplitude of sine wave (deg)
        vel_off : `float`
            Velocity offset between start and center of sine wave (deg/sec).
            The difference is split evenly between current and commanded
            (start and end).
        period : `float`
            Period of sine wave (sec)
        frac_phase : `float`
            Phase to start sine wave, as a fraction of period
        cmd_interval : `float`
            Interval at which to call `set_cmd` (sec)
        min_pos : `float`
            Minimum allowed position (deg)
        max_pos : `float`
            Maximum allowed position (deg)
        max_vel : `float`
            Maximum allowed velocity (deg/sec)
        max_accel : `float`
            Maximum allowed acceleration (deg/sec^2)
        nsettle : `int`
            Number of consective tracking updates before
            actuator.kind(t) reports tracking.
        max_nslew : `int`
            Maximum allowed number of slew iterations.
        max_pos_err : `float`
            Maximum allowed position error (deg).
            This is checked at times within +/- one interval
            of each call to `TrackingActuator.set_cmd` while
            ``TrackingActuator.kind(t)`` is tracking.
        max_vel_err : `float`
            Maximum allowed velocity error (deg).
            This is checked at times within +/- one interval
            of each call to `TrackingActuator.set_cmd` while
            ``TrackingActuator.kind(t)`` is tracking.
        """
        ncmd = int(period/cmd_interval)
        dtmax_track = cmd_interval*2  # > cmd_interval to allow tracking

        sinfunc = SinFunctor(
            pos_off=pos_off/2,
            pos_ampl=pos_ampl,
            vel_off=vel_off/2,
            period=period,
            t0=frac_phase*period,
        )
        pos0, vel0 = sinfunc(0)
        actuator = simactuators.TrackingActuator(
            min_pos=min_pos, max_pos=max_pos,
            max_vel=max_vel, max_accel=max_accel,
            dtmax_track=dtmax_track, nsettle=nsettle,
            t0=-cmd_interval)  # so we can start with tracking, of possible
        self.assertEqual(len(actuator.curr), 1)
        actuator.curr[0].pos0 = pos0 - pos_off
        actuator.curr[0].vel0 = vel0 - vel_off
        self.assertEqual(actuator.curr[0].accel0, 0)
        self.assertEqual(actuator.curr.kind, actuator.Kind.Stopped)
        actuator.curr[0].t0 = 0
        pos_errors = []
        vel_errors = []
        curr_vels = []
        curr_accels = []
        # count of number o times actuator.curr.kind is slew or track
        nslew = 0
        ntrack = 0

        for cmd_t in np.linspace(start=0, stop=period, num=ncmd):
            pos, vel = sinfunc(cmd_t)
            actuator.set_cmd(pos=pos, vel=vel, t=cmd_t)

            # Check that actuator.kind(t) transitions
            # from Slewing to Tracking after nsettle instances of
            # actuator.curr.kind being Tracking
            if actuator.curr.kind == actuator.Kind.Tracking:
                ntrack += 1
                self.assertEqual(actuator._ntrack, ntrack)
                if actuator._ntrack > actuator.nsettle:
                    self.assertEqual(actuator.kind(cmd_t), actuator.Kind.Tracking)
                else:
                    self.assertEqual(actuator.kind(cmd_t), actuator.Kind.Slewing)
            elif actuator.curr.kind == actuator.Kind.Slewing:
                nslew += 1
                if ntrack > 0:
                    self.fail(f"Slew found after tracking: nslew={nslew}; ntrack={ntrack}")

            # check that the commanded PVT was properly recorded
            cmd_pos, cmd_vel, cmd_accel = actuator.cmd.pva(cmd_t)
            self.assertAlmostEqual(cmd_pos, pos)
            self.assertAlmostEqual(cmd_vel, vel)
            self.assertAlmostEqual(cmd_accel, 0)

            # Check tracking error, once we have settled
            # (when actuator.kind(t) says we are tracking)
            if actuator.kind(cmd_t) == actuator.Kind.Tracking:
                for frac_dt in (-0.6, -0.3, 0, 0.3, 0.6):
                    test_t = cmd_interval*frac_dt + cmd_t
                    cmd_pos, cmd_vel, cmd_accel = actuator.cmd.pva(test_t)
                    curr_pos, curr_vel, curr_accel = actuator.curr.pva(test_t)
                    pos_errors.append(curr_pos - cmd_pos)
                    vel_errors.append(curr_vel - cmd_vel)
                    curr_vels.append(curr_vel)
                    curr_accels.append(curr_accel)
        pos_err = np.abs(pos_errors).max()
        vel_err = np.abs(vel_errors).max()
        max_vel = np.abs(curr_vels).max()
        max_accel = np.abs(curr_accels).max()
        print(f"pos_err={pos_err:0.1e}; vel_err={vel_err:0.1e}; nslew={nslew}")
        self.assertLessEqual(nslew, max_nslew)
        self.assertLess(pos_err, max_pos_err)
        self.assertLess(vel_err, max_vel_err)
        self.assertLessEqual(max_vel, max_vel)
        # use a fudge factor for accel because it will typically
        # be at the limit
        self.assertLessEqual(max_accel, max_accel*1.000001)


if __name__ == '__main__':
    unittest.main()
