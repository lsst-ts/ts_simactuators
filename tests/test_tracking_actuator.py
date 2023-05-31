# This file is part of ts_simactuators.
#
# Developed for the Rubin Observatory Telescope and Site System.
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
import typing
import unittest

import numpy as np
import pytest
from lsst.ts import simactuators


class SinFunctor:
    """Functor to compute sine wave position and associated velocity

    position = pos_off + vel_off*(tai-start_tai) + pos_ampl*sin(theta)
    velocity = vel_off + vel_amp*sin(theta)
    where theta = 2 pi (tai-start_tai) / period
    and vel_ampl = 2 pi ampl / period
    """

    def __init__(
        self,
        pos_off: float,
        pos_ampl: float,
        vel_off: float,
        period: float,
        start_tai: float,
    ) -> None:
        assert period > 0
        self.pos_off = pos_off
        self.pos_ampl = pos_ampl
        self.vel_off = vel_off
        self.period = period
        self.start_tai = start_tai
        self.vel_ampl = self.pos_ampl * 2 * math.pi / self.period

    def __call__(self, tai: float) -> typing.Tuple[float, float]:
        """Compute position, velocity at the specified time."""
        dt = tai - self.start_tai
        theta = 2 * math.pi * dt / self.period
        return (
            self.pos_off + self.vel_off * dt + self.pos_ampl * math.sin(theta),
            self.vel_off + self.vel_ampl * math.cos(theta),
        )


class TestTrackingActuator(unittest.TestCase):
    def test_constructor(self) -> None:
        min_position = -1
        max_position = 2
        max_velocity = 3
        max_acceleration = 4
        dtmax_track = 0.5
        nsettle = 1
        tai = 0.5

        actuator = simactuators.TrackingActuator(
            min_position=min_position,
            max_position=max_position,
            max_velocity=max_velocity,
            max_acceleration=max_acceleration,
            dtmax_track=dtmax_track,
            nsettle=nsettle,
            tai=tai,
        )
        assert actuator.min_position == min_position
        assert actuator.max_position == max_position
        assert actuator.max_velocity == max_velocity
        assert actuator.max_acceleration == max_acceleration
        assert actuator.dtmax_track == dtmax_track
        assert actuator.nsettle == nsettle
        assert actuator.target.tai == tai
        assert actuator.target.position == 0
        assert actuator.target.velocity == 0
        assert actuator.target.acceleration == 0
        assert actuator.target.jerk == 0
        assert actuator.path[0].tai == tai
        assert actuator.path[0].position == 0
        assert actuator.path[0].velocity == 0
        assert actuator.path[0].acceleration == 0
        assert actuator.path[0].jerk == 0
        assert len(actuator.path) == 1
        assert actuator.path.kind == actuator.Kind.Stopped

        for start_position in (
            min_position,
            (min_position + max_position) / 2,
            max_position,
        ):
            actuator = simactuators.TrackingActuator(
                min_position=min_position,
                max_position=max_position,
                max_velocity=max_velocity,
                max_acceleration=max_acceleration,
                dtmax_track=dtmax_track,
                nsettle=nsettle,
                tai=tai,
                start_position=start_position,
            )
            assert actuator.path[-1].position == pytest.approx(start_position)
            assert actuator.target.position == pytest.approx(start_position)

        # Check that initial position is 0 if in range
        # [min_position, max_position), else min_position.
        for min_position, max_position in itertools.product(
            (-10, 0, 10),
            (-9, -1, 9),
        ):
            if max_position <= min_position:
                continue
            if min_position <= 0 <= max_position:
                expected_p0 = 0
            else:
                expected_p0 = min_position
            actuator = simactuators.TrackingActuator(
                min_position=min_position,
                max_position=max_position,
                max_velocity=max_velocity,
                max_acceleration=max_acceleration,
                dtmax_track=dtmax_track,
                nsettle=nsettle,
                tai=tai,
            )
            assert actuator.path[-1].position == pytest.approx(expected_p0)
            assert actuator.target.position == pytest.approx(expected_p0)

    def test_constructor_errors(self) -> None:
        min_position = -1
        max_position = 2
        max_velocity = 3
        max_acceleration = 4
        dtmax_track = 0.5
        nsettle = 1
        tai = 0.5

        # require min_position < max_position
        with pytest.raises(ValueError):
            simactuators.TrackingActuator(
                min_position=max_position,
                max_position=max_position,
                max_velocity=max_velocity,
                max_acceleration=max_acceleration,
                dtmax_track=dtmax_track,
                nsettle=nsettle,
                tai=tai,
            )
        with pytest.raises(ValueError):
            simactuators.TrackingActuator(
                min_position=max_position + 1,
                max_position=max_position,
                max_velocity=max_velocity,
                max_acceleration=max_acceleration,
                dtmax_track=dtmax_track,
                nsettle=nsettle,
                tai=tai,
            )

        # require max_velocity > 0
        with pytest.raises(ValueError):
            simactuators.TrackingActuator(
                min_position=min_position,
                max_position=max_position,
                max_velocity=0,
                max_acceleration=max_acceleration,
                dtmax_track=dtmax_track,
                nsettle=nsettle,
                tai=tai,
            )
        with pytest.raises(ValueError):
            simactuators.TrackingActuator(
                min_position=min_position,
                max_position=max_position,
                max_velocity=-1,
                max_acceleration=max_acceleration,
                dtmax_track=dtmax_track,
                nsettle=nsettle,
                tai=tai,
            )

        # require max_acceleration > 0
        with pytest.raises(ValueError):
            simactuators.TrackingActuator(
                min_position=min_position,
                max_position=max_position,
                max_velocity=max_velocity,
                max_acceleration=0,
                dtmax_track=dtmax_track,
                nsettle=nsettle,
                tai=tai,
            )
        with pytest.raises(ValueError):
            simactuators.TrackingActuator(
                min_position=min_position,
                max_position=max_position,
                max_velocity=max_velocity,
                max_acceleration=-1,
                dtmax_track=dtmax_track,
                nsettle=nsettle,
                tai=tai,
            )

    def test_matched_start(self) -> None:
        """Test slew_cmd for a sine path where initial position and velocity
        for current and target match.

        Follow a sine wave path for one period. There should be only tracking.
        """
        print("\ntest_matched_start")
        period = 40  # period of sine wave (sec)
        min_position = -100  # large enough to not be relevant
        max_position = 100  # ditto
        max_velocity = 7  # maximum allowed velocity (deg/sec)
        max_acceleration = 10  # maximum allowed acceleration (deg/sec^2)
        cmd_interval = 0.05  # interval between commanded points (sec)
        for pos_ampl, vel_off, frac_phase in itertools.product(
            (0.1, 0.3),
            (0, 0.1, -0.2),
            (0, 0.2, 0.7),
        ):
            with self.subTest(
                pos_ampl=pos_ampl, vel_off=vel_off, frac_phase=frac_phase
            ):
                self.check_sin_path(
                    pos_off=0,
                    pos_ampl=pos_ampl,
                    vel_off=0,
                    period=period,
                    frac_phase=frac_phase,
                    cmd_interval=cmd_interval,
                    min_position=min_position,
                    max_position=max_position,
                    max_velocity=max_velocity,
                    max_acceleration=max_acceleration,
                    nsettle=1,
                    max_nslew=0,
                    max_position_err=0.00001,
                    max_velocity_err=0.001,
                )

    def test_mismatched_start(self) -> None:
        """Test slew_cmd for a sine path where initial position and velocity
        for current and target do not match.

        Follow a sine wave path for one period. There should be an initial
        slew followed by only tracking.
        """
        print("\ntest_mismatched_start")
        period = 40
        min_position = -100  # large enough to not be relevant
        max_position = 100  # ditto
        max_velocity = 5
        max_acceleration = 5
        cmd_interval = 0.05
        for pos_off, pos_ampl, vel_off, frac_phase in itertools.product(
            (1, -30, 30),
            (0.1, 0.3),
            (0, -0.1, 0.2),
            (0, 0.2, 0.7),
        ):
            with self.subTest(
                pos_off=pos_off,
                pos_ampl=pos_ampl,
                vel_off=vel_off,
                frac_phase=frac_phase,
            ):
                self.check_sin_path(
                    pos_off=pos_off,
                    pos_ampl=pos_ampl,
                    vel_off=vel_off,
                    period=period,
                    frac_phase=frac_phase,
                    cmd_interval=cmd_interval,
                    min_position=min_position,
                    max_position=max_position,
                    max_velocity=max_velocity,
                    max_acceleration=max_acceleration,
                    nsettle=1,
                    max_nslew=200,
                    max_position_err=0.00001,
                    max_velocity_err=0.001,
                )

    def test_stop(self) -> None:
        for pos_off, vel_off in itertools.product(
            (30, -30),
            (0, 1, -1),
        ):
            cmd_interval = 0.05
            with self.subTest(pos_off=pos_off, vel_off=vel_off):
                actuator = simactuators.TrackingActuator(
                    min_position=-100,
                    max_position=100,
                    max_velocity=5,
                    max_acceleration=10,
                    dtmax_track=0.05,
                    nsettle=1,
                    tai=0,
                )
                assert actuator.kind(0) == actuator.Kind.Stopped
                for i in range(5):
                    tai = i * cmd_interval + 1  # +1 to avoid 0
                    cmd_pos = pos_off + tai * vel_off
                    actuator.set_target(position=cmd_pos, velocity=vel_off, tai=tai)
                assert actuator.kind(tai) == actuator.Kind.Slewing

                # expected starting position and velocity for the halt
                tai_start_halt = tai + cmd_interval
                segment_before_stop = actuator.path.at(tai_start_halt)
                actuator.stop(tai=tai_start_halt)
                assert actuator.kind(tai_start_halt) == actuator.Kind.Stopping
                segment_start_halt = actuator.path.at(tai_start_halt)
                assert segment_before_stop.position == pytest.approx(
                    segment_start_halt.position
                )
                assert segment_before_stop.velocity == pytest.approx(
                    segment_start_halt.velocity
                )

                # the ending velocity and acceleration are, of course, 0
                # and the current and target positions should match at the end
                tai_end_halt = actuator.path[-1].tai
                segment_end_halt = actuator.path.at(tai_end_halt)
                assert segment_end_halt.velocity == 0
                assert segment_end_halt.acceleration == 0
                cmd_segment_end_halt = actuator.target.at(tai=tai_end_halt)
                assert cmd_segment_end_halt.position == pytest.approx(
                    segment_end_halt.position
                )
                assert cmd_segment_end_halt.velocity == 0
                assert cmd_segment_end_halt.acceleration == 0
                assert tai_end_halt > tai_start_halt
                # the end kind should be Stopped;
                # add a margin to the time  to avoid roundoff error
                assert actuator.kind(tai_end_halt + 0.0001) == actuator.Kind.Stopped

    def test_abort(self) -> None:
        for pos_off, vel_off in itertools.product(
            (30, -30),
            (0, 1, -1),
        ):
            cmd_interval = 0.05
            with self.subTest(pos_off=pos_off, vel_off=vel_off):
                actuator = simactuators.TrackingActuator(
                    min_position=-100,
                    max_position=100,
                    max_velocity=5,
                    max_acceleration=10,
                    dtmax_track=0.05,
                    nsettle=1,
                    tai=0,
                )
                assert actuator.kind(0) == actuator.Kind.Stopped
                for i in range(5):
                    tai = i * cmd_interval + 1  # +1 to avoid 0
                    cmd_pos = pos_off + tai * vel_off
                    actuator.set_target(position=cmd_pos, velocity=vel_off, tai=tai)
                assert actuator.kind(tai) == actuator.Kind.Slewing

                # expected starting position and velocity for the halt
                t_abort = tai + cmd_interval
                segment_before_abort = actuator.path.at(t_abort)
                actuator.abort(tai=t_abort)
                assert actuator.kind(t_abort) == actuator.Kind.Stopped
                segment_abort = actuator.path.at(t_abort)
                assert segment_before_abort.position == pytest.approx(
                    segment_abort.position
                )
                assert segment_abort.velocity == pytest.approx(0)
                assert len(actuator.path) == 1

                # abort does not change actuator.target
                cmd_segment_abort = actuator.target.at(t_abort)
                desired_cmd_pos = pos_off + t_abort * vel_off
                assert cmd_segment_abort.position == pytest.approx(desired_cmd_pos)
                assert cmd_segment_abort.velocity == pytest.approx(vel_off)

    def check_sin_path(
        self,
        pos_off: float,
        pos_ampl: float,
        vel_off: float,
        period: float,
        frac_phase: float,
        cmd_interval: float,
        min_position: float,
        max_position: float,
        max_velocity: float,
        max_acceleration: float,
        nsettle: int,
        max_nslew: int,
        max_position_err: float,
        max_velocity_err: float,
    ) -> None:
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
            Interval at which to call `set_target` (sec)
        min_position : `float`
            Minimum allowed position (deg)
        max_position : `float`
            Maximum allowed position (deg)
        max_velocity : `float`
            Maximum allowed velocity (deg/sec)
        max_acceleration : `float`
            Maximum allowed acceleration (deg/sec^2)
        nsettle : `int`
            Number of consective tracking updates before
            actuator.kind(tai) reports tracking.
        max_nslew : `int`
            Maximum allowed number of slew iterations.
        max_position_err : `float`
            Maximum allowed position error (deg).
            This is checked at times within +/- one interval
            of each call to `TrackingActuator.set_target` while
            ``TrackingActuator.kind(tai)`` is tracking.
        max_velocity_err : `float`
            Maximum allowed velocity error (deg).
            This is checked at times within +/- one interval
            of each call to `TrackingActuator.set_target` while
            ``TrackingActuator.kind(tai)`` is tracking.
        """
        ncmd = int(period / cmd_interval)
        dtmax_track = cmd_interval * 2  # > cmd_interval to allow tracking

        sinfunc = SinFunctor(
            pos_off=pos_off / 2,
            pos_ampl=pos_ampl,
            vel_off=vel_off / 2,
            period=period,
            start_tai=frac_phase * period,
        )
        position, velocity = sinfunc(0)
        actuator = simactuators.TrackingActuator(
            min_position=min_position,
            max_position=max_position,
            max_velocity=max_velocity,
            max_acceleration=max_acceleration,
            dtmax_track=dtmax_track,
            nsettle=nsettle,
            tai=-cmd_interval,
        )  # so we can start with tracking, of possible
        assert len(actuator.path) == 1
        actuator.path[0].position = position - pos_off
        actuator.path[0].velocity = velocity - vel_off
        assert actuator.path[0].acceleration == 0
        assert actuator.path.kind == actuator.Kind.Stopped
        actuator.path[0].tai = 0
        pos_errors = []
        vel_errors = []
        curr_vels = []
        curr_accels = []
        # count of number o times actuator.path.kind is slew or track
        nslew = 0
        ntrack = 0

        for cmd_t in np.linspace(start=0, stop=period, num=ncmd):
            position, velocity = sinfunc(cmd_t)
            actuator.set_target(position=position, velocity=velocity, tai=cmd_t)

            # Check that actuator.kind(tai) transitions
            # from Slewing to Tracking after nsettle instances of
            # actuator.path.kind being Tracking
            if actuator.path.kind == actuator.Kind.Tracking:
                ntrack += 1
                assert actuator._ntrack == ntrack
                if actuator._ntrack > actuator.nsettle:
                    assert actuator.kind(cmd_t) == actuator.Kind.Tracking
                else:
                    assert actuator.kind(cmd_t) == actuator.Kind.Slewing
            elif actuator.path.kind == actuator.Kind.Slewing:
                nslew += 1
                if ntrack > 0:
                    self.fail(
                        f"Slew found after tracking: nslew={nslew}; ntrack={ntrack}"
                    )

            # check that the commanded PVT was properly recorded
            target_segment = actuator.target.at(cmd_t)
            assert target_segment.position == pytest.approx(position)
            assert target_segment.velocity == pytest.approx(velocity)
            assert target_segment.acceleration == pytest.approx(0)

            # Check tracking error, once we have settled
            # (when actuator.kind(tai) says we are tracking)
            if actuator.kind(cmd_t) == actuator.Kind.Tracking:
                for frac_dt in (-0.6, -0.3, 0, 0.3, 0.6):
                    test_t = cmd_interval * frac_dt + cmd_t
                    target_segment = actuator.target.at(test_t)
                    current_segment = actuator.path.at(test_t)
                    pos_errors.append(
                        current_segment.position - target_segment.position
                    )
                    vel_errors.append(
                        current_segment.velocity - target_segment.velocity
                    )
                    curr_vels.append(current_segment.velocity)
                    curr_accels.append(current_segment.acceleration)
        pos_err = np.abs(pos_errors).max()
        vel_err = np.abs(vel_errors).max()
        max_velocity = np.abs(curr_vels).max()
        max_acceleration = np.abs(curr_accels).max()
        print(f"pos_err={pos_err:0.1e}; vel_err={vel_err:0.1e}; nslew={nslew}")
        assert nslew <= max_nslew
        assert pos_err < max_position_err
        assert vel_err < max_velocity_err
        assert max_velocity <= max_velocity
        # use a fudge factor for acceleration because it will typically
        # be at the limit
        assert max_acceleration <= max_acceleration * 1.000001
