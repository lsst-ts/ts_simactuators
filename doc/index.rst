.. py:currentmodule:: lsst.ts.simactuators

.. _lsst.ts.simactuators:

####################
lsst.ts.simactuators
####################

Python actuator simulators intended for the simulations in Commandable SAL Components (CSCs).

.. _lsst.ts.simactuators-using:

Using lsst.ts.simactuators
==========================

The primary classes are:

* `PointToPointActuator` simulates a simple actuator that moves to a give position at a constant speed, then stops.
  Examples include hexapod struts and telescope dome shutters.
* `CircularPointToPointActuator` a point to point actuator that moves in a circle with no limits.
  Examples include some telescope dome azimuth actuators.
* `TrackingActuator` simulates an actuator that slews to and tracks a path defined by a series of position, velocity, time triplets.
  Examples include telescope azimuth, altitude and rotator axes.
* `CircularTrackingActuator` a tracking actuator that moves in a circle with no limits.

Contributing
============

``lsst.ts.simactuators`` is developed at https://github.com/lsst-ts/ts_simactuators.
You can find Jira issues for this module using `labels=ts_simactuators <https://jira.lsstcorp.org/issues/?jql=project%20%3D%20DM%20AND%20labels%20%3D%20ts_simactuators>`_.


Python API reference
====================

.. automodapi:: lsst.ts.simactuators
   :no-main-docstr:
   :no-inheritance-diagram:

Version History
===============

.. toctree::
    version_history
    :maxdepth: 1
