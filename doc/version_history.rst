.. py:currentmodule:: lsst.ts.simactuators

.. _lsst.ts.simactuators.version_history:

###############
Version History
###############

v2.0.0
======

Breaking changes:

* Overhauled `PointToPointActuator` to use TAI dates (absolute times). Changes:
    * Replaced the ``current_position`` property with a ``position`` method that takes an optional TAI time.
    * Replaced the ``moving`` property with a ``moving`` method that takes an optional TAI time
    * Replaced the ``remaining_time`` property with a ``remaining`` method that takes an optional TAI time
    * Made the ``start_position`` constructor argument optional.
      The default value matches `TrackingActuator`.
    * Modified ``set_position`` to return the move duration.
    * Added an optional ``start_tai`` argument to ``set_position``.
      The default is the current time (the same behavior as before).
    * Added ``start_tai`` and ``end_tai`` properties.
    * Added a ``velocity`` method.

Other changes:

* Added `CircularPointToPointActuator`
* Added the optional ``initial_position`` constructor argument to `TrackingActuator`.
  The default gives the same behavior as before.
* Added `BasePointToPointActuator`
* Added `Direction` enumeration.

v1.0.1
======

Changes:

* Add ``tests/test_black.py`` to verify that files are formatted with black.
  This requires ts_salobj 5.11 or later.
* Fix flake8 warnings about f strings with no {}.
* Update ``.travis.yml`` to remove ``sudo: false`` to github travis checks pass once again.

Requires:

* ts_salobj 5.11

v1.0.0
======

Changes:

* Format with black.
* Add a revision history.

Requires:

* ts_salobj 5

v0.2.1
======

Add setuptools and conda build compatibility.

Requires:

* ts_salobj 5

v0.2.0
======

Additional cleanups that I forgot to make for the first version, plus fixing a bug in slew.py.

Requires:

* ts_salobj 5

v0.1.0
======

First release.

Requires:

* ts_salobj 5
