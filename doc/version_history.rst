.. py:currentmodule:: lsst.ts.simactuators

.. _lsst.ts.simactuators.version_history:

###############
Version History
###############

v2.2.11
------

* Update the version of ts-conda-build to 0.5 in the conda recipe.

v2.2.10
------

* Update the version of ts-conda-build to 0.4 in the conda recipe.

v2.2.9
------

* Fix a unit test warning.
* Modernize unit tests to use bare asserts.
* Use ts_pre_commit_config.
* Jenkinsfile: use the shared library.
* Remove scons support.

v2.2.8
------

* Modernize Jenkinsfile and do not run as root.
* pre-commit: update black to 23.1.0, isort to 5.12.0, mypy to 1.0.0, and pre-commit-hooks to v4.4.0.

Requires:

* ts_utils 1

v2.2.7
------

* Run isort.
* Add isort and mypy to pre-commit and update other pre-commit tasks.
* Modernize conda/meta.yaml.

Requires:

* ts_utils 1

v2.2.6
------

* Build with pyproject.toml.

Requires:

* ts_utils 1

v2.2.5
------

* Fix a black formatting issue.

Requires:

* ts_utils 1

v2.2.4
------

* Fix a new mypy error by not checking DM's `lsst/__init__.py` files.

Requires:

* ts_utils 1

v2.2.3
------

* Use ts_utils instead of ts_salobj.
* Add type annotations and run mypy.
* Add a Jenkinsfile.

Requires:

* ts_utils 1

v2.2.2
------

* Use `unittest.IsolatedAsyncioTestCase` instead of the abandoned asynctest package.
* Modernize doc/conf.py for documenteer 0.6.
* Format the code with black 20.8b1.

Requires:

* ts_salobj 6

v2.2.1
------

Changes:

* `RampGenerator`: use a safer way to copy an array.

Requires:

* ts_salobj 6

v2.2.0
------

Changes:

* Add `CosineGenerator` and `RampGenerator` generator functors.
  These are designed to be used in CSC commanders.

Requires:

* ts_salobj 6

v2.1.1
------

Changes:

* Update Jenkinsfile.conda to use the shared library.
* Pin the versions of ts_idl and ts_salobj in conda/meta.yaml.

Requires:

* ts_salobj 5.15 or 6

v2.1.0
------

Changes:

* Add ``tai`` argument to `BasePointToPointActuator.stop`.
* Pin the version of ``black`` in ``conda/meta.yaml``.

Requires:

* ts_salobj 5.15 or 6

v2.0.0
------

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

Requires:

* ts_salobj 5.15

v1.0.1
------

Changes:

* Add ``tests/test_black.py`` to verify that files are formatted with black.
  This requires ts_salobj 5.11 or later.
* Fix flake8 warnings about f strings with no {}.
* Update ``.travis.yml`` to remove ``sudo: false`` to github travis checks pass once again.

Requires:

* ts_salobj 5.11

v1.0.0
------

Changes:

* Format with black.
* Add a revision history.

Requires:

* ts_salobj 5

v0.2.1
------

Add setuptools and conda build compatibility.

Requires:

* ts_salobj 5

v0.2.0
------

Additional cleanups that I forgot to make for the first version, plus fixing a bug in slew.py.

Requires:

* ts_salobj 5

v0.1.0
------

First release.

Requires:

* ts_salobj 5
