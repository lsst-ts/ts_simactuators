###############
ts_simactuators
###############

``ts_simactuators`` is a Python package in the `LSST Science Pipelines <https://pipelines.lsst.io>`_
that provides simulators for actuators.

This is primarily intended to support simulation in Commandable SAL Components (CSCs).

`Documentation <https://ts-simactuators.lsst.io>`_

The package is compatible with ``setuptools``, as well as the `eups <https://github.com/RobertLuptonTheGood/eups>`_ package management system and ``scons`` build system.
Assuming you have the basic Vera C. Rubin LSST DM stack installed you can do the following, from within the package directory:

* ``setup -r .`` to setup the package and dependencies, at which point the unit tests can be run and the package can be used "in place".
* ``pytest`` to run the unit tests.
* ``python setup.py install`` to install the software.
* ``package-docs build`` to build the documentation.
  This requires ``documenteer``; see `building single package docs <https://developer.lsst.io/stack/building-single-package-docs.html>`_ for installation instructions.

This code is automatically formatted by ``black`` using a git pre-commit hook.
To enable this:

* Install the ``black`` Python package.
* Run ``git config core.hooksPath .githooks`` once in this repository.
