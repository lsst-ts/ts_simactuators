from setuptools import setup, find_namespace_packages

install_requires = []
tests_require = ["pytest", "pytest-cov", "pytest-flake8", "asynctest", "dds"]
dev_requires = install_requires + tests_require + ["documenteer[pipelines]"]
scm_version_template = """# Generated by setuptools_scm
__all__ = ["__version__"]

__version__ = "{version}"
"""

setup(
    name="ts_simactuators",
    description="Simulated actuators to simulation mode in Python CSCs.",
    use_scm_version={
        "write_to": "python/lsst/ts/simactuators/version.py",
        "write_to_template": scm_version_template,
    },
    setup_requires=["setuptools_scm", "pytest-runner"],
    install_requires=install_requires,
    package_dir={"": "python"},
    packages=find_namespace_packages(where="python"),
    package_data={"": ["*.rst", "*.yaml"]},
    data_files=[],
    scripts=[],
    tests_require=tests_require,
    extras_require={"dev": dev_requires},
    license="GPL",
    project_urls={
        "Bug Tracker": "https://jira.lsstcorp.org/secure/Dashboard.jspa",
        "Source Code": "https://github.com/lsst-ts/ts_simactuators",
    },
)
