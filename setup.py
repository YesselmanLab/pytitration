#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()


with open("README.md", "r", encoding="UTF-8") as f:
    readme = f.read()

with open("requirements.txt", "r", encoding="UTF-8") as f:
    requirements = f.read().splitlines()

setup(
    name="pytitration",
    version="0.1.0",
    description="a short script to compute kd and hill coefficients from chemical mapping data",
    long_description=readme,
    long_description_content_type="test/markdown",
    author="Joe Yesselman",
    author_email="jyesselm@unl.edu",
    url="https://github.com/jyesselm/pytitration",
    packages=[
        "pytitration",
    ],
    package_dir={"pytitration": "pytitration"},
    py_modules=["pytitration/cli"],
    include_package_data=True,
    #install_requires=requirements,
    zip_safe=False,
    keywords="pytitration",
    classifiers=[
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: Implementation :: PyPy",
    ],
    entry_points={"console_scripts": ["pytitration = pytitration.cli:cli"]},
)
