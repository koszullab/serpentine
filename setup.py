#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Serpentine binning package for Hi-C contact maps
"""

from setuptools import setup, find_packages

CLASSIFIERS = [
    "Development Status :: 3 - Alpha",
    "Intended Audience :: Science/Research",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: Artistic License",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3 :: Only",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Scientific/Engineering :: Visualization",
    "Operating System :: Microsoft :: Windows",
    "Operating System :: POSIX",
    "Operating System :: Unix",
    "Operating System :: MacOS",
]

name = "serpentine"

MAJOR = 0
MINOR = 1
MAINTENANCE = 1
VERSION = "{}.{}.{}".format(MAJOR, MINOR, MAINTENANCE)

LICENSE = "Artistic License 2.0"
URL = "https://github.com/koszullab/serpentine"

DESCRIPTION = __doc__.strip("\n")

with open("README.md") as f:
    LONG_DESCRIPTION = f.read()

with open("requirements.txt", "r") as f:
    REQUIREMENTS = f.read().splitlines()

with open("serpentine/version.py", "w") as f:
    f.write("__version__ = '{}'\n".format(VERSION))

setup(
    name=name,
    author="vittore.scolari@pasteur.fr",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    version=VERSION,
    license=LICENSE,
    classifiers=CLASSIFIERS,
    url=URL,
    packages=find_packages(exclude=["demos"]),
    install_requires=REQUIREMENTS,
    include_package_data=True,
    python_requires=">=3.4",
    long_description_content_type="text/markdown",
    entry_points={
        "console_scripts": ["serpentine=serpentine.serpentine:_main"]
    },
)
