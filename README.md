![Public domain serpentine logo](https://publicdomainvectors.org/photos/rygle_Snake_Colour_Outline.png)

# Serpentine binning


[![Build Status](https://travis-ci.org/koszullab/serpentine.svg?branch=master)](https://travis-ci.org/koszullab/serpentine)
[![Appveyor Status](https://ci.appveyor.com/api/projects/status/github/koszullab/serpentine?svg=true)](https://ci.appveyor.com/project/baudrly/serpentine)
[![codecov](https://codecov.io/gh/koszullab/serpentine/branch/master/graph/badge.svg)](https://codecov.io/gh/koszullab/serpentine)
[![Read the docs](https://readthedocs.org/projects/serpentine/badge)](https://serpentine.readthedocs.io)
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/koszullab/serpentine/master?filepath=doc%2Fnotebooks%2Fdemo_yeast.ipynb)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic%202.0-0298c3.svg)](https://opensource.org/licenses/Artistic-2.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

Locally smearing noisy regions in Hi-C contact maps as a prelude to differential analyses

#### Table of contents

1. [Synopsis](#synopsis)
2. [Installation](#installation)
3. [Documentation](#documentation)
4. [Authors](#authors)
5. [Copyright and license](#copyright-and-license)

## Synopsis

Use it as a Python 3 library:
```python
   import numpy as np
   import serpentine as sp

   A = np.loadtxt('./demos/A.csv')
   B = np.loadtxt('./demos/B.csv')
   trend, threshold = sp.MDbefore(A, B, show=False);

   sA, sB, sK = sp.serpentin_binning(A, B, threshold, threshold / 5)
```

Or as a standalone UNIX tool:
```
$ serpentine --help
   Serpentine binning

   An implementation of the so-called 'serpentine binning' procedure described
   in Scolari et al.

   Command line::

    Usage:
        serpentine.py [<matrixA>] [<matrixB>] [--threshold=auto] [--verbose]
                      [--min-threshold=auto] [--trend=high] [--triangular]
                      [--limit=3] [--demo] [--demo-size=500]

    Arguments:
        matrixA                         The first input matrix, in plain text
                                        CSV format. Optional in demo mode.
        matrixB                         The second input matrix, in plain text
                                        CSV format. Optional in demo mode or
                                        single binning mode.

    Options:
        -h, --help                      Display this help message.
        --version                       Display the program's current version.
        -t auto, --threshold auto       Threshold value to trigger binning.
                                        [default: auto]
        -m auto, --min-threshold auto   Minimum value to force trigger binning
                                        in either matrix. [default: auto]
        --trend high                    Trend to subtract to the differential
                                        matrix, possible values are "mean":
                                        equal amount of positive and negative
                                        differences, and "high": normalize
                                        at the regions with higher coverage.
                                        [default: high]
        --triangular                    Treat the matrix as triangular,
                                        useful when plotting matrices adjacent
                                        to the diagonal. [default: False]
        --limit 3                       Set the z-axis limit on the
                                        plot of the differential matrix.
                                        [default: 3]
        --demo                          Run a demo on randomly generated
                                        matrices. [default: False]
        --demo-size 500                 Size of the test matrix for the demo.
                                        [default: 500]
        -v, --verbose                   Show verbose output. [default: False]
```


## Installation

```sh
   sudo pip3 install -e git+https://github.com/koszullab/serpentine.git@master#egg=serpentine
```

## Documentation

Executing the command `serpentine  --help` will give you a brief help of the command line tool. For a detailed reference to the python library functions, please 
read the documentation on [readthedocs website](https://serpentine.readthedocs.io/en/latest/)

## Authors

Cluster Buster ([scovit](https://github.com/scovit), a.k.a. Vittore F. Scolari),
Lyamovich ([baudrly](https://github.com/baudrly), a.k.a. Lyam Baudry)

## Copyright and license

Copyright © 2017 Institut Pasteur, this software has been developed in
the Regulation Spatiale des Chromosomes team of Pasteur Institut,
Paris, France.

This library is free software; you can redistribute it and/or modify
it under the Artistic License.
