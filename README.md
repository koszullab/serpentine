![Public domain serpentine logo](https://publicdomainvectors.org/photos/rygle_Snake_Colour_Outline.png)

# Serpentine binning


[![Build Status](https://travis-ci.org/koszullab/serpentine.svg?branch=master)](https://travis-ci.org/koszullab/serpentine)
[![codecov](https://codecov.io/gh/koszullab/serpentine/branch/master/graph/badge.svg)](https://codecov.io/gh/koszullab/serpentine)
[![Read the docs](https://readthedocs.org/projects/serpentine/badge)](https://serpentine.readthedocs.io)
[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/koszullab/serpentine/master?filepath=notebooks%2Fdemo_yeast.ipynb)
[![License: Artistic-1.0](https://img.shields.io/badge/License-Artistic%201.0-0298c3.svg)](https://opensource.org/licenses/Artistic-1.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

Locally smearing noisy regions in Hi-C contact maps as a prelude to differential analyses

## Synopsis

```python
   import numpy as np
   import serpentine as sp

   A = np.loadtxt('./demos/A.csv')
   B = np.loadtxt('./demos/B.csv')
   trend, threshold = sp.MDbefore(A, B, show=False);

   sA, sB, sK = sp.serpentin_binning(A, B, threshold, threshold / 5)
```

## Installation

```sh
   sudo pip3 install -e git+https://github.com/koszullab/serpentine.git@master#egg=serpentine
```

## Documentation

Please read the documentation on [readthedocs website](https://serpentine.readthedocs.io/en/latest/)

## Authors

Cluster Buster (scovit, a.k.a. Vittore F. Scolari),
Lyamovich (baudrly, a.k.a. Lyam Baudry)

## Copyright and license

Copyright © 2017 Institut Pasteur, this software has been developed in
the Regulation Spatiale des Chromosomes team of Pasteur Institut,
Paris, France.

This library is free software; you can redistribute it and/or modify
it under the Artistic License.
