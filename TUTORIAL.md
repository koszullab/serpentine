# Serpentine binning tutorial

This tutorial aims at demonstrating use cases and for improving Hi-C contact maps with distribution-aware binning, helping readers reproduce the steps in Scolari et al., 2018 and documenting readers with the implementation. For a detailed step-by-step analysis, see the [serpentine notebook](notebooks/demo_yeast.ipynb).

## Dependencies

Python 3 with the following libraries:

* numpy
* matplotlib

## Testing

Run the following:

```bash
python serpentine.py --test --threshold 30 --min-threshold 3 --size 100
```

This randomly generates two 100x100 matrices of pixels between 0 and 10, binning both such that serpentines are not below 3 on average value in each of them and not below 30 in total value. It then plots both matrices as well as their log-ratio, binned and unbinned. Tweak the parameters to see how the matrices evolve.

## Binning prepared Hi-C datasets

Run the following:

```bash
wget https://github.com/koszullab/serpentine/dataset{1..2}.dat
python serpentines.py --inputs dataset1.dat dataset2.dat
```

This performs the same binning on prepared datasets from *Escherichia coli*. Currently supported formats for Hi-C datasets generated with [DADE](https://github.com/koszullab/DADE). They are essentially triangular, dense matrix files with the first row and the first column being used to indicate genomic position.

## Using the library

You can directly use the library's functions if you are already manipulating numpy contact maps in your analysis. Open a Python 3 console, and run the following:

```python
import numpy as np
from serpentine import serpentine_binning
from matplotlib import pyplot as plt

matrix3 = np.loadtxt("https://github.com/koszullab/serpentine/dataset3.dat", dtype=np.float64)
matrix4 = np.loadtxt("https://github.com/koszullab/serpentine/dataset4.dat", dtype=np.float64)
_, _, binned_matrix3, binned_matrix4, log_ratio, _ = serpentine_binning(matrix3, matrix4)

plt.imshow(log_ratio, cmap='seismic')
plt.show()

```

This is useful if your contact maps come from other sources and aren't in a DADE format.
