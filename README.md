# Random Wavelet Series (RWS)

These files reproduces all figures from the paper:

> Title: **Besov Regularity of Random Wavelet Series **

> Authors: **Andreas Horst , Thomas Jahn , Felix Voigtlaender**
> Link to preprint: ** **
## Requirements
To run the code you will need the following python packages
1. Numpy
2. Scipy
3. Pywavelets
4. Matplotlib
5. Jupyter
6. Seaborn

All of the packages can be installed via pip, e.g., pywavelets can be installed by calling pip install PyWavelets in a terminal.
For further installation guidelines please go to each package website and follow their instructions.

## Quickstart

1. Clone or download the repository
2. Navigate to the directory
3. Import packages as in the example code below

```python
# Imports
from RWS import Besov_Bernoulli_prior
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import seaborn as sns

```
Then you are ready to go!
## Example code 

See the file `Paper_figures.ipynb` for the code generating the figures in the paper. The code also serves as an example of how to use the code in RWS.py
