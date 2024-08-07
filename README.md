# MEEGLET
> Morlet wavelets for M/EEG analysis, [ˈmiːglɪt]

This package provides a lean implementation of Morlet wavelets *designed for power-spectral analysis of M/EEG resting-state signals*.

- Distinct __frequency-domain parametrization of Morlet wavelets__
- Established __spectral M/EEG metrics__ share same wavelet convolutions
- Harmonized & tested __Python__ and __MATLAB__ implementation __numerically equivalent__
- Comprehensive __mathematical documentation__


```python
import matplotlib.pyplot as plt
from meeglet import define_frequencies, define_wavelets, plot_wavelet_family

foi, sigma_time, sigma_freq, bw_oct, qt = define_frequencies(
    foi_start=1, foi_end=32, bw_oct=1, delta_oct=1
)

wavelets = define_wavelets(
    foi=foi, sigma_time=sigma_time, sfreq=1000., density='oct'
)

plot_wavelet_family(wavelets, foi, fmax=64)
plt.gcf().set_size_inches(9, 3)
```

## Documentation
|                          |                                                  |
|:-------------------------|:-------------------------------------------------|
| __Background__           | overview on scope, rationale & design choices    |
| __Python tutorials__     | M/EEG data analysis examples                     |
| __Python API__           | Documentation of Python functions and unit tests |
| __MATLAB functionality__ | MATLAB documentation and data analysis example   |

Use the left sidebar for navigating conveniently!

## Installation

## from PyPi

In your environment of choice, use pip to install meeglet:

```bash
pip install meeglet
```

### from the sources

Please clone the software, consider installing the dependencies listed in the \`environment.yml.

Then do in your conda/mamba environment of choice:

``` bash
pip install -e .
```

## Citation

When using our package, please cite our two reference articles:

Python implementation and covariance computation.

``` bibtex
@article{bomatter2024,
	author = {Bomatter, Philipp and Paillard, Joseph and Garces, Pilar and Hipp, J{\"o}rg and Engemann, Denis-Alexander},
	title = {Machine learning of brain-specific biomarkers from EEG},
	year = {2024},
	journal = {eBioMedicine},
	url = {https://doi.org/10.1016/j.ebiom.2024.105259},
	date = {2024/08/05},
	publisher = {Elsevier},
	isbn = {2352-3964},
	month = {2024/08/06},
	volume = {106},
}
```

General methodology, MATLAB implementation and power-envelope correlations.

``` bibtex
@article{hipp2012large,
  title={Large-scale cortical correlation structure of spontaneous oscillatory activity},
  author={Hipp, Joerg F and Hawellek, David J and Corbetta, Maurizio and Siegel, Markus and Engel, Andreas K},
  journal={Nature neuroscience},
  volume={15},
  number={6},
  pages={884--890},
  year={2012},
  publisher={Nature Publishing Group US New York}
}
```

## Related software

M/EEG features based on Morlet wavelets using the more familiar time-domain parametrization can be readily computed is sevaral major software packages for M/EEG analysis:

- [FieldTrip](https://www.fieldtriptoolbox.org/)
- [BrainStorm](https://neuroimage.usc.edu/brainstorm/)
- [MNE](https://mne.tools/stable/index.html)
- [MNE-Connectivity](https://mne.tools/mne-connectivity/stable/index.html)
