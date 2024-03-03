# Getting started

> meeglet: Morlet wavelets for M/EEG analysis

Analysis of M/EEG signals in the frequency domain has proven very useful
and is widely used. This reflects a key contribution of periodic
processes to electrophysiological signals, which are well represented in
the frequency domain. Electrophysiological signal (as many other
naturally occurring signals) show a *log-frequency behavior*, i.e., the
magnitude and spectral width of oscillatory signals scale logarithmic
rather than linear with frequency (Buzsáki and Mizuseki 2014). The
commonly used Fourier transform is optimal for signals with a linear
frequency behavior. Frequency transformations that account for
*log-frequency behavior* may be better suited for electrophysiological
signals. Morlet wavelets, that are widely used for time-frequency
analyses (Tallon-Baudry et al. 1996), can be used to this end.

This package provides a lean implementation of Morlet wavelets using a
frequency-domain parametrization initially developed to facilitate the
spectral analysis of M/EEG resting-state signals (Hipp et al. 2012). The
package derives besides spectral power also the covariance matrix and
several other bi-variate and summary measures (coherence, imaginary
coherence, phase-locking value, phase-lag index, de-biased weighted
phase-lag index, power correlations, orthoganlized power correlations
and the global interaction measure). These representations can be used
for illustration and for further statistical analyses (e.g. for cluster
permutation statistics or as features for machine learning).

## Introduction

The implementation uses a log-frequency grid with a wavelet design that
increases spectral smoothing log-linearly with frequency. The
implementation parameterizes spectral spacing and smoothing in units of
*octaves* (`delta_oct` and `bw_oct`). We provide a python and a Matlab
implementation. In the following we describe the python version.

``` python
import numpy as np
import matplotlib.pyplot as plt

import meeglet
from meeglet import define_frequencies, define_wavelets

bw_oct = 0.5 # half an octave of standard deviation around center frequency
delta_oct = 1 # one octave spacing between frequencies of interest (foi)

foi, sigma_time, sigma_freq, bw_oct, qt = define_frequencies(
    foi_start=1, foi_end=32, delta_oct=delta_oct, bw_oct=bw_oct
)
```

The first step prepares frequencies of interest, `foi`, for which Morlet
wavelets are constructed.

``` python
wavelets  = define_wavelets(foi, sigma_time, sfreq=1000)

fig, axes = plt.subplots(len(wavelets), sharex=True)
axes = list(reversed(axes))
colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(wavelets)))

for ii, (w, *_) in enumerate(wavelets):
    w_range = np.arange(len(w)) - len(w)/2
    axes[ii].plot(w_range, w.real, color=colors[ii])
    axes[ii].plot(w_range, w.imag, color=colors[ii], linestyle='--', alpha=0.7)
    axes[ii].set_title(f'{foi[ii]} Hz', y=0.5, x=0.1)

plt.xlabel("Time (ms)")
```

The Morlet wavelets are widest in lower frequencies and shorter in
higher frequencies. Note also that the number of cycles is constant
across wavelets, which is in line with other implementations of Morlet
wavelets.

In addition, our Morlet wavelet family does not only have a log-linearly
increasing spectral width and a log-linearly decreasing temporal width,
but kernel width and frequency sampling are *coordinated*, i.e., fewer
wavelets are deployed at higher frequencies. The spectral distance
between two wavelets, i.e. spacing, is expressed in *octaves*, hence,
log-linear too.

``` python
delta_foi = np.c_[
    2 ** (np.log2(foi) - delta_oct / 2), 
    2 ** (np.log2(foi) + delta_oct / 2)
]
plt.figure()
plt.loglog(foi, sigma_freq, marker='o', base=2, label=r'$\sigma_f$', color='orange')
plt.loglog(foi, sigma_time, marker='o', base=2, label=r'$\sigma_t$', color='steelblue')
plt.loglog(delta_foi.T,  np.c_[sigma_freq, sigma_freq].T, color='orange')
plt.loglog(delta_foi.T, np.c_[sigma_time, sigma_time].T, color='steelblue')

plt.legend()
plt.xticks(ticks=foi, labels=foi)
plt.xlabel('Frequency (Hz)')
plt.ylabel(r'Bandwidth of Wavelet (temporal: $\sigma_t$, spectral: $\sigma_f$)')
plt.grid(True)
```

## Power spectral density in units of µ²/oct

The implementation provides the power spectrum in units of $µV^2/oct$.
While this is an unusual normalization for power spectral density, it is
a logical continuation embracing the logarithmic nature of
electrophysiological signals. It automatically leads to larger values at
higher frequencies compared to the traditional $µV^2/Hz$ normalization
accounting for lower amplitudes at higher frequencies.

## Key features of this Morlet wavelet implementation

- The log-linear parametrization in octaves enables intuitive reasoning
  about frequencies of interest.
- Motivated by resting-state EEG and proved useful in past
  pharmacological research & clinical applications.
- Helps avoid thinking in (*arbitrary*) frequency bands while expressing
  prior knowledge about log-linear scaling of brain structure and
  function (Buzsáki and Mizuseki 2014).
- A simple time-domain convolution is implemented, which is not
  optimized for speed or efficacy (Cohen 2019; Gramfort et al. 2013). In
  our experience, computation times are not an issue when using modern
  computers.
- The main reason for using a frequency domain representation based on
  Morlet wavelets over widely used Fourier transform based approaches is
  : Fourier transform implies a linear grid and constant spectral width
  and is therefor not optimal for log-frequency scaling of
  electrophysiological signals.
- Our implementation allows ignoring bad segments labeled as missing
  (NaN), propagating the number of effective samples to down-stream
  statistics like power or covariance.
- Being Morlet wavelets, the kernel families that can be obtained from
  our package can be used with other convolution implementations, if
  desired.

## Installation

Please clone the software, consider installing the dependencies listed
in the \`environment.yml.

Then do in your conda/mamba environment of choice:

``` bash

pip install -e .
```

## Citation

Please cite our reference articles when usinng this package.

For the software library:

``` bibtex
@article {bomatter2023,
    author = {Philipp Bomatter and Joseph Paillard and Pilar Garces and Joerg F Hipp and Denis A Engemann},
    title = {Machine learning of brain-specific biomarkers from EEG},
    elocation-id = {2023.12.15.571864},
    year = {2023},
    doi = {10.1101/2023.12.15.571864},
    publisher = {Cold Spring Harbor Laboratory},
    URL = {https://www.biorxiv.org/content/early/2023/12/21/2023.12.15.571864},
    eprint = {https://www.biorxiv.org/content/early/2023/12/21/2023.12.15.571864.full.pdf},
    journal = {bioRxiv}
}
```

And for the general methodology:

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

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-buzsaki2014log" class="csl-entry">

Buzsáki, György, and Kenji Mizuseki. 2014. “The Log-Dynamic Brain: How
Skewed Distributions Affect Network Operations.” *Nature Reviews
Neuroscience* 15 (4): 264–78.

</div>

<div id="ref-cohen2019better" class="csl-entry">

Cohen, Michael X. 2019. “A Better Way to Define and Describe Morlet
Wavelets for Time-Frequency Analysis.” *NeuroImage* 199: 81–86.

</div>

<div id="ref-gramfort2013meg" class="csl-entry">

Gramfort, Alexandre, Martin Luessi, Eric Larson, Denis A Engemann,
Daniel Strohmeier, Christian Brodbeck, Roman Goj, et al. 2013. “MEG and
EEG Data Analysis with MNE-Python.” *Frontiers in Neuroscience*, 267.

</div>

<div id="ref-hipp2012large" class="csl-entry">

Hipp, Joerg F, David J Hawellek, Maurizio Corbetta, Markus Siegel, and
Andreas K Engel. 2012. “Large-Scale Cortical Correlation Structure of
Spontaneous Oscillatory Activity.” *Nature Neuroscience* 15 (6): 884–90.

</div>

<div id="ref-tallon1996stimulus" class="csl-entry">

Tallon-Baudry, Catherine, Olivier Bertrand, Claude Delpuech, and Jacques
Pernier. 1996. “Stimulus Specificity of Phase-Locked and
Non-Phase-Locked 40 Hz Visual Responses in Human.” *Journal of
Neuroscience* 16 (13): 4240–49.

</div>

</div>
