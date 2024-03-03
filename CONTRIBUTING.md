# Contributing to meeglet

## Types of contributions

If you want to contribute to meeglet, you can:

- report bugs
- contribute documentation
- request, propose, and contirbute features

## Bug reports

If you face an issue or find a bug, please open an issue and provide a minimal example to reproduce the bug.

If you are familiar with using git and Github, you can also directly propose a fix and include a unit test that would have failed before your fix.

We can then either address the issue directly or provide guidance on how to move forward.

## Contributing documentation

To contribute documentation you must be able to build documentation.

A great stating point is to taking an existing tutorial notebook, copy, rename and adjust it to display the content you would like to show.

But before embarking on this adventure, ask yourself whether you really want to provide new documentation specifically to meeglet or if you have a more general topic that is better documented in some of the downstream libraries or related projects, e.g. [MNE-Python](https://mne.tools) or [MNE-connectivity](https://mne.tools/mne-connectivity/stable/index.html).

For documenting API code, please open a pull request with the proposed changes.

In general, for all proposed documentation changes please provide evidence of successful rendering of the changes.

## Contributing source code

This section should contain information on

### How to access the project source code

Clone the repository, ideally create a fork on your GitHub account.

### General project layout

All source code lives in the ```./nbs``` folder.
All library code lives in the ```./nbs/api``` subfolder.
Everyhthing else is documentation.
The python library is generated from the exported code in the api folder.

### Any requirements to the development environment

Pleasse check out our requirements file.

In general, use Python >= 3.9 and a recent MNE version.

You must have a [Quarto](https://quarto.org/) and [nbdev](https://nbdev.fast.ai/) installed and you should have familiarized yourself with the usage.

### Code formatting guidelines

The code should be approximately pep8 formatted with the exception of doc strings, which follow the [nbdev style](https://nbdev.fast.ai/tutorials/best_practices.html).


### How to run the test suite

To test the the core functionality use:

```nbdev_test --do_print --path nbs/api```

To test all notebooks, including documentation, use:

```nbdev_test --do_print```

## Feature requests

If you wish to make feature requests, please open an issue and describe the desired functionality.
Ideally, provide some code that shows how the API could be used would your envisioned feature be implemented.

But before embarking on this adventure, ask yourself whether you really want to provide a new feature specifically to meeglet or if your feature request is better suited in some of the downstream libraries or related projects, e.g. [MNE-Python](https://mne.tools) or [MNE-connectivity](https://mne.tools/mne-connectivity/stable/index.html).
