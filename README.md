# SpatialLeiden

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Checked with mypy](https://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://github.com/pre-commit/pre-commit)
[![Docs](https://app.readthedocs.org/projects/spatialleiden/badge/?version=latest)](https://spatialleiden.readthedocs.io)
[![PyPI](https://img.shields.io/pypi/v/spatialleiden)](https://pypi.org/project/spatialleiden)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/spatialleiden/README.html)


``SpatialLeiden`` is an implementation of
[Multiplex Leiden clustering](https://leidenalg.readthedocs.io/en/stable/multiplex.html)
that can be used to cluster spatially resolved omics data.

``SpatialLeiden`` integrates with the [scverse](https://scverse.org/) by leveraging
[anndata](https://anndata.readthedocs.io/) but can also be used independently.

## Installation

`spatialleiden` is available on [PyPI](https://pypi.org/project/spatialleiden/) and
[bioconda](https://bioconda.github.io/recipes/spatialleiden/README.html).

```sh
# PyPI
pip install spatialleiden
```

```sh
# or conda
conda install bioconda::spatialleiden
```

For detailed installation instructions please refer to the
[documentation](https://spatialleiden.readthedocs.io/page/installation.html).

## Documentation

For documentation of the package please refer to the
[ReadTheDocs page](https://spatialleiden.readthedocs.io/).

## Citations

If you are using `spatialleiden` for your research please cite

Müller-Bötticher, N., Sahay, S., Eils, R., & Ishaque, N. (2025).
SpatialLeiden: spatially aware Leiden clustering.
*Genome Biology*, 26(1), 24. https://doi.org/10.1186/s13059-025-03489-7

```
@article{spatialleiden2025,
	author = {Müller-Bötticher, Niklas and Sahay, Shashwat and Eils, Roland and Ishaque, Naveed},
	title = {SpatialLeiden: spatially aware Leiden clustering},
	journal = {Genome Biology},
	year = {2025},
	month = {Feb},
	day = {07},
	volume = {26},
	number = {1},
	pages = {24},
	doi = {10.1186/s13059-025-03489-7},
	url = {https://doi.org/10.1186/s13059-025-03489-7}
}
```

## Versioning

This project follows the [SemVer](https://semver.org/) guidelines for versioning.

## License

This project is licensed under the MIT License - for details please refer to the
[LICENSE](./LICENSE) file.
