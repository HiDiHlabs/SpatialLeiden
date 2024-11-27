# SpatialLeiden

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Code style: Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)
[![Checked with mypy](https://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://github.com/pre-commit/pre-commit)

``SpatialLeiden`` is an implementation of
[Multiplex Leiden clustering](https://leidenalg.readthedocs.io/en/stable/multiplex.html)
that can be used to cluster spatially resolved omics data.

``SpatialLeiden`` integrates with the [scverse](https://scverse.org/) by leveraging
[scanpy](https://scanpy.readthedocs.io/) and [anndata](https://anndata.readthedocs.io/)
but can also be used independently.

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

Müller-Bötticher, N., Sahay, S., Eils, R., and Ishaque, N.
"SpatialLeiden - Spatially-aware Leiden clustering"
bioRxiv (2024) https://doi.org/10.1101/2024.08.23.609349

```
@article {spatialleiden2024,
	author = {Müller-Bötticher, Niklas and Sahay, Shashwat and Eils, Roland and Ishaque, Naveed},
	title = {SpatialLeiden - Spatially-aware Leiden clustering},
	year = {2024},
	doi = {10.1101/2024.08.23.609349},
	journal = {bioRxiv},
	publisher = {Cold Spring Harbor Laboratory}
}
```

## Versioning

This project follows the [SemVer](https://semver.org/) guidelines for versioning.

## License

This project is licensed under the MIT License - for details please refer to the
[LICENSE](./LICENSE) file.
