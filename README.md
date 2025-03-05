# dcd-map: Map MaveDB data to computable and interoperable variant objects

[![image](https://img.shields.io/pypi/v/dcd_mapping.svg)](https://pypi.python.org/pypi/dcd_mapping)
[![image](https://img.shields.io/pypi/l/dcd_mapping.svg)](https://pypi.python.org/pypi/dcd_mapping)
[![image](https://img.shields.io/pypi/pyversions/dcd_mapping.svg)](https://pypi.python.org/pypi/dcd_mapping)
[![Actions status](https://github.com/ave-dcd/dcd_mapping/actions/workflows/checks.yaml/badge.svg)](https://github.com/ave-dcd/dcd_mapping/actions/checks.yaml)
[![DOI](https://zenodo.org/badge/472473437.svg)](https://zenodo.org/doi/10.5281/zenodo.11406657)

<!-- description -->

This library implements a novel method for mapping [MaveDB scoreset data](https://mavedb.org/) to [GA4GH Variation Representation Specification (VRS 2.0)](https://vrs.ga4gh.org/en/2.x/) objects, enhancing interoperability for genomic medicine applications. See [Arbesfeld et. al. (2024)](https://www.biorxiv.org/content/10.1101/2023.06.20.545702) for a preprint edition of the mapping manuscript, or [download the resulting mappings directly](https://mavedb-mapping.s3.us-east-2.amazonaws.com/mappings_20250220.tar.gz).

If you make use of this software or the resultant mappings, please cite the manuscript:

> Jeremy A. Arbesfeld, Estelle Y. Da, James S. Stevenson, Kori Kuzma, Anika Paul, Tierra Farris, Benjamin J. Capodanno, Sally B. Grindstaff, Kevin Riehle, Nuno Saraiva-Agostinho, Jordan F. Safer, Aleksandar Milosavljevic, Julia Foreman, Helen V. Firth, Sarah E. Hunt, Sumaiya Iqbal, Melissa S. Cline, Alan F. Rubin, Alex H. Wagner. bioRxiv 2023.06.20.545702; doi: [https://doi.org/10.1101/2023.06.20.545702](https://doi.org/10.1101/2023.06.20.545702)

<!-- /description -->

## Prerequisites

* Universal Transcript Archive (UTA): see [README](https://github.com/biocommons/uta?tab=readme-ov-file#installing-uta-locally) for setup instructions. Users with access to Docker on their local devices can use the available Docker image; otherwise, start a relatively recent (version 14+) PostgreSQL instance and add data from the available database dump.
* SeqRepo: see [README](https://github.com/biocommons/biocommons.seqrepo?tab=readme-ov-file#requirements) for setup instructions. The SeqRepo data directory must be writeable; see specific instructions [here](https://github.com/biocommons/biocommons.seqrepo/blob/main/docs/store.rst) for more.
* Gene Normalizer: see [documentation](https://gene-normalizer.readthedocs.io/0.3.0-dev1/install.html) for data setup instructions.
* blat: Must be available on the local PATH and executable by the user. Otherwise, its location can be set manually with the `BLAT_BIN_PATH` env var. See the [UCSC Genome Browser FAQ](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3) for download instructions.


## Installation

Install from [PyPI](https://pypi.python.org/pypi/dcd_mapping):

```
python3 -m pip install dcd-mapping
```

## Usage

Use the `dcd-map` command with a scoreset URN, eg

```shell
$ dcd-map urn:mavedb:00000083-c-1
```

Output is saved in the format `<URN>_mapping_results_<ISO datetime>.json` in the directory specified by the environment variable `MAVEDB_STORAGE_DIR`, or `~/.local/share/dcd-mapping` by default.

Use `dcd-map --help` to see other available options.

## Notebooks

Notebooks for manuscript data analysis and figure generation are provided within `notebooks/analysis`. See [`notebooks/analysis/README.md`](notebooks/analysis/README.md) for more information.

## Development

Clone the repo

```
git clone https://github.com/ave-dcd/dcd_mapping
cd dcd_mapping
```

Create and activate a virtual environment

```
python3 -m virtualenv venv
source venv/bin/activate
```

Install as editable and with developer dependencies

```
python3 -m pip install -e '.[dev,tests]'
```

Add pre-commit hooks

```
pre-commit install
```

Run tests with `pytest`

```
pytest
```
