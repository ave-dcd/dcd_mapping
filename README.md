# dcd-map: Map MaveDB data to computable and interoperable variant objects

This library implements a novel method for mapping [MaveDB scoreset data](https://mavedb.org/) to [GA4GH Variation Representation Specification (VRS)](https://vrs.ga4gh.org/en/stable/) objects, enhancing interoperability for genomic medicine applications. See [Arbesfeld et. al. (2023)](https://www.biorxiv.org/content/10.1101/2023.06.20.545702v1) for more information, or [download the resulting mappings directly](https://mavedb-mapping.s3.us-east-2.amazonaws.com/mappings.tar.gz).

## Installation

Install from PyPI:

```
python3 -m pip install dcd-mapping
```

Also ensure the following data dependencies are available:

* Universal Transcript Archive (UTA): see [README](https://github.com/biocommons/uta?tab=readme-ov-file#installing-uta-locally) for setup instructions. Users with access to Docker on their local devices can use the available Docker image; otherwise, start a relatively recent (version 14+) PostgreSQL instance and add data from the available database dump.
* SeqRepo: see [README](https://github.com/biocommons/biocommons.seqrepo?tab=readme-ov-file#requirements) for setup instructions.
* Gene Normalizer: see [documentation](https://gene-normalizer.readthedocs.io/0.3.0-dev1/install.html) for data setup instructions.
* blat: Must be available on the local PATH and executable by the user. Otherwise, its location can be set manually with the `BLAT_BIN_PATH` env var. See the [UCSC Genome Browser FAQ](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3) for download instructions. For our experiments, we placed the binary in the same directory as these notebooks.

## Usage

Use the `dcd-map` command with a scoreset URN, eg

```shell
$ dcd-map urn:mavedb:00000083-c-1
```

Output is saved in the format `<URN>_mapping_results.json` in the directory specified by the environment variable `MAVEDB_STORAGE_DIR`, or `~/.local/share/dcd-mapping` by default.

## Notebooks

Notebooks for manuscript data analysis and figure generation are provided within `notebooks/analysis`. See `notebooks/analys/README.md` for more information.

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
