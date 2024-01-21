# dcd-mapping

## Results

Mapped MaveDB Scoresets can be downloaded here: https://mavedb-mapping.s3.us-east-2.amazonaws.com/mappings.tar.gz

## Tool usage

Use the `dcd-map` command with a scoreset URN, eg

```shell
$ dcd-map urn:mavedb:00000083-c-1
```

Output is saved in the format `<URN>_mapping_results.json` in the directory specified by the environment variable `MAVEDB_STORAGE_DIR`, or `~/.local/share/dcd-mapping` by default.

## Setup

Following installation instructions for [CoolSeqTool](https://coolseqtool.readthedocs.io/en/0.4.0-dev1/install.html) and [Gene Normalizer](https://gene-normalizer.readthedocs.io/en/latest/install.html) should take care of the external data dependencies.

Note that Gene Normalizer's `pg` dependency group must be installed to make use of the PostgreSQL-based backend:

```shell
python3 -m pip install 'gene-normalizer[pg]'
```

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

When debugging, use flags `-c` and `-d` to enable caching of BLAT alignment and debug logging, respectively.
