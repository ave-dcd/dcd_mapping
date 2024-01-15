# dcd-mapping

## Usage

Use the `dcd-map` command with a scoreset URN, eg

```shell
$ dcd-map urn:mavedb:00000083-c-1
```

Output is saved in the format `<URN>_mapping_results.json` in the directory specified by the environment variable `MAVEDB_STORAGE_DIR`, or `~/.local/share/dcd-mapping` by default.

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

## Data

Mapped MaveDB Scoresets can be downloaded here: https://mavedb-mapping.s3.us-east-2.amazonaws.com/mappings.tar.gz
