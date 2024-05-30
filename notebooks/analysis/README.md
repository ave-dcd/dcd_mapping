# Notebooks for "Mapping MAVE data for use in human genomics applications"

Code for data analysis and figure generation for ["Mapping MAVE data for use in human genomics applications" (Arbesfeld et. al.)](https://www.biorxiv.org/content/10.1101/2023.06.20.545702v1):

* `mavedb_mapping_analysis.ipynb`: This notebook applies the mapping algorithm to a set of 209 examined score sets from MaveDB, successfully creating mappings for ~2.5 million variant pairs across 207 score sets.
* `mapping_analysis.ipynb`: This notebook computes reference sequence concordance across the generated VRS mapping pairs. The notebook also computes the number of unique pre-mapped and post-mapped variants.
* `mavedb_scoreset_breakdown.ipynb`: This notebook generates the summary statistics that are described in the manuscript.

## Environment

A compatible Python environment can be generated using the included `requirements.txt` file.

First, create and activate a virtual environment of your preference. For example, using `virtualenv`:

```shell
python3 -m virtualenv venv
source venv/bin/activate
```

Then install all requirements in `requirements.txt`:

```shell
python3 -m pip install -r requirements.txt
```

## Layout

TODO
