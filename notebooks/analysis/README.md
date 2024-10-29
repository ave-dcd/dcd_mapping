# Notebooks for "Mapping MAVE data for use in human genomics applications"

Code for data analysis and figure generation for ["Mapping MAVE data for use in human genomics applications" (Arbesfeld et. al.)](https://www.biorxiv.org/content/10.1101/2023.06.20.545702v1):

* [`mavedb_mapping.ipynb`](mavedb_mapping.ipynb): This notebook applies the mapping algorithm to a set of 209 examined score sets from MaveDB, successfully creating mappings for ~2.5 million variant pairs across 207 score sets.
* [`mapping_analysis.ipynb`](mapping_analysis.ipynb): This notebook computes reference sequence concordance across the generated VRS mapping pairs. The notebook also computes the number of unique pre-mapped and post-mapped variants.
* [`mavedb_scoreset_breakdown.ipynb`](mavedb_scoreset_breakdown.ipynb): This notebook generates the summary statistics that are described in the manuscript.

**Note that these notebooks are using the [0.1.3 release](https://github.com/ave-dcd/dcd_mapping/releases/tag/0.1.3) of the `dcd_mapping` library** -- they are intended to reflect the state of the code at the time of artifact generation, without any features that have been added since. The included `requirements.txt` file should produce an environment matching this expectation.

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

## Directory layout

After executing mapping code, this directory will contain working and output data in the following locations:

```
├── README.md
├── analysis_files
│   ├── mappings
│   │   └── <mapping output files>
│   └── <mapping checkpoint files>
├── experiment_scoresets.txt
├── mapping_analysis.ipynb
├── mave_mapping_fig_3b.R
├── mavedb_files
│   └── <Scoreset records and metadata from MaveDB>
├── mavedb_mapping.ipynb
├── mavedb_scoreset_breakdown.ipynb
└── requirements.txt
```
