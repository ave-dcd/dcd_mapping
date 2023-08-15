from mavedb_mapping.metadata_process import metadata_obtain
from mavedb_mapping.blat_alignment import mave_to_blat
from mavedb_mapping.transcript_selection import main
from mavedb_mapping.vrs_mapping import vrs_mapping
from mavedb_mapping import data_file_path
from mavedb_mapping.output import output
import json


def main_map(urn, scores_csv):
    dat, scores = metadata_obtain(urn, scores_csv)
    blat_dict = mave_to_blat(dat)
    transcripts = main(blat_dict, dat)
    vrs = vrs_mapping(dat, transcripts, blat_dict, scores)
    out = output(dat, transcripts, vrs, blat_dict)
    return out
