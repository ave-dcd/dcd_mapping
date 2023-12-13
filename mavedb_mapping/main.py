from mavedb_mapping.metadata_process import metadata_obtain
from mavedb_mapping.blat_alignment import mave_to_blat
from mavedb_mapping.transcript_selection import main
from mavedb_mapping.vrs_mapping import vrs_mapping
from mavedb_mapping import data_file_path
from mavedb_mapping.output import output
from mavedb_mapping.find_start_pos import if_start_not_first
import json

"""Function to create a mapping given a scoreset and its scores"""


def main_map(scoreset, scores_csv):
    # Function to obtain expected inputs from scoreset and scores data

    dat, scores = metadata_obtain(scoreset, scores_csv)
    # dat and scores are expected inputs

    # Mapping process
    blat_dict = mave_to_blat(dat)
    transcripts = main(blat_dict, dat)
    # for seq where start pos is not first
    c = if_start_not_first(dat, scores)
    if c:
        transcripts["start"] = c
    vrs = vrs_mapping(dat, transcripts, blat_dict, scores)
    out = output(dat, transcripts, vrs, blat_dict)
    return out
