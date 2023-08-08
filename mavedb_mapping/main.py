from mavedb_mapping.metadata_process import metadata_obtain
from mavedb_mapping.blat_alignment import mave_to_blat
from mavedb_mapping.transcript_selection import main
from mavedb_mapping.vrs_mapping import vrs_mapping
from mavedb_mapping import data_file_path


def main_map(urn, scores_csv):
    dat = metadata_obtain(urn)
    print(dat)
    blat_dict = mave_to_blat(dat)
    print(blat_dict)
    transcripts = main(blat_dict, dat)
    print(transcripts)
    vrs = vrs_mapping(dat, transcripts, blat_dict, scores_csv)
    print(vrs)
    return vrs
