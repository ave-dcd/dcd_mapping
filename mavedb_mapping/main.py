from mavedb_mapping.metadata_process import metadata_obtain
from mavedb_mapping.blat_alignment import mave_to_blat
from mavedb_mapping.transcript_selection import main
from mavedb_mapping.vrs_mapping import vrs_mapping
from mavedb_mapping import data_file_path

def main_map(urn):
    dat = metadata_obtain(urn)
    blat_dict = mave_to_blat(dat)
    tra = main(blat_dict, dat)
    vrs = vrs_mapping(dat, tra, blat_dict)
    return vrs
