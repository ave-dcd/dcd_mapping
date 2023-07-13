from metadata_process import metadata_obtain
from blat_alignment import mave_to_blat
from transcript_selection import main
from vrs_mapping import vrs_mapping
from time import time

def main_map(urn):
    dat = metadata_obtain(urn)
    blat_dict = mave_to_blat(dat)
    tra = main(blat_dict,dat)
    vrs = vrs_mapping(dat,tra,blat_dict,"")
    return vrs
