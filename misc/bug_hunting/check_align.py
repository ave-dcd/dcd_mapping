import logging
import pickle

from cool_seq_tool.schemas import Strand

from dcd_mapping.align import align
from dcd_mapping.resources import get_scoreset_metadata

logging.basicConfig(
    filename="dcd-check-align.log",
    format="%(asctime)s %(levelname)s:%(name)s:%(message)s",
    level=logging.WARNING,
    force=True,
)

_logger = logging.getLogger(__name__)

with open("notebooks/analysis/results/mave_blat.pickle", "rb") as f:
    mave_blat_dict = pickle.load(f)

with open("misc/bug_hunting/human_urns.txt", "r") as f:
    urns = [line.strip() for line in f.readlines()]

strand_reformat = {1: Strand.POSITIVE, -1: Strand.NEGATIVE}


def chrom_reformat(chrom):
    return f"chr{chrom}"


for urn in urns:
    print(f"Checking {urn}...")
    try:
        metadata = get_scoreset_metadata(urn)
        alignment = align(metadata, False, True)
    except Exception as e:
        _logger.error("%s error: %s", urn, e)
        continue
    if urn not in mave_blat_dict:
        continue
    original = mave_blat_dict[urn]
    try:
        for name, actual, expected in [
            ("chromosome", alignment.chrom, chrom_reformat(original["chrom"])),
            ("strand", alignment.strand, strand_reformat[original["strand"]]),
        ]:
            if actual != expected:
                _logger.error(
                    "%s %s mismatch: %s (actual) vs %s (expected)",
                    urn,
                    name,
                    actual,
                    expected,
                )
    except Exception as e:
        _logger.error("%s exception: %s", urn, e)
