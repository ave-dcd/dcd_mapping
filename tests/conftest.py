import pytest
from mavedb_mapping.blat_alignment import mave_to_blat
from mavedb_mapping.metadata_process import metadata_obtain
from mavedb_mapping.transcript_selection import main
from mavedb_mapping.vrs_mapping import vrs_mapping
from mavedb_mapping import data_file_path

import requests
import pandas as pd


def blat_for_tests(dat):
    """
    Function that mocks the BLAT Alignment outputs
    """
    seq = dat["target_sequence"]
    type_ = dat["target_sequence_type"]
    database = "hg38"
    blat_url = f"https://genome.ucsc.edu/cgi-bin/hgBlat?userSeq={seq}&type={type_}&db={database}&output=json"
    response = requests.get(blat_url)
    blat_output = response.json()
    print(blat_output)
    hsp = blat_output["blat"][0]
    if "fix" in hsp[13]:
        hsp = blat_output["blat"][1]
    diff_list = hsp[-3].split(",")
    hit_starts = hsp[-1].split(",")
    query_starts = hsp[-2].split(",")

    cov = [int(x) for x in diff_list]
    coverage = (sum(cov) / len(seq)) * 100

    query_list = []
    for i in range(1, len(query_starts)):
        temp = f"[{query_starts[i-1]}:{query_starts[i]}]"
        query_list.append(temp)
    temp = f"[{query_starts[-1]}:{int(query_starts[-1])+int(diff_list[-1])}]"
    query_list.append(temp)
    hit_list = []
    for i in range(len(query_starts)):
        temp = f"[{hit_starts[i]}:{int(hit_starts[i])+int(diff_list[i])}]"
        hit_list.append(temp)
    data = {"query_ranges": query_list, "hit_ranges": hit_list}
    result_df = pd.DataFrame(data)
    if hsp[8] == "+":
        strand = 1
    else:
        strand = -1
    chr = hsp[13].replace("chr", "")

    mave_blat_dict = {
        "chrom": chr,
        "strand": strand,
        "hits": result_df,
        "coverage": coverage,
    }
    return mave_blat_dict


@pytest.fixture
def mock_fun1(monkeypatch):
    monkeypatch.setattr("mavedb_mapping.blat_alignment.mave_to_blat", blat_for_tests)


@pytest.fixture(
    scope="package", params=["urn:mavedb:00000041-a-1", "urn:mavedb:00000001-a-4"]
)
def full_mapping(request):
    scoreset_path = f"{data_file_path}{request.param}"
    scores_path = f"{data_file_path}scores-{(request.param)[11:]}"
    scores_csv = open(scores_path)
    with open(scoreset_path) as scoreset:
        mave_dat,scores = metadata_obtain(scoreset,scores_csv)
    mave_blat = blat_for_tests(mave_dat)
    mappings_dict = main(mave_blat, mave_dat)
    vrs_mapped = vrs_mapping(mave_dat, mappings_dict, mave_blat, scores)
    return vrs_mapped, mave_dat["urn"]


@pytest.fixture()
def obtain_transcripts(request):
    scoreset_path = f"{data_file_path}{request.param}"
    scores_path = f"{data_file_path}scores-{(request.param)[11:]}"
    scores_csv = open(scores_path)
    with open(scoreset_path) as scoreset:
        mave_dat,scores = metadata_obtain(scoreset,scores_csv)
    mave_blat = blat_for_tests(mave_dat)
    mappings_dict = main(mave_blat, mave_dat)
    return mappings_dict
