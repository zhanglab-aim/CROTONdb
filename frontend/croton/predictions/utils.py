from itertools import islice

import logging
import pysam
import pandas as pd
from django.conf import settings

from .table_summary import TableSummary
from utils.hostname import timeit
logger = logging.getLogger(__name__)


# Default difference for querying tabix between start and end

DEFAULT_POSITION_RANGE = 100

# TBD for the rest of chromosomes
CHROMOSOME_START_POSITION = {
    '1': 69045,
    '10': 47027,
    '11': 168941,
    '12': 66853,
    '13': 19173852,
    '14': 18601105,
    '15': 19964658,
    'Y': 2786973,
}

# TBD for all the columns
# Mapping between columns ids and correct column names
HEADER_MAP = {
    0: "chromosome",
    1: "pos",
    2: "id",
    3: "ref",
    4: "alt",
    5: "pamid",
    6: "genename",
    7: "num",
    8: "start",
    9: "end",
    10: "strand",
    11: "inx",
    12: "ref_seq",
    13: "alt_seq",
    14: "ref_preds",
    15: "alt_preds",
    16: "abs_diff",
    17: "ref_del_frq",
    18: "alt_del_frq",
    19: "diff_del_frq",
    20: "ref_1ins",
    21: "alt_1ins",
    22: "diff_1ins",
    23: "ref_1del",
    24: "alt_1del",
    25: "diff_1del",
    26: "ref_onemod3",
    27: "alt_onemod3",
    28: "diff_onemod3",
    29: "ref_twomod3",
    30: "alt_twomod3",
    31: "diff_twomod3",
    32: "ref_frameshift",
    33: "alt_frameshift",
    34: "diff_frameshift"
}


@timeit
def fetch_tabix_croton_predictions(
    tabix_file_or_url,
    chrom=None,
    start=None,
    end=None,
):
    """ Fetches the prediction from the tabix file for given chromosome and positions.
    Args:
        tabix_file_or_url ([string]): path to the tabix file
        chrom ([string], optional): Chromosome name. Defaults to None.
        start ([int], optional): Start positionto be read. Defaults to None.
        end ([int], optional): End postion to be read. Defaults to None.

    Returns:
        [pandas.DataFrame]: Data frame that contains all the data read given the input parameters.
    """
    try:
        with pysam.TabixFile(tabix_file_or_url, encoding='utf-8') as t:
            if chrom:
                tabix_result = list(
                    t.fetch(
                        'chr{}'.format(chrom),
                        start,
                        end,
                        parser=pysam.asTuple(),
                    )
                )
            else:
                tabix_result = list(t.fetch(parser=pysam.asTuple()))

    except ValueError:
        logger.debug('No data found for positions {}-{}'.format(start, end))
        return ([], len(t.header) > 0)

    if not len(tabix_result):
        logger.debug('No data found for positions {}-{}'.format(start, end))
        # the header might still be available even no results
        return ([], len(t.header) > 0)

    if len(t.header):
        headers = t.header[0].split('\t')
        return (pd.DataFrame(tabix_result, columns=headers), True)

    return (pd.DataFrame(tabix_result), False)


@timeit
def convert_tabix_results(tabix_results, header_map={}, start=0, end=None, filter_keys=None):
    """Converts tabix results in data frame into the list of dictionary objects"""
    results = []

    # print(header_map)
    # header columns as they are
    orig_cols = tabix_results.columns
    # for idx, row in tabix_results.iterrows():
    # end=100
    for idx, row in islice(tabix_results.iterrows(), start, end):
        prediction = {}
        for col_idx, orig_col_name in enumerate(tabix_results.columns):
            col_name = orig_col_name
            if header_map and col_idx in header_map:

                # replace column name with the one in header_map
                col_name = header_map[col_idx]
                if filter_keys and not col_name in filter_keys:
                    # skipping this col_name because it is not supposed to be added
                    continue

            prediction[col_name] = row[col_idx]
        results.append(prediction)
    return results


def build_response_json(tabix_results_df, header_map):
    """Builds JSON object with all the results

    Args:
        tabix_results_df ([pandas.DataFrame]): predictions read from the tabix file
        header_map ([dict]): header_map required to change the indexes from data frame into column names

    Returns:
        response [json]: JSON object with all the data and calculated statistics
    """
    results_len = len(tabix_results_df)
    if results_len == 0:
        response = {}
        response["results"] = {}
        response["results"]["total"] = 0
        return response

    # REF_INDEX = 3
    # ALT_INDEX = 4
    # REF1MOD3_INDEX = 26
    # REF2MOD3_INDEX = 29
    # ALT1INS_INDEX = 21
    # ALT1MOD3_INDEX = 27
    # ALT2MOD3_INDEX = 30
    # assert header_map[REF_INDEX] == 'ref'
    # assert header_map[ALT_INDEX] == 'alt'
    # assert header_map[REF1MOD3_INDEX] == 'ref_onemod3'
    # assert header_map[REF2MOD3_INDEX] == 'ref_twomod3'
    # assert header_map[ALT1INS_INDEX] == 'alt_1ins'
    # assert header_map[ALT1MOD3_INDEX] == 'alt_onemod3'
    # assert header_map[ALT2MOD3_INDEX] == 'alt_twomod3'

    REF_INDEX = 'ref'
    ALT_INDEX = 'alt'
    REF1MOD3_INDEX = 'ref_onemod3'
    REF2MOD3_INDEX = 'ref_twomod3'
    ALT1INS_INDEX = 'alt_1ins'
    ALT1MOD3_INDEX = 'alt_onemod3'
    ALT2MOD3_INDEX = 'alt_twomod3'

    # FREQUENCIES
    ref_onemod3_mean = tabix_results_df[REF1MOD3_INDEX].astype(
        float).mean()*100
    ref_twomod3_mean = tabix_results_df[REF2MOD3_INDEX].astype(
        float).mean()*100
    no_frameshift = 100 - (ref_onemod3_mean+ref_twomod3_mean)

    # VALUE COUNTS
    snp_value_counts = tabix_results_df.value_counts(
        subset=[REF_INDEX, ALT_INDEX], sort=False)
    # to get index (columns): snp_value_counts.index[0], e.g.g ('A', 'C')
    # to get value (row): snp_value_counts.get(0): 135
    matrix = [[0, 0, 0, 0] for i in range(4)]
    n_lst = ['A', 'C', 'G', 'T']
    # for stacked bar plot it's better to have "to" mutation, so
    # matrix is a list of toA, toC, toG and toT lists
    for i, tpl in enumerate(snp_value_counts.index):
        matrix[n_lst.index(tpl[1])][n_lst.index(tpl[0])
                                    ] = snp_value_counts.get(i)

    filter_keys = ["chromosome", "pos", "id", "ref", "alt", "pamid", "start", "end", "strand", "ref_seq", "alt_seq",
                   "abs_diff", "ref_1ins", "alt_1ins", "diff_1ins", "ref_onemod3", "alt_onemod3", "diff_onemod3", "ref_twomod3", "alt_twomod3", "diff_twomod3", "ref_frameshift", "alt_frameshift", "diff_frameshift"]

    predictions_list = convert_tabix_results(
        tabix_results_df, header_map=header_map, filter_keys=filter_keys)

    # TableSummary
    ts = TableSummary(settings.TEST_ECDF_FILE, header_map)

    # print(ts.get_table_summary(predictions_list))
    summary_dict = ts.get_summary(predictions_list)

    response = {
        'summary':  summary_dict,
        'stats': {
            'freq': {'ref_onemod3_mean': ref_onemod3_mean,
                     'ref_twomod3_mean': ref_twomod3_mean,
                     'no_frameshift': no_frameshift},
            'snp': {'toA': matrix[0],
                    'toC': matrix[1],
                    'toG': matrix[2],
                    'toT': matrix[3]},
            'histograms': {
                'alt1ins': list((tabix_results_df[ALT1INS_INDEX].astype(float)*100).round(2)),
                'alt1mod3': list((tabix_results_df[ALT1MOD3_INDEX].astype(float)*100).round(2)),
                'alt2mod3': list((tabix_results_df[ALT2MOD3_INDEX].astype(float)*100).round(2)),
            },
        },
        'results': {
            'total': results_len,
            'predictions': [] if results_len == 0 else predictions_list
        }
    }
    return response
