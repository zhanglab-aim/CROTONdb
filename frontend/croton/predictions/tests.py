from django.test import TestCase

import pprint
import statistics
from .utils import fetch_tabix_croton_predictions
from .utils import fetch_tabix_croton_predictions, convert_tabix_results, build_response_json
from .utils import HEADER_MAP, header_map_new
from .table_summary import TableSummary

test_resources_dir = 'predictions/test/'
test_ecdf_file = f'{test_resources_dir}test_gw_data.pkl'
test_no_header_tabix_file = f'{test_resources_dir}chrom1-test-no-header.tsv.gz'
# ds_model_name_default = 'test'
#
# Test tabix file taken from chromosome 1 file
# Header line starts with #
# Created based on VCF file using
# bgzip 10.tsv
# then indexing with tabix
# tabix -p vcf 10.tsv.gz
# chrom1 position from 69045 to 69072
"""
# chrom  pos     id      ref     alt     pamid   genename        num     start   end     strand  inx     ref_seq alt_seq ref_preds       alt_preds       abs_diff        ref_del_frq     alt_del_frq     diff_del_frq    ref_1ins        alt_1ins        diff_1ins       ref_1del        alt_1del        diff_1del       ref_onemod3     alt_onemod3     diff_onemod3    ref_twomod3     alt_twomod3     diff_twomod3    ref_frameshift  alt_frameshift  diff_frameshift
chr1    69045   rs1360507296    A       G       OR4F5|2 OR4F5   2       69024   69084   +       21      CTCTTCTTCAAGGTAACTGCAGAGGCTATTTCCTGGAATGAATCAACGAGTGAAACGAAT    CTCTTCTTCAAGGTAACTGCGGAGGCTATTTCCTGGAATGAATCAACGAGTGAAACGAAT    [0.9040266  0.03224605 0.18726781 0.27797917 0.40537828 0.6896229 ]     [0.90453994 0.03175539 0.18687397 0.29372668 0.407959   0.7080409 ]     0.018418014049530033    0.9040266275405884      0.904539942741394       -0.0005133152008056642  0.03224605321884155     0.03175538778305054     0.0004906654357910156   0.18726781010627747     0.18687397241592407     0.0003938376903533936   0.27797916531562805     0.2937266826629639      -0.015747517347335815   0.4053782820701599      0.4079590141773224      -0.002580732107162476   0.6896228790283203      0.7080408930778503      -0.018418014049530033
"""
test_tabix_file = f'{test_resources_dir}chrom1-test.tsv.gz'
test_no_header_tabix_file = f'{test_resources_dir}chrom1-test-no-header.tsv.gz'


#
# python manage.py test
# python manage.py test predictions
# python manage.py test predictions.tests.PredictionResultTestCase.test_histograms
#
#

class PredictionTestCase(TestCase):
    def test_fetch_tabix_all_results(self):
        pos = 69045
        start = pos - 1
        end = 69072
        result, is_header = fetch_tabix_croton_predictions(
            test_tabix_file,
            chrom=1,
            start=start,
            end=end,
        )
        self.assertTrue(is_header)
        self.assertTrue(isinstance(list(result.index), list))
        self.assertEqual(9, len(list(result.index)))

    def test_fetch_tabix_all_results_no_header(self):
        pos = 69045
        start = pos - 1
        end = 69072
        result, is_header = fetch_tabix_croton_predictions(
            test_no_header_tabix_file,
            chrom=1,
            start=start,
            end=end,
        )
        self.assertFalse(is_header)
        self.assertTrue(isinstance(list(result.index), list))
        self.assertEqual(9, len(list(result.index)))

    def test_fetch_tabix_few_results(self):
        pos = 69045
        start = pos - 1
        end = pos
        result, is_header = fetch_tabix_croton_predictions(
            test_tabix_file,
            chrom=1,
            start=start,
            end=end,
        )
        self.assertTrue(is_header)
        self.assertTrue(isinstance(list(result.index), list))
        self.assertEqual(3, len(list(result.index)))

    def test_fetch_tabix_few_results_no_header(self):
        pos = 69045
        start = pos - 1
        end = pos
        result, is_header = fetch_tabix_croton_predictions(
            test_no_header_tabix_file,
            chrom=1,
            start=start,
            end=end,
        )
        self.assertFalse(is_header)
        self.assertTrue(isinstance(list(result.index), list))
        self.assertEqual(3, len(list(result.index)))

    def test_fetch_tabix_no_results(self):
        pos = 1
        start = pos - 1
        end = pos
        result, is_header = fetch_tabix_croton_predictions(
            test_tabix_file,
            chrom=1,
            start=start,
            end=end,
        )
        self.assertTrue(is_header)
        self.assertTrue(isinstance(result, list))
        self.assertEqual(0, len(list(result)))

    def test_fetch_tabix_no_results_no_header(self):
        pos = 1
        start = pos - 1
        end = pos
        result, is_header = fetch_tabix_croton_predictions(
            test_no_header_tabix_file,
            chrom=1,
            start=start,
            end=end,
        )
        self.assertFalse(is_header)
        self.assertTrue(isinstance(result, list))
        self.assertEqual(0, len(list(result)))

    def test_fetch_tabix_all_results_no_pos(self):
        result, is_header = fetch_tabix_croton_predictions(
            test_tabix_file,
        )
        self.assertTrue(is_header)
        self.assertTrue(isinstance(list(result.index), list))
        self.assertEqual(9, len(list(result.index)))

    def test_fetch_tabix_all_results_no_pos(self):
        result, is_header = fetch_tabix_croton_predictions(
            test_no_header_tabix_file,
        )
        self.assertFalse(is_header)
        self.assertTrue(isinstance(list(result.index), list))
        self.assertEqual(9, len(list(result.index)))

    def test_fetch_tabix_predictions_ValueError(self):
        result, is_header = fetch_tabix_croton_predictions(
            test_tabix_file,
            chrom=2,
        )
        self.assertTrue(is_header)
        self.assertTrue(isinstance(result, list))
        self.assertEqual(0, len(result))

    def test_fetch_tabix_predictions_ValueError_no_header(self):
        result, is_header = fetch_tabix_croton_predictions(
            test_no_header_tabix_file,
            chrom=2,
        )
        self.assertFalse(is_header)
        self.assertTrue(isinstance(result, list))
        self.assertEqual(0, len(result))


class PredictionResultTestCase(TestCase):

    # chrom pos     ref_onemod3             ref_twomod3
    # chr1	69045	0.27797916531562805	0.4053782820701599
    # chr1	69045	0.28556448221206665	0.3480002284049988
    # chr1	69045	0.4016413986682892	0.3400116264820099

    ref_onemod3 = [0.27797916531562805,
                   0.28556448221206665, 0.4016413986682892]

    ref_twomod3 = [0.4053782820701599,	0.3480002284049988, 0.3400116264820099]

    # chr1	69045	A	G
    # chr1	69045	A	G
    # chr1	69045	A	G
    exp_to_A = [0, 0, 0, 0]
    exp_to_C = [0, 0, 0, 0]
    exp_to_G = [3, 0, 0, 0]
    exp_to_T = [0, 0, 0, 0]

    # HISTOGRAMS
    # chr1	69045	0.03175538778305054	0.2937266826629639	0.4079590141773224
    # chr1	69045	0.02373066544532776	0.29201287031173706	0.3529443144798279
    # chr1	69045	0.0942293107509613	0.4017235636711121	0.3391927480697632
    exp_hist1 = [3.18, 2.37, 9.42]
    exp_hist2 = [29.37, 29.2, 40.17]
    exp_hist3 = [40.80, 35.29, 33.92]

    def setUp(self):
        pos = 69045
        start = pos - 1
        end = 69072
        end = pos
        res1, is_header = fetch_tabix_croton_predictions(
            test_no_header_tabix_file,
            chrom=1,
            start=start,
            end=end,
        )
        
        '''
        try: self.response1 = build_response_json(res1, HEADER_MAP)
        except: self.response1 = build_response_json(res1, header_map_new)'''
        self.response1 = build_response_json(res1, HEADER_MAP)
        self.pp = pprint.PrettyPrinter(indent=4)
        self.pp.pprint(self.response1)

        # 2 (empty results)
        pos = 6904
        start = pos - 1
        end = 6907
        end = pos
        res2, is_header = fetch_tabix_croton_predictions(
            test_no_header_tabix_file,
            chrom=1,
            start=start,
            end=end,
        )
        print(res2)
        '''
        try: self.response2 = build_response_json(res2, HEADER_MAP)
        except: self.response2 = build_response_json(res2, header_map_new)'''
        self.response2 = build_response_json(res2, HEADER_MAP)
    def test(self):
        self.assertEqual(3, self.response1['results']['total'])

    def testEmpty(self):
        self.assertEqual(0, self.response2['results']['total'])

    def test_frequencies(self):
        print("===============test_frequencies")

        freqs = self.response1['stats']['freq']
        exp_ref_onemod3 = (statistics.mean(self.ref_onemod3))*100
        self.assertAlmostEqual(exp_ref_onemod3, freqs['ref_onemod3_mean'])

        exp_ref_twomod3 = (statistics.mean(self.ref_twomod3))*100
        self.assertAlmostEqual(exp_ref_twomod3, freqs['ref_twomod3_mean'])

        exp_no_frameshift = 100 - (exp_ref_onemod3+exp_ref_twomod3)
        self.assertAlmostEqual(exp_no_frameshift, freqs['no_frameshift'])

    def test_snp_counts(self):
        print("===============test_snp_counts")

        snp = self.response1['stats']['snp']
        self.assertEqual(self.exp_to_A, snp['toA'])
        self.assertEqual(self.exp_to_C, snp['toC'])
        self.assertEqual(self.exp_to_G, snp['toG'])
        self.assertEqual(self.exp_to_T, snp['toT'])

    def test_histograms(self):
        print("===============test_histograms")

        histograms = self.response1['stats']['histograms']
        # self.pp.pprint(histograms)

        self.assertListEqual(self.exp_hist1, histograms['alt1ins'])
        self.assertListEqual(self.exp_hist2, histograms['alt1mod3'])
        self.assertListEqual(self.exp_hist3, histograms['alt2mod3'])


class TableSummaryTestCase(TestCase):

    def setUp(self):
        pos = 69045
        start = pos - 1
        end = 69072
        res, is_header = fetch_tabix_croton_predictions(
            test_no_header_tabix_file,
            chrom=1,
            start=start,
            end=end,
        )
        header_map = HEADER_MAP
        print(res)
        res = convert_tabix_results(res, header_map=header_map)
        self.res = res

    def test_get_table_summary(self):
        ts = TableSummary(test_ecdf_file, HEADER_MAP)
        result_dict = ts.get_summary(self.res)
        print(result_dict)
        self.assertTrue(isinstance(result_dict, dict))
        keys = result_dict.keys()
        self.assertIn('frameshift', keys)
        self.assertIn('variant_impact', keys)
        self.assertIn('variant_load', keys)


class ConvertTestCase(TestCase):

    def setUp(self):
        pos = 69045
        start = pos - 1
        end = 69072
        end = pos
        self.res1, is_header = fetch_tabix_croton_predictions(
            test_no_header_tabix_file,
            chrom=1,
            start=start,
            end=end,
        )
        self.pp = pprint.PrettyPrinter(indent=4)
        self.pp.pprint(self.res1)

    def test_full_convert(self):
        predictions = convert_tabix_results(self.res1, header_map=HEADER_MAP)
        self.pp.pprint(predictions)
        self.assertEqual(3, len(predictions))

        # reversed_keys = {value: for key, value in HEADER_MAP.items()}
        reversed_dict = {v: k for k, v in HEADER_MAP.items()}
        for obj in predictions:
            self.assertEqual(reversed_dict.keys(), obj.keys())

    def test_limited_keys(self):

        # 0: "chromosome",
        # 1: "pos",
        # 2: "id",
        # 3: "ref",
        # 4: "alt",
        # 5: "diff_frameshift",

        required_keys = ["chromosome", "pos",
                         "id", "ref", "alt", "diff_frameshift"]
        predictions = convert_tabix_results(
            self.res1, header_map=HEADER_MAP, filter_keys=required_keys)
        self.pp.pprint(predictions)
        self.assertEqual(3, len(predictions))

        # reversed_keys = {value: for key, value in HEADER_MAP.items()}
        reversed_dict = {v: k for k, v in HEADER_MAP.items()}
        for obj in predictions:
            self.assertEqual(required_keys, list(obj.keys()))

    def test_limited_keys_and_values(self):

        # 0: "chromosome",
        # 1: "pos",
        # 2: "id",
        # 3: "ref",
        # 4: "alt",
        # 5: "diff_frameshift",
        full_predictions = convert_tabix_results(
            self.res1, header_map=HEADER_MAP)
        required_keys = ["chromosome", "pos",
                         "id", "ref", "alt", "diff_frameshift", "alt_1ins"]

        required_keys = ["chromosome", "pos", "id", "ref", "alt", "pamid", "start", "end", "strand", "ref_seq", "alt_seq",
                         "abs_diff", "ref_1ins", "alt_1ins", "diff_1ins", "ref_onemod3", "alt_onemod3", "diff_onemod3", "ref_twomod3", "alt_twomod3", "diff_twomod3", "ref_frameshift", "alt_frameshift", "diff_frameshift"]

        limited_predictions = convert_tabix_results(
            self.res1, header_map=HEADER_MAP, filter_keys=required_keys)
        self.assertEqual(3, len(limited_predictions))

        # reversed_keys = {value: for key, value in HEADER_MAP.items()}
        reversed_dict = {v: k for k, v in HEADER_MAP.items()}
        for idx, pred in enumerate(limited_predictions):
            self.assertEqual(set(required_keys), set(pred.keys()))
            for key in required_keys:
                self.assertEqual(pred[key], full_predictions[idx][key])
