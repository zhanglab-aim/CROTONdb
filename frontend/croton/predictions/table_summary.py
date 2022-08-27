"""Generate a summary text card for any given genomic interval of CROTON ref/alt predictions
"""

import os
import pandas as pd
import numpy as np
import scipy.stats as ss
import pickle
from django.conf import settings
from utils.hostname import timeit


#ECDF = pickle.load(open('/mnt/ceph/users/zzhang/croton/data/052121-variant-2/gw_data.pkl', 'rb'))


class TableSummary:
    # class namespace CONSTANTS
    REF_PRED_POSITIONS = [
        20,  # ref 1ins
        26,  # ref onemod3
        29,  # ref twomod3
        32,  # ref frameshift
    ]
    DIFF_PRED_POSITIONS = [
        22,  # diff 1ins
        28,  # diff onemod3
        31,  # diff twomod3
        34,  # diff frameshift
    ]

    REF_PRED_POSITIONS_new = [
        'ref_1ins',
        'ref_onemod3',
        'ref_twomod3',
        'ref_frameshift'
    ]

    DIFF_PRED_POSITIONS_new = [
        'diff_1ins',
        'diff_onemod3',
        'diff_twomod3',
        'diff_frameshift'
    ]

    SUMMARY_SEPARATOR = '\n----------\n'

    def __init__(self, ecdf_fp, header_map):
        assert os.path.isfile(ecdf_fp)
        self.ecdf = pickle.load(open(ecdf_fp, 'rb'))
        self.header_map = header_map

    @staticmethod
    def _unify_table(res):
        if type(res) is list:
            res_ = pd.DataFrame.from_dict(res, orient='columns')
        elif type(res) is pd.DataFrame:
            res_ = res.copy()
        else:
            raise TypeError(
                'cannot convert type for query table, got %s' % type(res))
        return res_

    @timeit
    def get_table_summary(self, res):
        res = self._unify_table(res)
        ref_preds, ref_cdfs = self.refpred_percentile(res)
        vareff_preds, vareff_cdfs = self.vareff_percentile(res)
        obs_varcnt, exp_varcnt, varcnt_cdf = self.varcnt_percentile(res)
        s1 = self.pred_score_summary(
            pred=ref_preds['ref_frameshift'], cdf=ref_cdfs['ref_frameshift'], label='Frameshift')
        s3 = self.pred_score_summary(
            pred=vareff_preds['diff_frameshift'], cdf=vareff_cdfs['diff_frameshift'], label='Variant impact to frameshift')
        s2 = self.varcnt_summary(obs_varcnt, exp_varcnt, varcnt_cdf)
        return self.SUMMARY_SEPARATOR.join([s1, s2, s3])

    @staticmethod
    def cdf_to_label(ecdf_value):
        if ecdf_value > 0.9:
            return 'Substantially Higher'
        elif ecdf_value > 0.7:
            return 'Higher'
        elif ecdf_value <= 0.7 and ecdf_value > 0.3:
            return 'Typical'
        elif ecdf_value > 0.1:
            return 'Lower'
        else:
            return 'Substantially Lower'

    @staticmethod
    def category(ecdf_value):
        if ecdf_value > 0.9:
            return 'substantially higher'
        elif ecdf_value > 0.7:
            return 'higher'
        elif ecdf_value <= 0.7 and ecdf_value > 0.3:
            return 'typical'
        elif ecdf_value > 0.1:
            return 'lower'
        else:
            return 'substantially lower'

    @staticmethod
    def pred_score_summary(pred, cdf, label):
        """
        Paramaters
        ----------
        pred : float
            a single prediction score for the query interval
        cdf : float
            the percentile in the genome-wide prediction scores given the prediction
        label : str
            the string label for what the score is
        """
        str_formatter = "The query interval has {category} {label}\nPredicted: {" \
                        "score}\nPercentile over genome: {ecdf} "
        category = TableSummary.cdf_to_label(cdf)
        s = str_formatter.format(
            category=category,
            score='%.2f%%' % (pred*100),
            ecdf='%.2f%%' % (cdf*100),
            label=label
        )
        return s

    @staticmethod
    def varcnt_summary(obs_varcnt, exp_varcnt, varcnt_cdf):
        str_formatter = "Variant load for this query is {category}\nObserved variants : {obs}\nExpected variants : {exp}"
        category = TableSummary.cdf_to_label(varcnt_cdf)
        s = str_formatter.format(
            category=category,
            obs=obs_varcnt,
            exp=int(exp_varcnt)
        )
        return s

    @timeit
    def get_summary(self, res):
        result_dict = {}
        res = self._unify_table(res)
        ref_preds, ref_cdfs = self.refpred_percentile(res)

        # FRAMESHIFT
        result_dict["frameshift"] = {}
        result_dict["frameshift"]["predicted"] = '%.2f%%' % (
            ref_preds['ref_frameshift']*100)
        result_dict["frameshift"]["perc_over_genome"] = '%.2f%%' % (
            ref_cdfs['ref_frameshift']*100)
        result_dict["frameshift"]["category"] = TableSummary.cdf_to_label(
            ref_cdfs['ref_frameshift'])

        # VARIANT impact to frameshift
        vareff_preds, vareff_cdfs = self.vareff_percentile(res)
        result_dict["variant_impact"] = {}
        result_dict["variant_impact"]["predicted"] = '%.2f%%' % (
            vareff_preds['diff_frameshift']*100)
        result_dict["variant_impact"]["perc_over_genome"] = '%.2f%%' % (
            vareff_cdfs['diff_frameshift']*100)
        result_dict["variant_impact"]["category"] = TableSummary.cdf_to_label(
            vareff_cdfs['diff_frameshift'])

        # variant load
        obs_varcnt, exp_varcnt, varcnt_cdf = self.varcnt_percentile(res)
        result_dict["variant_load"] = {}
        result_dict["variant_load"]["observed"] = obs_varcnt
        result_dict["variant_load"]["expected"] = int(exp_varcnt)
        result_dict["variant_load"]["category"] = TableSummary.cdf_to_label(
            varcnt_cdf)

        return result_dict

    def get_pred_percentile(self, res, pred_pos, ecdf_key, take_abs=False):
        try: 
            pred_idx = [self.header_map[i] for i in pred_pos]
            interval_preds = res[['pamid'] + pred_idx].drop_duplicates()
            if take_abs is True:
                preds = interval_preds[pred_idx].astype(
                    'float').abs().mean().to_dict()
            else:
                preds = interval_preds[pred_idx].astype('float').mean().to_dict()
            preds_cdf = {k: self.ecdf[ecdf_key]['%.3f' % v][k]
                        for k, v in preds.items()}
        except:
            preds = preds_cdf = {'ref_1ins':1, 'ref_onemod3':1, 'ref_twomod3':1, 'ref_frameshift':1, 'diff_1ins': 1, 'diff_onemod3':1, 'diff_twomod3':1, 'diff_frameshift':1} 

        return preds, preds_cdf

    def refpred_percentile(self, res):
        try: 
            preds, preds_cdf = self.get_pred_percentile(
                res, pred_pos=self.REF_PRED_POSITIONS, ecdf_key='ref_cdf')
        except:
            preds, preds_cdf = self.get_pred_percentile(
                res, pred_pos=self.REF_PRED_POSITIONS_new, ecdf_key='ref_cdf')
        return preds, preds_cdf

    def vareff_percentile(self, res):
        try: 
            preds, preds_cdf = self.get_pred_percentile(
                res, pred_pos=self.DIFF_PRED_POSITIONS, ecdf_key='diff_cdf', take_abs=True)
        except:
            preds, preds_cdf = self.get_pred_percentile(
                res, pred_pos=self.DIFF_PRED_POSITIONS_new, ecdf_key='diff_cdf', take_abs=True)
        return preds, preds_cdf

    def varcnt_percentile(self, res):
        obs_varcnt = res.shape[0]
        uniq_df = res[['pamid', 'start', 'end']].drop_duplicates()
        obs_len = (uniq_df['end'].astype('int') -
                   uniq_df['start'].astype('int')).sum()
        exp_varcnt = obs_len * self.ecdf['perbp_varcnt_cdf']
        cdf = ss.poisson.cdf(obs_varcnt, mu=exp_varcnt)
        return obs_varcnt, exp_varcnt, cdf
