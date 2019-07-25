# Purpose: Store scores for the locus so they can be passed to the ensemble
# method and/or referened later

class ScoreContainer:
    def __init__(self):
        self.finemap_clpp = None
        self.finemap_clpp_mod = None
        self.coloc_h4 = None
        self.rtc_neg_log_pval = None
        self.smr_neg_log_pval = None
        self.gsmr_neg_log_pval = None
        self.twas_neg_log_pval = None
        self.baseline_neg_log_pval = None
        self.smart_baseline_neg_log_pval = None
        self.ensemble_score = None
