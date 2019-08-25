# Author: Mike Gloudemans
#
# Class: TestLocus A locus to test for colocalization, along with all of the
# relevant data needed to conduct the tests.

import finemap
import traceback
import ecaviar
import caviarbf
import coloc
import rtc
import twas
import plot_loci as plot
import math
import baseline
import smr
import gsmr
import ensemble
from ScoreContainer import ScoreContainer

class TestLocus:
    def __init__(self, data, settings, basedir, tmpdir, gene, snp, gwas_file, eqtl_file, trait):
        self.basedir = basedir
        self.data = data
        self.settings = settings
        self.gene = gene
        self.chrom = snp.chrom
        self.pos = snp.pos
        self.pval = snp.pval
        self.gwas_suffix = gwas_file.split("/")[-1].replace(".", "_")
        self.eqtl_suffix = eqtl_file.split("/")[-1].replace(".", "_")
        self.gwas_file = gwas_file
        self.eqtl_file = eqtl_file
        self.tmpdir = tmpdir
        self.conditional_level = 0      # Currently serves no purpose, but may be implemented later
        self.trait = trait
        self.scores = ScoreContainer()

    # Run colocalization tests. Which ones to do will depend
    # on the settings file.
    def run(self):

        # Note: Might eventually create a wrapper function for finemap and ecaviar that dispatches the two 
        # as necessary, depending on overlaps.

        plotworthy = False

        if "finemap" in self.settings["methods"]:
            try:
                clpps = finemap.run_finemap(self)

                if not isinstance(clpps, basestring):
                    clpp, clpp_mod = clpps
                    if clpp_mod > 0.3: 
                        plotworthy = True
                    
                    self.scores.finemap_clpp = clpp
                    self.scores.finemap_clpp_mod = clpp_mod

            except Exception as e:
                error = str(e)
                error = error + "\t" + traceback.format_exc().replace("\n", "NEWLINE").replace("\t", "TAB")
                with open("{0}/ERROR_variants.txt".format(self.basedir),"a") as a:
                    a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tfinemap\n".format(self.gwas_file, self.eqtl_file, self.chrom, self.pos, self.gene, self.trait, error))

        if "ecaviar" in self.settings["methods"]:
            try:
                clpp = ecaviar.run_ecaviar(self)
                if clpp > 0.02: 
                    plotworthy = True

            except Exception as e:
                error = str(e)
                error = error + "\t" + traceback.format_exc().replace("\n", "NEWLINE").replace("\t", "TAB")
                with open("{0}/ERROR_variants.txt".format(self.basedir),"a") as a:
                    a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tecaviar\n".format(self.gwas_file, self.eqtl_file, self.chrom, self.pos, self.gene, self.trait, error))

        if "caviarbf" in self.settings["methods"]:
            try:
                clpp = caviarbf.run_caviarbf(self)
                if clpp > 0.02:
                    plotworthy = True

            except Exception as e:
                error = str(e)
                error = error + "\t" + traceback.format_exc().replace("\n", "NEWLINE").replace("\t", "TAB")
                with open("{0}/ERROR_variants.txt".format(self.basedir),"a") as a:
                    a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tcaviarbf\n".format(self.gwas_file, self.eqtl_file, self.chrom, self.pos, self.gene, self.trait, error))

        # This isn't a very clean way to purge shared temporary files, 
        # though technically it should usually work? Fix this later
        if "finemap" in self.settings["methods"] or "ecaviar" in self.settings["methods"] or "caviarbf" in self.settings["methods"]:
            if not self.settings["debug"]:
                finemap.purge_tmp_files(self)
        
        if "coloc" in self.settings["methods"]:
            try:
                h4pp = coloc.run_coloc(self)
                if h4pp > 0.5:
                    plotworthy = True

                self.scores.coloc_h4 = h4pp

            except Exception as e:
                error = str(e)
                error = error + "\t" + traceback.format_exc().replace("\n", "NEWLINE").replace("\t", "TAB")
                with open("{0}/ERROR_variants.txt".format(self.basedir),"a") as a:
                    a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tcoloc\n".format(self.gwas_file, self.eqtl_file, self.chrom, self.pos, self.gene, self.trait, error))

        if "rtc" in self.settings["methods"]:
            try:
                rtc_score = rtc.run_rtc(self)
                if rtc_score > 0.8:
                    plotworthy = True

                self.scores.rtc_neg_log_pval = rtc_score

            except Exception as e:
                error = str(e)
                error = error + "\t" + traceback.format_exc().replace("\n", "NEWLINE").replace("\t", "TAB")
                with open("{0}/ERROR_variants.txt".format(self.basedir),"a") as a:
                    a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\trtc\n".format(self.gwas_file, self.eqtl_file, self.chrom, self.pos, self.gene, self.trait, error))

        if "twas" in self.settings["methods"]:
            try:
                twas_p = twas.run_twas(self)
                if twas_p > 5:
                    plotworthy = True

                self.scores.twas_neg_log_pval = twas_p

            except Exception as e:
                error = str(e)
                error = error + "\t" + traceback.format_exc().replace("\n", "NEWLINE").replace("\t", "TAB")
                with open("{0}/ERROR_variants.txt".format(self.basedir),"a") as a:
                    a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\ttwas\n".format(self.gwas_file, self.eqtl_file, self.chrom, self.pos, self.gene, self.trait, error))

        if "baseline" in self.settings["methods"]:
            try:
                pval, pval5 = baseline.run_baseline(self)

                self.scores.baseline_neg_log_pval = pval
                self.scores.smart_baseline_neg_log_pval = pval5

            except Exception as e:
                error = str(e)
                error = error + "\t" + traceback.format_exc().replace("\n", "NEWLINE").replace("\t", "TAB")
                with open("{0}/ERROR_variants.txt".format(self.basedir),"a") as a:
                    a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tbaseline\n".format(self.gwas_file, self.eqtl_file, self.chrom, self.pos, self.gene, self.trait, error))

        if "smr" in self.settings["methods"]:
            try:
                pval = smr.run_smr(self)
                self.scores.smr_neg_log_pval = pval

            except Exception as e:
                error = str(e)
                error = error + "\t" + traceback.format_exc().replace("\n", "NEWLINE").replace("\t", "TAB")
                with open("{0}/ERROR_variants.txt".format(self.basedir),"a") as a:
                    a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tsmr\n".format(self.gwas_file, self.eqtl_file, self.chrom, self.pos, self.gene, self.trait, error))

        if "gsmr" in self.settings["methods"]:
            try:
                pval = gsmr.run_gsmr(self)
                self.scores.gsmr_neg_log_pval = pval

            except Exception as e:
                error = str(e)
                error = error + "\t" + traceback.format_exc().replace("\n", "NEWLINE").replace("\t", "TAB")
                with open("{0}/ERROR_variants.txt".format(self.basedir),"a") as a:
                    a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tgsmr\n".format(self.gwas_file, self.eqtl_file, self.chrom, self.pos, self.gene, self.trait, error))

        ''' 
        if "enloc" in self.settings.keys():
            try:
                enloc = enloc.run_enloc(self)

            except Exception as e:
                error = str(e)
                error = error + "\t" + traceback.format_exc().replace("\n", "NEWLINE").replace("\t", "TAB")
                with open("{0}/ERROR_variants.txt".format(self.basedir),"a") as a:
                    a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tenloc\n".format(self.gwas_file, self.eqtl_file, self.chrom, self.pos, self.gene, self.trait, error))
        '''

        if "ensemble" in self.settings["methods"]:
            try:
                score = ensemble.run_ensemble(self)
                self.scores.ensemble_score = score

            except Exception as e:
                error = str(e)
                error = error + "\t" + traceback.format_exc().replace("\n", "NEWLINE").replace("\t", "TAB")
                with open("{0}/ERROR_variants.txt".format(self.basedir),"a") as a:
                    a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tensemble\n".format(self.gwas_file, self.eqtl_file, self.chrom, self.pos, self.gene, self.trait, error))


        # Plot the result if it's significant.
        if (plotworthy or self.settings["plot_all"] == True) and ("plot_none" not in self.settings or self.settings["plot_none"] != True):
            plotted = plot.locus_compare(self)

            if isinstance(plotted, basestring):
                with open("{0}/unplotted_variants.txt".format(self.self.basedir),"a") as a:
                    a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(self.gwas_file, self.eqtl_file, self.chrom, self.pos, self.gene, plotted, self.trait))


