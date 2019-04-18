# Author: Mike Gloudemans
#
# Class: TestLocus
# A locus to test for colocalization, along with all of the
# relevant data needed to conduct the tests.
#

import finemap
import ecaviar
import caviarbf
import coloc
import rtc
import twas_real_data
import plot_loci as plot
import math
import baseline

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

    # Run colocalization tests. Which ones to do will depend
    # on the settings file.
    def run(self):

        #print "Analyzing {0} {6} {1} {2} {3} {4} {5}".format(self.gwas_suffix, self.eqtl_suffix, self.gene, self.chrom, self.pos, self.pval, self.trait)

        # Note: Might eventually create a wrapper function for finemap and ecaviar that dispatches the two 
        # as necessary, depending on overlaps.

        plotworthy = False

        if "finemap" in self.settings["methods"]:
            clpp = finemap.run_finemap(self)

            if not isinstance(clpp, basestring) and clpp > 0.3: 
                plotworthy = True

            # NOTE: Temporary; for debugging
            if isinstance(clpp, basestring):
                print "FAILED Analyzing {0} {6} {1} {2} {3} {4} {5} {7}".format(self.gwas_suffix, self.eqtl_suffix, self.gene, self.chrom, self.pos, self.pval, self.trait, clpp)

        if "ecaviar" in self.settings["methods"]:
            clpp = ecaviar.run_ecaviar(self)
            if clpp > 0.02: 
                plotworthy = True

        if "caviarbf" in self.settings["methods"]:
            clpp = caviarbf.run_caviarbf(self)
            if clpp > 0.02:
                plotworthy = True

        if "finemap" in self.settings["methods"] or "ecaviar" in self.settings["methods"] or "caviarbf" in self.settings["methods"]:
            finemap.purge_tmp_files(self)

        if "coloc" in self.settings["methods"]:
            h4pp = coloc.run_coloc(self)
            if h4pp > 0.5:
                plotworthy = True

        if "rtc" in self.settings["methods"]:
            rtc_score = rtc.run_rtc(self)
            if rtc_score > 0.8:
                plotworthy = True

        if "twas" in self.settings["methods"]:
            twas_p = twas_real_data.run_twas(self)
            if twas_p > 5:
                plotworthy = True

        if "baseline" in self.settings["methods"]:
            pval = baseline.run_baseline(self)

            print pval

        ''' 
        if "enloc" in self.settings.keys():
            run_enloc(basedir, data, settings)

        if "pics" in self.settings.keys():
            run_pics(basedir, data, settings)

        if "smr" in self.settings.keys():
            run_smr(basedir, data, settings)

        if "gwaspw" in self.settings.keys():
            run_gwaspw(basedir, data, settings)

        '''
        #plotworthy=True

        # Plot the result if it's significant.
        if (plotworthy or self.settings["plot_all"] == True) and self.settings["plot_none"] == False:
            plot.locus_compare(self)

