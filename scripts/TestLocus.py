# Author: Mike Gloudemans
#
# Class: TestLocus
# A locus to test for colocalization, along with all of the
# relevant data needed to conduct the tests.
#

import finemap
import ecaviar
import coloc
import plot_loci as plot

class TestLocus:
    def __init__(self, data, settings, basedir, gene, snp, gwas_file, eqtl_file):
        self.basedir = basedir
        self.data = data
        self.settings = settings
        self.gene = gene
        self.chrom = snp.chrom
        self.pos = snp.pos
        self.gwas_suffix = gwas_file.split("/")[-1].replace(".", "_")
        self.eqtl_suffix = eqtl_file.split("/")[-1].replace(".", "_")
        

    # Run colocalization tests. Which ones to do will depend
    # on the settings file.
    def run(self):

        print "Analyzing {0} {1} {2} {3} {4}".format(self.gwas_suffix, self.eqtl_suffix, self.gene, self.chrom, self.pos)

        # Note: Might eventually create a wrapper function for finemap and ecaviar that dispatches the two 
        # as necessary, depending on overlaps.

        plotworthy = False

        if "finemap" in self.settings["methods"]:
            clpp = finemap.run_finemap(self)

            # TODO: Increase threshold for significance here.
            if clpp > 0.001: 
                plotworthy = True


        if "ecaviar" in self.settings["methods"]:
            clpp = ecaviar.run_ecaviar(self)

            if clpp > 0.001: 
                plotworthy = True

        if "finemap" in self.settings["methods"] or "ecaviar" in self.settings["methods"]:
            finemap.purge_tmp_files(self)

        if "coloc" in self.settings["methods"]:
            coloc.run_coloc(self)

        ''' 

        if "enloc" in self.settings.keys():
            run_enloc(basedir, data, settings)

        if "pics" in self.settings.keys():
            run_pics(basedir, data, settings)

        if "smr" in self.settings.keys():
            run_smr(basedir, data, settings)

        if "gwaspw" in self.settings.keys():
            run_gwaspw(basedir, data, settings)

        if "twas" in self.settings.keys():
            run_twas(basedir, data, settings)
        '''

        # Plot the result if it's significant.
        if plotworthy:
            plot.locus_zoom_plot(self, clpp)
            plot.pvalue_plot(self, clpp)


