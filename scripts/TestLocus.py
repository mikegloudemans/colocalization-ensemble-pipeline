# Author: Mike Gloudemans
#
# Class: TestLocus
# A locus to test for colocalization, along with all of the
# relevant data needed to conduct the tests.
#

import finemap
import plot_loci as plot

class TestLocus:
    def __init__(self, data, settings, basedir, gene, snp, gwas_file, eqtl_file):
        self.basedir = basedir
        self.data = data
        self.settings = settings
        self.gene = gene
        self.snp = snp
        self.gwas_suffix = gwas_file.split("/")[-1].replace(".", "_")
        self.eqtl_suffix = eqtl_file.split("/")[-1].replace(".", "_")
        

    # Run colocalization tests. Which ones to do will depend
    # on the settings file.
    def run(self):

        print "Analyzing {0} {1} {2} {3} {4}".format(self.gwas_suffix, self.eqtl_suffix, self.gene, self.snp[0], self.snp[1])

        if "finemap" in self.settings["methods"]:
            clpp = finemap.run_finemap(self)

        # Plot the result if it's significant.
        # TODO: Increase threshold for significance here.
        if clpp > 0.00001: 
            plot.locus_zoom_plot(self, clpp)
            plot.pvalue_plot(self, clpp)


        '''
        if "ecaviar" in self.settings.keys():
            pass    

        if "coloc" in self.settings.keys():
            run_coloc(basedir, data, settings)

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
