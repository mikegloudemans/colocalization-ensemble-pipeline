# Author: Mike Gloudemans
#
# Class: TestLocus
# A locus to test for colocalization, along with all of the
# relevant data needed to conduct the tests.
#

class TestLocus:
    def __init__(self, data, settings, basedir):
        self.basedir = basedir
        self.data = data
        self.settings = settings

    # Run colocalization tests. Which ones to do will depend
    # on the settings file.
    def run(self):
        if "finemap" in self.settings.keys():
            run_finemap(basedir, data, settings)

        if "ecaviar" in self.settings.keys():
            run_ecaviar(basedir, data, settings)

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

