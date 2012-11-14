'''
Created on 05.10.2012

@author: mhoefli
'''

import unittest, tempfile, shutil, os
import FRETUtils.Ensemble as Ensemble

def createInvalidProbabilityClassFile(fname, classes, prbs):
    with open(fname, "w") as f:
        f.write("\n\n")
        for i in range(len(classes)):
            f.write(".*%s.* %s %6.4f\n" % (classes[i], classes[i], prbs[i]))
        f.write("a b c 0.7\n")

class testEnsembleFunctions(unittest.TestCase):
    def setUp(self):
        self.workdir = tempfile.mkdtemp()
        self.prevdir = os.curdir
        os.chdir(self.workdir)
        createInvalidProbabilityClassFile("broken.prb", ("all", "none"), (0.5, 0.5))

    def tearDown(self):
        os.chdir(self.prevdir)
        shutil.rmtree(self.workdir)

    def test_readInvalidEnsembleFile(self):
        self.assertRaises(ValueError, Ensemble.readProbabilities, "broken.prb")
