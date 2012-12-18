'''
Created on 18.12.2012

@author: mhoefli
'''
import unittest
from FRETUtils.TransferMatrices import GlobalAVGKappaTransferMatrix, DistanceAVGKappaTransferMatrix, DistanceKappaTransferMatrix
from FRETUtils.Helpers import rToEff, effToR
import tempfile, os, shutil
import numpy

def getConstantBurstsize(burstsize):
    return burstsize


class TransferMatrixTests(unittest.TestCase):

    def setUp(self):
        self.constant500Burstgen = lambda :getConstantBurstsize(500)


    def tearDown(self):
        pass


    def testGlobal(self):
        tm = GlobalAVGKappaTransferMatrix(20, 11, 5, self.constant500Burstgen, 5.475, [5, 6])
        self.assertEqual(tm.getMatrix().shape, (20, 11))
        self.assertEqual(tm.getMatrix()[9][5], 1.)

    def testDistanceAVG(self):
        R = numpy.linspace(5, 6, 1100)
        kappa = [ 2. / 3 ] * 1100
        weights = [ 1. ] * 1100
        tm = DistanceAVGKappaTransferMatrix(20, 11, 5, self.constant500Burstgen, 5.475, R, kappa, weights)
        tm.generateMatrix()
        self.assertEqual(tm.RRange[0], 5)
        self.assertEqual(tm.RRange[1], 6)

    def testDistanceAVGLocalKappavsGlobal(self):
        R = numpy.linspace(5.05, 5.95, 20)
        kappa = [ 2. / 3 ] * 20
        weights = [ 1. ] * 20
        tm = DistanceAVGKappaTransferMatrix(20, 11, 5, self.constant500Burstgen, 5.475, R, kappa, weights, RRange = (5, 6))
        tm.generateMatrix()
        self.assertEqual(tm.RRange[0], 5)
        self.assertEqual(tm.RRange[1], 6)
        print tm.getMatrix()
        self.assertEqual(tm.getMatrix()[9][5], 1.)
        tmref = GlobalAVGKappaTransferMatrix(20, 11, 5, self.constant500Burstgen, 5.475, [5, 6])
        self.assertAlmostEqual((tm.getMatrix() - tmref.getMatrix()).sum(), 0.0, delta = 1e-10)

    def testPlotting(self):
        tm = GlobalAVGKappaTransferMatrix(20, 11, 5, self.constant500Burstgen, 5.475, [5, 6])
        odir = tempfile.mkdtemp(suffix = "plottest")
        ofl = os.path.join(odir, "plot.png")
        tm.plotToFile(ofl)
        self.assertEqual(os.path.exists(ofl), True)
        shutil.rmtree(odir)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
