'''
Created on 18.12.2012

@author: mhoefli
'''

from __future__ import division

import unittest
from FRETUtils.TransferMatrices import TransferMatrix, GlobalAVGKappaTransferMatrix, DistanceAVGKappaTransferMatrix, DistanceKappaTransferMatrix, genEffIndex, generateBinMid, modifyR0
from FRETUtils.Helpers import rToEff, genRandomVec, getKappa
from FRETUtils.Efficiencies import efficiencyDeltaFix
import tempfile, os, shutil, random
import numpy
from math import sqrt
from matplotlib import pyplot

def matDiff(mat1, mat2):
    diff = sqrt(((mat1 - mat2) ** 2).max())
    return diff

def getConstantBurstsize(burstsize):
    return burstsize

def getShotNoise(myeff, effbins, samples, samplesize):
    effs = numpy.zeros(effbins)
    for _ in range(samples):
        rnds = numpy.random.random(samplesize)
        accs = (rnds < myeff).sum()
        eff = efficiencyDeltaFix((accs / samplesize,))
        ndx = genEffIndex(eff[0], effbins)
        effs[ndx] += 1 / samples
    return effs


class TransferMatrixTests(unittest.TestCase):

    def setUp(self):
        self.constant5000Burstgen = lambda :getConstantBurstsize(5000)
        self.constant5Burstgen = lambda :getConstantBurstsize(5)
        self.constant50Burstgen = lambda :getConstantBurstsize(50)
        self.constant500Burstgen = lambda :getConstantBurstsize(500)

    def tearDown(self):
        pass

    def assertMatrixAlmostEqual(self, m1, m2, delta = 0.05):
        diff = matDiff(m1, m2)
        mmax = m1.max()
        if m2.max() > max:
            mmax = m2.max()
        diff /= mmax
        print "Matrix difference is %f percent" % (diff * 100)
        print "Matrix1 sum", m1.sum()
        print "Matrix1 sum", m2.sum()
        if diff > delta:
            pyplot.figure()
            pyplot.subplot(221)
            pyplot.imshow(m1.T, origin = "lower", interpolation = 'nearest', aspect = "auto")
            pyplot.title("Matrix1")
            pyplot.subplot(222)
            pyplot.imshow(m2.T, origin = "lower", interpolation = 'nearest', aspect = "auto")
            pyplot.title("Matrix2")
            pyplot.subplot(223)
            pyplot.imshow((numpy.sqrt((m1 - m2) ** 2)).T, origin = "lower", interpolation = 'nearest', aspect = "auto")
            pyplot.title("abs differences")
            pyplot.subplot(224)
            pyplot.imshow((((m1 - m2))).T, origin = "lower", interpolation = 'nearest', aspect = "auto")
            pyplot.title("differences")

            pyplot.show()

        self.assertAlmostEqual(diff, 0.0, delta = delta)


    def testEfficiencyIndexForEff1(self):
        self.assertEqual(genEffIndex(1., 10), 9)

    def testAbstract(self):
        tm = TransferMatrix(20, 11, [5, 6])
        self.assertRaises(NotImplementedError, tm.generateMatrix)

    def testSamplesOutsideRange(self):
        R = numpy.linspace(4.95, 6.05, 22)
        kappa = [ 2. / 3 ] * 22
        weights = [ 1. ] * 22
        tm = DistanceAVGKappaTransferMatrix(20, 11, 200, self.constant5000Burstgen, 5.475, R, kappa, weights, RRange = (5., 6))
        tm.generateMatrix()
        self.assertEqual(tm.RRange[0], 5)
        self.assertEqual(tm.RRange[1], 6)
        self.assertAlmostEqual(tm.getMatrix()[9][5], 1., delta = 0.01)
        tmref = GlobalAVGKappaTransferMatrix(20, 11, 200, self.constant5000Burstgen, 5.475, (5, 6))
        self.assertMatrixAlmostEqual(tm.getMatrix() , tmref.getMatrix(), delta = 0.15)

    def testBinWithZeroSamples(self):
        R = numpy.linspace(4.95, 6.05, 22)
        kappa = [ 2. / 3 ] * 22
        weights = [ 1. ] * 22
        tm = DistanceAVGKappaTransferMatrix(20, 11, 5, self.constant5000Burstgen, 5.475, R, kappa, weights, RRange = (5., 7))
        self.assertRaises(ValueError, tm.generateMatrix)
        tm = DistanceKappaTransferMatrix(20, 11, 5, self.constant5000Burstgen, 5.475, R, kappa, weights, RRange = (5., 7))
        self.assertRaises(ValueError, tm.generateMatrix)

    def testGlobal(self):
        tm = GlobalAVGKappaTransferMatrix(20, 11, 5, self.constant5000Burstgen, 5.475, (5, 6))
        self.assertEqual(tm.getMatrix().shape, (20, 11))
        self.assertEqual(tm.RRange[0], 5)
        self.assertEqual(tm.RRange[1], 6)
        self.assertEqual(tm.getMatrix()[9][5], 1.)

    def testGlobalShotNoise(self):
        for sn, sng in ((5, self.constant5Burstgen), (50, self.constant50Burstgen), (500, self.constant500Burstgen), (5000, self.constant5000Burstgen)):
            tm = GlobalAVGKappaTransferMatrix(20, 11, 1000, sng, 5.475, (5, 6))
            tmx = tm.getMatrix()
            print "Burstsizes %d" % sn
            for _ in range(7):
                mbin = random.choice(range(20))
                binmid = generateBinMid(1, mbin, 20, 5)
                eff = rToEff(binmid, R0 = 5.475)
                print "Testing mbin %d at pos %f with efficiency %f" % (mbin, binmid, eff)
                shotnoise = getShotNoise(eff, 11, 1000, sn)
                tmxvec = tmx[mbin, :]
                tmxvec.shape = (11, 1)
                self.assertAlmostEqual((tmxvec - shotnoise).sum(), 0.0, delta = 0.005)


    def testDistanceAVG(self):
        R = numpy.linspace(5, 6, 1100)
        kappa = [ 2. / 3 ] * 1100
        weights = [ 1. ] * 1100
        tm = DistanceAVGKappaTransferMatrix(20, 11, 5, self.constant5000Burstgen, 5.475, R, kappa, weights)
        tm.generateMatrix()
        self.assertEqual(tm.getMatrix().shape, (20, 11))
        self.assertEqual(tm.RRange[0], 5)
        self.assertEqual(tm.RRange[1], 6)
        self.assertEqual(tm.getMatrix()[9][5], 1.)

    def testDistanceAVGAlteredKappa(self):
        R = numpy.linspace(5, 6, 1100)
        kappa = [ 2. / 3 ] * 200 + [ 1. / 3 ] * 400 + [ 2. / 3 ] * 500
        weights = [ 1. ] * 1100
        tm = DistanceAVGKappaTransferMatrix(20, 11, 5, self.constant5000Burstgen, 5.475, R, kappa, weights)
        tm.generateMatrix()
        self.assertEqual(tm.getMatrix().shape, (20, 11))
        self.assertEqual(tm.RRange[0], 5)
        self.assertEqual(tm.RRange[1], 6)
        R0neu = modifyR0(5.475, 1. / 3)
        print "New R0 is %f" % R0neu
        eff = rToEff(5.475, R0 = R0neu)
        effndx = int(eff * 11)
        self.assertAlmostEqual(tm.getMatrix()[9][effndx], 1, delta = 0.01)

    def testDistanceAVGShotNoise(self):
        R = numpy.linspace(5, 6, 1100)
        kappa = [ 2. / 3 ] * 1100
        weights = [ 1. ] * 1100
        for sn, sng in ((5, self.constant5Burstgen), (50, self.constant50Burstgen), (500, self.constant500Burstgen), (5000, self.constant5000Burstgen)):
            tm = DistanceAVGKappaTransferMatrix(20, 11, 1000, sng, 5.475, R, kappa, weights, RRange = (5, 6))
            tmx = tm.getMatrix()
            print "Burstsizes %d" % sn
            for _ in range(7):
                mbin = random.choice(range(20))
                binmid = generateBinMid(1, mbin, 20, 5)
                eff = rToEff(binmid, R0 = 5.475)
                print "Testing mbin %d at pos %f with efficiency %f" % (mbin, binmid, eff)
                shotnoise = getShotNoise(eff, 11, 1000, sn)
                tmxvec = tmx[mbin, :]
                tmxvec.shape = (11, 1)
                self.assertAlmostEqual((tmxvec - shotnoise).sum(), 0.0, delta = 0.005)

    def testDistanceKappa(self):
        R = numpy.linspace(5, 6, 1100)
        kappa = [ 2. / 3 ] * 1100
        weights = [ 1. ] * 1100
        tm = DistanceKappaTransferMatrix(20, 11, 5, self.constant5000Burstgen, 5.475, R, kappa, weights)
        tm.generateMatrix()
        self.assertEqual(tm.getMatrix().shape, (20, 11))
        self.assertEqual(tm.RRange[0], 5)
        self.assertEqual(tm.RRange[1], 6)
        self.assertEqual(tm.getMatrix()[9][5], 1.)

    def testDistanceAVGLocalKappavsGlobal(self):
        R = numpy.linspace(5.05, 5.95, 20)
        kappa = [ 2. / 3 ] * 20
        weights = [ 1. ] * 20
        tm = DistanceAVGKappaTransferMatrix(20, 11, 5000, self.constant5000Burstgen, 5.475, R, kappa, weights, RRange = (5, 6))
        tm.generateMatrix()
        self.assertEqual(tm.RRange[0], 5)
        self.assertEqual(tm.RRange[1], 6)
        self.assertAlmostEqual(tm.getMatrix()[9][5], 1., delta = 0.01)
        tmref = GlobalAVGKappaTransferMatrix(20, 11, 5000, self.constant5000Burstgen, 5.475, [5, 6])
        self.assertMatrixAlmostEqual(tm.getMatrix(), tmref.getMatrix(), delta = 0.05)

    def testGlobalVsGlobal(self):
        R = numpy.linspace(5.05, 5.95, 20)
        kappa = [ 2. / 3 ] * 20
        weights = [ 1. ] * 20
        tm = DistanceKappaTransferMatrix(20, 11, 5000, self.constant5000Burstgen, 5.475, R, kappa, weights, RRange = (5, 6))
        tm.generateMatrix()
        self.assertEqual(tm.RRange[0], 5)
        self.assertEqual(tm.RRange[1], 6)
        self.assertAlmostEqual(tm.getMatrix()[9][5], 1., delta = 0.01)
        tmref = GlobalAVGKappaTransferMatrix(20, 11, 5000, self.constant5000Burstgen, 5.475, [5, 6])
        tmref2 = GlobalAVGKappaTransferMatrix(20, 11, 5000, self.constant5000Burstgen, 5.475, [5, 6])
        self.assertMatrixAlmostEqual(tmref2.getMatrix(), tmref.getMatrix(), delta = 0.05)

#    def testDistanceKappaTransferMatrixvsGlobal(self):
#        R = numpy.linspace(5.05, 5.95, 20)
#        kappa = [ 2. / 3 ] * 20
#        weights = [ 1. ] * 20
#        tm = DistanceKappaTransferMatrix(20, 11, 5000, self.constant5000Burstgen, 5.475, R, kappa, weights, RRange = (5, 6))
#        tm.generateMatrix()
#        self.assertEqual(tm.RRange[0], 5)
#        self.assertEqual(tm.RRange[1], 6)
#        self.assertAlmostEqual(tm.getMatrix()[9][5], 1., delta = 0.01)
#        tmref = GlobalAVGKappaTransferMatrix(20, 11, 5000, self.constant5000Burstgen, 5.475, [5, 6])
#        self.assertMatrixAlmostEqual(tm.getMatrix(), tmref.getMatrix(), delta = 0.05)

    def testDistanceAVGLocalKappaAveraging(self):
        R = numpy.concatenate((numpy.linspace(5.05, 5.95, 20), numpy.linspace(5.05, 5.95, 20), numpy.linspace(5.05, 5.95, 20)))
        kappa = [ 1.  ] * 20 + [ 1. ] * 20 + [ 0. ] * 20
        weights = [ 1. ] * 60
        tm = DistanceAVGKappaTransferMatrix(20, 11, 5, self.constant5000Burstgen, 5.475, R, kappa, weights, RRange = (5, 6))
        self.assertEqual(tm.getKappaAVG()[2], 2. / 3)
        self.assertEqual(tm.RRange[0], 5)
        self.assertEqual(tm.RRange[1], 6)
        self.assertEqual(tm.getMatrix()[9][5], 1.)

    def testPlotting(self):
        tm = GlobalAVGKappaTransferMatrix(20, 11, 5, self.constant5000Burstgen, 5.475, [5, 6])
        odir = tempfile.mkdtemp(suffix = "plottest")
        ofl = os.path.join(odir, "plot.png")
        tm.plotToFile(ofl)
        self.assertEqual(os.path.exists(ofl), True)
        shutil.rmtree(odir)

    def testRandomUncorr(self):
        length = 100000
        startdist = 4
        enddist = 7
        rbins = 10
        ebins = 20
        bursts = 10000
        R0 = 5.475
        bgen = self.constant50Burstgen
        R = numpy.random.random(length) * (enddist - startdist) + startdist
        kappa2 = numpy.array(list(getKappa(genRandomVec(), genRandomVec(), genRandomVec()) ** 2 for _ in range(length)))
        print "Kappa^2 mean is ", kappa2.mean()
        R0mod = modifyR0(R0, kappa2.mean())
        prbs = numpy.ones(length)
        globaltm = GlobalAVGKappaTransferMatrix(rbins, ebins, bursts, bgen, R0mod, (startdist, enddist))
        localtm = DistanceAVGKappaTransferMatrix(rbins, ebins, bursts, bgen, R0, R, kappa2, prbs, RRange = (startdist, enddist))
        self.assertMatrixAlmostEqual(globaltm.getMatrix(), localtm.getMatrix(), delta = 0.10)

#        nonetm = DistanceKappaTransferMatrix(rbins, ebins, bursts, bgen, R0, R, kappa2, prbs, RRange = (startdist, enddist))
#        self.assertMatrixAlmostEqual(globaltm.getMatrix(), nonetm.getMatrix(), delta = 0.05)

#    def testKappaDistanceCorr(self):
#        length = 100000
#        startdist = 4
#        enddist = 7
#        rbins = 10
#        ebins = 11
#        bursts = 5000
#        R0 = 5.475
#        bgen = self.constant500Burstgen
#        R = numpy.random.random(length) * (enddist - startdist) + startdist
#        ktmp = []
#        for i in range(length):
#            ktmp.append(getKappa(genRandomVec(), genRandomVec(), genRandomVec()) ** 2 * startdist / R[i])
#        kappa2 = numpy.array(ktmp)
#        prbs = numpy.ones(length)
#        globaltm = GlobalAVGKappaTransferMatrix(rbins, ebins, bursts, bgen, R0, (startdist, enddist))
#        localtm = DistanceAVGKappaTransferMatrix(rbins, ebins, bursts, bgen, R0, R, kappa2, prbs, RRange = (startdist, enddist))
#        self.assertNotAlmostEqual(matDiff(globaltm.getMatrix(), localtm.getMatrix()), 0.0, delta = 0.1)
#        nonetm = DistanceKappaTransferMatrix(rbins, ebins, bursts, bgen, R0, R, kappa2, prbs, RRange = (startdist, enddist))
#
#        self.assertMatrixAlmostEqual(localtm.getMatrix(), nonetm.getMatrix(), delta = 0.05)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
