'''
Created on 04.10.2012

@author: mhoefli
'''

import unittest
import shutil
import tempfile
import os
from RunTestsHelpers import FakeOptions, FakeOptionsRK, createConstantDummyRKTraj, createProbabilityClassFile, writeConfigFiles, writeRKConfigFiles, writeBurstFile, writeInvalidBurstFile, createStepFunctionDummyRKTrajWithL, writeInvalidDummyFile
from FRETUtils.Run import runMCFRET, runTrajPrbAdd
import numpy

class testFullRun(unittest.TestCase):

    def createTestTrajectories(self):
        # high, no and average FRET
        createConstantDummyRKTraj("high.npz", 10000, 10, 1, 4, fformat = "numpy")
        createConstantDummyRKTraj("low.npz", 10000, 10, 10, 0, fformat = "numpy")
        createConstantDummyRKTraj("aver.npz", 10000, 10, 4, 2. / 3, fformat = "numpy")
        createStepFunctionDummyRKTrajWithL("step.npz", 10000, 10, 0.5, 4., 2. / 3, fformat = "numpy")

        createConstantDummyRKTraj("high.dat", 10000, 10, 1, 4, fformat = "plain")
        createConstantDummyRKTraj("low.dat", 10000, 10, 10, 0, fformat = "plain")
        createConstantDummyRKTraj("aver.dat", 10000, 10, 4, 2. / 3, fformat = "plain")
        createStepFunctionDummyRKTrajWithL("step.dat", 10000, 10, 0.5, 4., 2. / 3, fformat = "plain")

    def createProbabilityClassFiles(self):
        createProbabilityClassFile("high.prb", ("high", "low", "aver", "step"), (1.0, 0.0, 0.0, 0.0))
        createProbabilityClassFile("low.prb", ("high", "low", "aver", "step"), (0.0, 1.0, 0.0, 0.0))
        createProbabilityClassFile("aver.prb", ("high", "low", "aver", "step"), (0.0, 0.0, 1.0, 0.0))
        createProbabilityClassFile("mixed.prb", ("high", "low", "aver", "step"), (0.25, 0.25, 0.5, 0.0))
        createProbabilityClassFile("step.prb", ("high", "low", "aver", "step"), (0.0, 0.0, 0.0, 1.0))
        createProbabilityClassFile("invalid.prb", ("high",), (1.0,))

    def setUp(self):
        self.workdir = tempfile.mkdtemp()
        self.prevdir = os.curdir
        os.chdir(self.workdir)
        self.createTestTrajectories()
        self.createProbabilityClassFiles()
        writeConfigFiles()
        writeBurstFile()
        writeInvalidBurstFile()
        writeInvalidDummyFile()

    def tearDown(self):
        os.chdir(self.prevdir)
        shutil.rmtree(self.workdir)

    def test_output(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.binaryofile = "binary.pkl"
        options.burstcompofile = "burstcomp.txt"
        options.burstsizeofile = "burstsizes.txt"
        options.decaytimeofile = "decaytimes.txt"
        options.endtimeofile = "endtimes.txt"
        options.pbfile = "high.prb"
        options.configfilename = "standard.conf"
        runMCFRET(options)
        self.assertAlmostEqual(numpy.loadtxt("effs.txt").mean(), 1., delta = 0.05)

    def test_invalidDir(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "high.prb"
        options.configfilename = "standard.conf"
        os.mkdir(os.path.join(self.workdir, "empty"))
        options.trajdirectory = "empty"
        self.assertRaises(ValueError, runMCFRET, options)

    def test_profiling(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "high.prb"
        options.configfilename = "standard.conf"
        options.prffile = "profile.log"
        runMCFRET(options)
        self.assertTrue(os.path.exists("profile.log"))

    def test_highEff(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "high.prb"
        options.configfilename = "standard.conf"
        runMCFRET(options)
        self.assertAlmostEqual(numpy.loadtxt("effs.txt").mean(), 1., delta = 0.05)

    def test_lowEff(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "low.prb"
        options.configfilename = "standard.conf"
        runMCFRET(options)
        self.assertAlmostEqual(numpy.loadtxt("effs.txt").mean(), 0., delta = 0.05)

    def test_averEff(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "aver.prb"
        options.configfilename = "standard.conf"
        runMCFRET(options)
        self.assertAlmostEqual(numpy.loadtxt("effs.txt").mean(), 0.5, delta = 0.05)

    def test_samespeciesEff(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "aver.prb"
        options.configfilename = "same-species.conf"
        runMCFRET(options)
        self.assertAlmostEqual(numpy.loadtxt("effs.txt").mean(), 0.5, delta = 0.05)

    def test_invalidBurstGenMethod(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "aver.prb"
        options.configfilename = "invalidburstacc.conf"
        self.assertRaises(ValueError, runMCFRET, options)

    def test_alltrajEff(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "aver.prb"
        options.configfilename = "all.conf"
        runMCFRET(options)
        self.assertAlmostEqual(numpy.loadtxt("effs.txt").mean(), 0.5, delta = 0.05)

    def test_thermalPhot(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "aver.prb"
        options.configfilename = "thermal.conf"
        runMCFRET(options)
        self.assertAlmostEqual(numpy.loadtxt("effs.txt").mean(), 0.5, delta = 0.05)

    def test_correctedPhot(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "aver.prb"
        options.configfilename = "corrected.conf"
        runMCFRET(options)
        self.assertAlmostEqual(numpy.loadtxt("effs.txt").mean(), 0.5, delta = 0.05)

    def test_invalidBSDGen(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "aver.prb"
        options.configfilename = "invalidbsd.conf"
        self.assertRaises(ValueError, runMCFRET, options)

    def test_invalidBSCutOff(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "aver.prb"
        options.configfilename = "invalidbscutoff.conf"
        self.assertRaises(ValueError, runMCFRET, options)

    def test_expBurst(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "aver.prb"
        options.configfilename = "expbst.conf"
        options.expbfile = "expbursts.bst"
        runMCFRET(options)
        self.assertAlmostEqual(numpy.loadtxt("effs.txt").mean(), 0.5, delta = 0.05)

    def test_expBurstCorrected(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "aver.prb"
        options.configfilename = "expbstcorr.conf"
        options.expbfile = "expbursts.bst"
        runMCFRET(options)
        self.assertAlmostEqual(numpy.loadtxt("effs.txt").mean(), 0.5, delta = 0.05)

    def test_cannotAssignClass(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "invalid.prb"
        options.configfilename = "expbstcorr.conf"
        options.expbfile = "expbursts.bst"
        self.assertRaises(ValueError, runMCFRET, options)

    def test_invalidMinMaxRange(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "aver.prb"
        options.configfilename = "wrongminmaxstart.conf"
        options.expbfile = "expbursts.bst"
        self.assertRaises(IndexError, runMCFRET, options)

    def test_invalidMinMaxStartTraj(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "aver.prb"
        options.configfilename = "wrongminmaxstart2.conf"
        options.expbfile = "expbursts.bst"
        self.assertRaises(ValueError, runMCFRET, options)

    def test_rejectPhotonTest(self):
        for i in ("trajectory", "same-species", "all"):
            options = FakeOptions()
            options.efficiencyofile = "effs.txt"
            options.pbfile = "aver.prb"
            options.configfilename = "rejecttest-%s.conf" % i
            options.expbfile = "expbursts.bst"
            self.assertRaises(ValueError, runMCFRET, options)

    def test_stepFunctionBelowReject(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "step.prb"
        options.configfilename = "stepreject.conf"
        options.expbfile = "expbursts.bst"
        runMCFRET(options)
        self.assertAlmostEqual(numpy.loadtxt("effs.txt").mean(), 0.5, delta = 0.05)

    def test_trajFormat(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.trajformat = "dat"
        options.pbfile = "step.prb"
        options.configfilename = "stepreject.conf"
        options.expbfile = "expbursts.bst"
        runMCFRET(options)
        self.assertAlmostEqual(numpy.loadtxt("effs.txt").mean(), 0.5, delta = 0.05)

    def test_invalidTrajFormat(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.trajformat = "inv"
        options.pbfile = "step.prb"
        options.configfilename = "stepreject.conf"
        options.expbfile = "expbursts.bst"
        self.assertRaises(ValueError, runMCFRET, options)


    def test_invalidPhotongenerator(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "step.prb"
        options.configfilename = "invalid_photongenerator.conf"
        options.expbfile = "expbursts.bst"
        self.assertRaises(ValueError, runMCFRET, options)


    def test_cythonPhotongenerator(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "step.prb"
        options.configfilename = "cython.conf"
        options.expbfile = "expbursts.bst"
        try:
#            from FRETUtils.PhotonGenerator import tryGetCythonPhoton #@UnusedImport
            runMCFRET(options)
            self.assertAlmostEqual(numpy.loadtxt("effs.txt").mean(), 0.5, delta = 0.05)
        except ImportError:
            pass
#            from FRETUtils.PhotonGenerator import tryGetCythonPhoton
#            self.assertRaises(ImportError,runMCFRET,options)

    def test_cextensionPhotongenerator(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "step.prb"
        options.configfilename = "cextension.conf"
        options.expbfile = "expbursts.bst"
        try:
#            from FRETUtils.PhotonGenerator import tryGetCPhoton #@UnusedImport
            runMCFRET(options)
            self.assertAlmostEqual(numpy.loadtxt("effs.txt").mean(), 0.5, delta = 0.05)
        except ImportError:
            pass
#            from FRETUtils.PhotonGenerator import tryGetCPhoton
#            self.assertRaises(ImportError,runMCFRET,options)

class testFullRunMultiprocessing(unittest.TestCase):

    def createTestTrajectories(self):
        # high, no and average FRET
        createConstantDummyRKTraj("high.npz", 10000, 10, 1, 4, fformat = "numpy")
        createConstantDummyRKTraj("low.npz", 10000, 10, 10, 0, fformat = "numpy")
        createConstantDummyRKTraj("aver.npz", 10000, 10, 4, 2. / 3, fformat = "numpy")
        createStepFunctionDummyRKTrajWithL("step.npz", 10000, 10, 0.5, 4., 2. / 3, fformat = "numpy")

        createConstantDummyRKTraj("high.dat", 10000, 10, 1, 4, fformat = "plain")
        createConstantDummyRKTraj("low.dat", 10000, 10, 10, 0, fformat = "plain")
        createConstantDummyRKTraj("aver.dat", 10000, 10, 4, 2. / 3, fformat = "plain")
        createStepFunctionDummyRKTrajWithL("step.dat", 10000, 10, 0.5, 4., 2. / 3, fformat = "plain")

    def createProbabilityClassFiles(self):
        createProbabilityClassFile("high.prb", ("high", "low", "aver", "step"), (1.0, 0.0, 0.0, 0.0))
        createProbabilityClassFile("low.prb", ("high", "low", "aver", "step"), (0.0, 1.0, 0.0, 0.0))
        createProbabilityClassFile("aver.prb", ("high", "low", "aver", "step"), (0.0, 0.0, 1.0, 0.0))
        createProbabilityClassFile("mixed.prb", ("high", "low", "aver", "step"), (0.25, 0.25, 0.5, 0.0))
        createProbabilityClassFile("step.prb", ("high", "low", "aver", "step"), (0.0, 0.0, 0.0, 1.0))
        createProbabilityClassFile("invalid.prb", ("high",), (1.0,))

    def setUp(self):
        self.workdir = tempfile.mkdtemp()
        self.prevdir = os.curdir
        os.chdir(self.workdir)
        self.createTestTrajectories()
        self.createProbabilityClassFiles()
        writeConfigFiles()
        writeBurstFile()
        writeInvalidBurstFile()
        writeInvalidDummyFile()

    def tearDown(self):
        os.chdir(self.prevdir)
        shutil.rmtree(self.workdir)

    def test_dualCPU(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "aver.prb"
        options.configfilename = "dual.conf"
        runMCFRET(options)
        self.assertAlmostEqual(numpy.loadtxt("effs.txt").mean(), 0.5, delta = 0.05)

    def test_multiCPU(self):
        options = FakeOptions()
        options.efficiencyofile = "effs.txt"
        options.pbfile = "aver.prb"
        options.configfilename = "multi.conf"
        runMCFRET(options)
        self.assertAlmostEqual(numpy.loadtxt("effs.txt").mean(), 0.5, delta = 0.05)


class testFullRKPrbConversion(unittest.TestCase):
    def createTestTrajectories(self):
        # high, no and average FRET
        createConstantDummyRKTraj("high.npz", 1000, 10, 1, 4, fformat = "numpy")
        createConstantDummyRKTraj("low.npz", 1000, 10, 10, 0, fformat = "numpy")
        createConstantDummyRKTraj("aver.npz", 1000, 10, 4, 2. / 3, fformat = "numpy")
        createStepFunctionDummyRKTrajWithL("step.npz", 1000, 10, 0.5, 4., 2. / 3, fformat = "numpy")

        createConstantDummyRKTraj("high.dat", 1000, 10, 1, 4, fformat = "plain")
        createConstantDummyRKTraj("low.dat", 1000, 10, 10, 0, fformat = "plain")
        createConstantDummyRKTraj("aver.dat", 1000, 10, 4, 2. / 3, fformat = "plain")
        createStepFunctionDummyRKTrajWithL("step.dat", 1000, 10, 0.5, 4., 2. / 3, fformat = "plain")

    def createProbabilityClassFiles(self):
        createProbabilityClassFile("high.prb", ("high", "low", "aver", "step"), (1.0, 0.0, 0.0, 0.0))
        createProbabilityClassFile("low.prb", ("high", "low", "aver", "step"), (0.0, 1.0, 0.0, 0.0))
        createProbabilityClassFile("aver.prb", ("high", "low", "aver", "step"), (0.0, 0.0, 1.0, 0.0))
        createProbabilityClassFile("mixed.prb", ("high", "low", "aver", "step"), (0.25, 0.25, 0.5, 0.0))
        createProbabilityClassFile("step.prb", ("high", "low", "aver", "step"), (0.0, 0.0, 0.0, 1.0))
        createProbabilityClassFile("invalid.prb", ("high",), (1.0,))

    def setUp(self):
        self.workdir = tempfile.mkdtemp()
        self.prevdir = os.curdir
        os.chdir(self.workdir)
        self.createTestTrajectories()
        self.createProbabilityClassFiles()
        writeRKConfigFiles()
        writeBurstFile()
        writeInvalidBurstFile()
        writeInvalidDummyFile()

    def tearDown(self):
        os.chdir(self.prevdir)
        shutil.rmtree(self.workdir)

    def test_trajClassAssignOnly(self):
        options = FakeOptionsRK()
        options.pbfile = "mixed.prb"
        runTrajPrbAdd(options)
        outdata = numpy.loadtxt(options.outtrajfile, comments = "#")
        self.assertAlmostEqual(outdata[0:1000, 3].sum(), 0.5, delta = 0.1)
        self.assertAlmostEqual(outdata[1000:2000, 3].sum(), 0.25, delta = 0.1)
        self.assertAlmostEqual(outdata[2000:3000, 3].sum(), 0.25, delta = 0.1)
        self.assertAlmostEqual(outdata[3000:4000, 3].sum(), 0.0, delta = 0.1)

    def test_trajPhotonFlooding(self):
        options = FakeOptionsRK()
        options.pbfile = "mixed.prb"
        options.configfilename = "standardRK.conf"
        runTrajPrbAdd(options)
        outdata = numpy.loadtxt(options.outtrajfile, comments = "#")
        self.assertAlmostEqual(outdata[:, 3].sum(), 8100. / 1000, delta = 100. / 1000)

    def test_trajClippingPhotonFlooding(self):
        options = FakeOptionsRK()
        options.pbfile = "mixed.prb"
        options.configfilename = "standardRKclip.conf"
        runTrajPrbAdd(options)
        outdata = numpy.loadtxt(options.outtrajfile, comments = "#")
        self.assertEqual(len(outdata[:, 3]), 3760)

    def test_trajPhotonFloodingRejection(self):
        options = FakeOptionsRK()
        options.pbfile = "mixed.prb"
        options.configfilename = "rejectRK.conf"
        runTrajPrbAdd(options)
        outdata = numpy.loadtxt(options.outtrajfile, comments = "#")
        self.assertAlmostEqual(outdata[:, 3].sum(), 0.0, delta = 100.0 / 1000)


class testFullRKPrbConversionMultiprocessing(unittest.TestCase):
    def createTestTrajectories(self):
        # high, no and average FRET
        createConstantDummyRKTraj("high.npz", 1000, 10, 1, 4, fformat = "numpy")
        createConstantDummyRKTraj("low.npz", 1000, 10, 10, 0, fformat = "numpy")
        createConstantDummyRKTraj("aver.npz", 1000, 10, 4, 2. / 3, fformat = "numpy")
        createStepFunctionDummyRKTrajWithL("step.npz", 1000, 10, 0.5, 4., 2. / 3, fformat = "numpy")

        createConstantDummyRKTraj("high.dat", 1000, 10, 1, 4, fformat = "plain")
        createConstantDummyRKTraj("low.dat", 1000, 10, 10, 0, fformat = "plain")
        createConstantDummyRKTraj("aver.dat", 1000, 10, 4, 2. / 3, fformat = "plain")
        createStepFunctionDummyRKTrajWithL("step.dat", 1000, 10, 0.5, 4., 2. / 3, fformat = "plain")

    def createProbabilityClassFiles(self):
        createProbabilityClassFile("high.prb", ("high", "low", "aver", "step"), (1.0, 0.0, 0.0, 0.0))
        createProbabilityClassFile("low.prb", ("high", "low", "aver", "step"), (0.0, 1.0, 0.0, 0.0))
        createProbabilityClassFile("aver.prb", ("high", "low", "aver", "step"), (0.0, 0.0, 1.0, 0.0))
        createProbabilityClassFile("mixed.prb", ("high", "low", "aver", "step"), (0.25, 0.25, 0.5, 0.0))
        createProbabilityClassFile("step.prb", ("high", "low", "aver", "step"), (0.0, 0.0, 0.0, 1.0))
        createProbabilityClassFile("invalid.prb", ("high",), (1.0,))

    def setUp(self):
        self.workdir = tempfile.mkdtemp()
        self.prevdir = os.curdir
        os.chdir(self.workdir)
        self.createTestTrajectories()
        self.createProbabilityClassFiles()
        writeRKConfigFiles()
        writeBurstFile()
        writeInvalidBurstFile()
        writeInvalidDummyFile()

    def test_trajPhotonFlooding(self):
        options = FakeOptionsRK()
        options.pbfile = "mixed.prb"
        options.configfilename = "standardRKmulti.conf"
        runTrajPrbAdd(options)
        outdata = numpy.loadtxt(options.outtrajfile, comments = "#")
        self.assertAlmostEqual(outdata[:, 3].sum(), 16200. / 1000, delta = 100.0 / 1000)



