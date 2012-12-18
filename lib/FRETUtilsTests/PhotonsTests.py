'''
Created on 02.10.2012

@author: mhoefli
'''

import unittest
from FRETUtils import Photons
from FRETUtils.Config import SecureConfigParser


class DummyConfigParser(SecureConfigParser):
    def __init__(self):
        SecureConfigParser.__init__(self)
        self.setdefault("System", "photongenerator", "holla", str, lambda x: x in ["cython", "python", "cextension", "holla"] , "must be one of cython, python, holla or cextension")

class testPhotonFunctions(unittest.TestCase):
    def setUp(self):
        self.donorphoton = Photons.Photon(True, 1000, 4000, thermal = False)
        self.acceptorphoton = Photons.Photon(False, 1000, 4000, thermal = False)

    def test_donor_thermal(self):
        self.donorphoton.checkThermal(0.0, 0.0)
        self.assertEqual(self.donorphoton.thermal, True)
    def test_donor_fluorescence(self):
        self.donorphoton.checkThermal(1.0, 0.0)
        self.assertEqual(self.donorphoton.thermal, False)
    def test_acceptor_thermal(self):
        self.acceptorphoton.checkThermal(0.0, 0.0)
        self.assertEqual(self.acceptorphoton.thermal, True)
    def test_acceptor_fluorescence(self):
        self.acceptorphoton.checkThermal(0.0, 1.0)
        self.assertEqual(self.acceptorphoton.thermal, False)
    def test_InvalidPhotonIndex(self):
        conf = "dummy"
        self.assertRaises(IndexError, Photons.getPhoton, [0, 1, 2, 3, 4, 5], conf, -1)
    def test_InvalidPhotonGeneratorSet(self):
        cfg = DummyConfigParser()
        self.assertRaises(ValueError, Photons.setPhotonGenerator, cfg)
