'''
Created on 02.10.2012

@author: mhoefli
'''

import unittest
from FRETUtils import Photons


class testPhotonFunctions(unittest.TestCase):
    def setUp(self):
        self.donorphoton = Photons.Photon(True,1000,4000,thermal=False)
        self.acceptorphoton = Photons.Photon(False,1000,4000,thermal=False)
    
    def test_donor_thermal(self):
        self.donorphoton.checkThermal(0.0, 0.0)
        self.assertEqual(self.donorphoton.thermal,True)
    def test_donor_fluorescence(self):        
        self.donorphoton.checkThermal(1.0, 0.0)
        self.assertEqual(self.donorphoton.thermal,False)
    def test_acceptor_thermal(self):
        self.acceptorphoton.checkThermal(0.0, 0.0)
        self.assertEqual(self.acceptorphoton.thermal,True)
    def test_acceptor_fluorescence(self):
        self.acceptorphoton.checkThermal(0.0, 1.0)
        self.assertEqual(self.acceptorphoton.thermal,False)        