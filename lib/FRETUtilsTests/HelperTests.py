'''
Created on 04.10.2012

@author: martin
'''

import unittest
import FRETUtils.Helpers as Helpers

class testHelperFunctions(unittest.TestCase):
    def testJtoNm(self):
        self.assertAlmostEqual(Helpers.JouleInNm(42.), 0., 4)

    def testRToEff(self):
        self.assertAlmostEqual(Helpers.rToEff(100., 5.), 0., 4)
        self.assertAlmostEqual(Helpers.rToEff(0.1, 5.), 1., 4)
        self.assertAlmostEqual(Helpers.rToEff(5., 5.), 0.5, 4)

    def testEffToR(self):
        self.assertGreater(Helpers.effToR(0.001, 5.), 10.)
        self.assertLess(Helpers.effToR(0.999999, 5.), 1.)
        self.assertAlmostEqual(Helpers.effToR(0.5, 5.), 5., 4)
