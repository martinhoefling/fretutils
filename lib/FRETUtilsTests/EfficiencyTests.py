'''
Created on 09.10.2012

@author: mhoefli
'''

import unittest
import numpy
from FRETUtils.Efficiencies import efficiencyDeltaFix

class testEfficiencyFunctions(unittest.TestCase):
    def test_deltaFix(self):
        arr=numpy.random.rand(1000)
        aavg=arr.mean()
        self.assertAlmostEqual(aavg, efficiencyDeltaFix(arr).mean(), delta=0.01)
        