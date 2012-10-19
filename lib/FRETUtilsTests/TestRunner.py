'''
Created on 04.10.2012

@author: martin
'''

import unittest,os

from FRETUtilsTests.EfficiencyTests import testEfficiencyFunctions
from FRETUtilsTests.EnsembleTests import testEnsembleFunctions
from FRETUtilsTests.HelperTests import testHelperFunctions
from FRETUtilsTests.PhotonsTests import testPhotonFunctions
from FRETUtilsTests.RunTests import testFullRun, testFullRunMultiprocessing, testFullRKPrbConversion, testFullRKPrbConversionMultiprocessing

def runAllTests(directory="."):
    print os.path.abspath(os.curdir)
    suite = unittest.TestLoader().discover(directory,pattern='*Tests.py')
    unittest.TextTestRunner(verbosity=2).run(suite)

def runSelectedTests():
    pass

if __name__ == "__main__":
    runAllTests()
else:
    runSelectedTests()