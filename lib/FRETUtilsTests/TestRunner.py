'''
Created on 04.10.2012

@author: martin
'''

import unittest, os

from FRETUtilsTests.EfficiencyTests import testEfficiencyFunctions  # @UnusedImport
from FRETUtilsTests.EnsembleTests import testEnsembleFunctions  # @UnusedImport
from FRETUtilsTests.HelperTests import testHelperFunctions  # @UnusedImport
from FRETUtilsTests.PhotonsTests import testPhotonFunctions  # @UnusedImport
from FRETUtilsTests.ConfigTests import testSecureConfigParser  # @UnusedImport
from FRETUtilsTests.RunTests import testFullRun, testFullRunMultiprocessing, testFullRKPrbConversion, testFullRKPrbConversionMultiprocessing  # @UnusedImport

def runAllTests(directory = "."):
    print os.path.abspath(os.curdir)
    suite = unittest.TestLoader().discover(directory, pattern = '*Tests.py')
    unittest.TextTestRunner(verbosity = 2).run(suite)

if __name__ == "__main__":
    runAllTests()
