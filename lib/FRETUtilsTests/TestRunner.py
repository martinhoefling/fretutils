'''
Created on 04.10.2012

@author: martin
'''

import unittest, os


from FRETUtilsTests.EfficiencyTests import testEfficiencyFunctions  # @UnusedImport IGNORE:W0611
from FRETUtilsTests.EnsembleTests import testEnsembleFunctions  # @UnusedImport IGNORE:W0611
from FRETUtilsTests.HelperTests import testHelperFunctions  # @UnusedImport IGNORE:W0611
from FRETUtilsTests.PhotonsTests import testPhotonFunctions  # @UnusedImport IGNORE:W0611
from FRETUtilsTests.ConfigTests import testSecureConfigParser  # @UnusedImport IGNORE:W0611
from FRETUtilsTests.RunTests import testFullRun, testFullRunMultiprocessing, testFullRKPrbConversion, testFullRKPrbConversionMultiprocessing  # @UnusedImport IGNORE:W0611
from FRETUtilsTests.TransferMatrixTests import TransferMatrixTests  # @UnusedImport IGNORE:W0611

def runAllTests(directory = "."):
    print os.path.abspath(os.curdir)
    suite = unittest.TestLoader().discover(directory, pattern = '*Tests.py')
    unittest.TextTestRunner(verbosity = 2).run(suite)

if __name__ == "__main__":
    runAllTests()
