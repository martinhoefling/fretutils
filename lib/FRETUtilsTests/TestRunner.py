'''
Created on 04.10.2012

@author: martin
'''

import unittest

def runAllTests(directory="."):
    suite = unittest.TestLoader().discover(directory,pattern='*Tests.py')
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    runAllTests()
