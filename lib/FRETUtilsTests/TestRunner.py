'''
Created on 04.10.2012

@author: martin
'''

import unittest,os

def runAllTests(directory="."):
    print os.path.abspath(os.curdir)
    suite = unittest.TestLoader().discover(directory,pattern='*Tests.py')
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    runAllTests()
else:
    runAllTests()