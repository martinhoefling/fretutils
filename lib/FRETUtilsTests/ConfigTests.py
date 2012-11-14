'''
Created on 29.10.2012

@author: mhoefli
'''

import unittest,tempfile

from FRETUtils import Config #@UnresolvedImport

class testSecureConfigParser(unittest.TestCase):
    def testAddWithoutDefault(self):
        config = Config.SecureConfigParser()
        self.assertRaises(ValueError,config.set,"Invalid","Invalid","a")
    
    def testAddTwice(self):
        config = Config.SecureConfigParser()
        config.setdefault("Invalid","Invalid","a",str, lambda x: True,"must be string")
        self.assertRaises(ValueError,config.setdefault,"Invalid","Invalid","a",str, lambda x: True,"must be string")
        
    def testSetInvalidType(self):
        config = Config.SecureConfigParser()
        config.setdefault("Invalid","Invalid",1,int,lambda x:x>0,"must be larger zero")
        self.assertRaises(ValueError,config.set,"Invalid","Invalid","a")
        
    def testSetInvalidValue(self):
        config = Config.SecureConfigParser()
        config.setdefault("Invalid","Invalid",1,int,lambda x:x>0,"must be larger zero")
        self.assertRaises(ValueError,config.set,"Invalid","Invalid",-1)
        
    def testReadWriteReconstruction(self):
        config1 = Config.ReconstructionConfigParser()
        config2 = Config.ReconstructionConfigParser()
        config3 = Config.ReconstructionConfigParser()
        
        tf = tempfile.mktemp()
        with open(tf,"w") as tfh:
            config1.write(tfh)
        with open(tf) as tfh2:
            config2.readfp(tfh2)
            
        config3.read(tf)
        
    