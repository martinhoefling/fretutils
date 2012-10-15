'''
Created on 10.06.2010

@author: martin
'''
import Constants

def JouleInNm(joule):
    return Constants.h_planck*Constants.c_light/joule

class FunctionIterator():
    def __init__(self,bFunc,count=-1):
        self.bFunc=bFunc
        self.count=count
        self.ccount=0
    
    def __iter__(self):
        return self
    
    def next(self):
        if self.count==self.ccount:
            raise StopIteration
        self.ccount += 1
        return self.bFunc()

def effToR(eff,R0=5.4):
    return R0*pow((1./eff-1),1./6)

def rToEff(R,R0=5.4):
    return 1./(1.+pow(R/R0,6.))