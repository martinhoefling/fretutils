'''
Created on 10.06.2010

@author: martin
'''
import Constants
import numpy
from math import sqrt

def JouleInNm(joule):
    return Constants.h_planck * Constants.c_light / joule

class FunctionIterator():
    def __init__(self, bFunc, count = -1):
        self.bFunc = bFunc
        self.count = count
        self.ccount = 0

    def __iter__(self):
        return self

    def next(self):
        if self.count == self.ccount:
            raise StopIteration
        self.ccount += 1
        return self.bFunc()

def effToR(eff, R0 = 5.4):
    return R0 * pow((1. / eff - 1), 1. / 6)

def rToEff(R, R0 = 5.4):
    return 1. / (1. + pow(R / R0, 6.))

def getKappa(donor, acceptor, distance):
    return (donor * acceptor).sum() - 3 * (donor * distance).sum() * (acceptor * distance).sum()

def genRandomVec():
    gauss = numpy.random.normal(size = 3)
    return gauss / sqrt((gauss ** 2).sum())
