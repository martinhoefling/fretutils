# -*- coding: utf-8 -*-
'''
Created on 10.06.2010

@author: martin
'''

from Helpers import FunctionIterator
from numpy import random, arange, fromiter
import numpy

class Burst:
    def __init__(self,burstsize):
        self.bsize=burstsize
        self.photons=[]
        self.donorphot=0
        self.donortherm=0
        self.acceptorphot=0
        self.acceptortherm=0
        
    def appendPhoton(self,photon):
        self.photons.append(photon)
        if photon.donor:
            if photon.thermal:
                self.donortherm+=1
            else:
                self.donorphot+=1
        else:
            if photon.thermal:
                self.acceptortherm+=1
            else:
                self.acceptorphot+=1
    
    def checkSizeReached(self,QD,QA,QYcorrected=False):
        if QYcorrected and float(self.acceptorphot)/QA+float(self.donorphot)/QD>self.bsize:
            return True
        elif (not QYcorrected) and self.acceptorphot+self.donorphot==self.bsize:
            return True

    def bsizeCorr(self,QD,QA):
        return float(self.acceptorphot)/QA+float(self.donorphot)/QD
     
    def getEfficiency(self,QD,QA):
        return float(self.acceptorphot)/QA/(float(self.acceptorphot)/QA+float(self.donorphot)/QD)
    

def genPowerlawTable(mina,maxa,powerl):
    """generates a lookup table according to a power law"""
    return arange(mina,maxa+1)**powerl

def getAcceptRejectBurst(table,mina,maxa):
    """generates a random number based on a table with given probabilities"""
    maxval=table.max()
    while True:
        randint=random.random_integers(mina,maxa)
        randfloat=random.uniform(0.0,maxval)
        if randfloat < table[randint-mina]:
            return randint

def readBurstSizes(fname,corrected=True):
    """Reads experimental burst sizes from FRET comma separated file with 4 columns. 1st & 2nd raw donor and acceptor count, 3rd and 4th corrected count"""
    expbs=[]
    with open(fname) as fh:

        for line in fh:
            spl = line.split(",")
            if len(spl)==4:
                if corrected:
                    acceptor=float(spl[2])
                    donor=float(spl[3])   
                else:
                    acceptor=float(spl[0])
                    donor=float(spl[1])
                
                expbs.append(int(acceptor+donor))

    return numpy.array(expbs)

def getBurstSizes(bcount,burstFunc):
    """ Returns bcount bursts. The generator function and its parameters can be set via burstFunc.
    >>> a=getBurstSizes(100)
    >>> len(a)
    100
    >>> random.seed(1)
    >>> getBurstSizes(10)
    array([ 92,  27, 100,  46,  23,  74,  44,  73,  28,  70])
    """
    
    return fromiter(FunctionIterator(burstFunc,count=bcount),numpy.int)


