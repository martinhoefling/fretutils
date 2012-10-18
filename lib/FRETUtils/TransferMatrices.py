'''
Created on 29.10.2010

@author: mhoefli
'''

from __future__ import division

from FRETUtils.Bursts import getBurstSizes
    
import numpy as np
from FRETUtils.Helpers import rToEff
from FRETUtils.Efficiencies import efficiencyDeltaFix



def generateBinMid(self, myrange, i,RBins, RStart):
    binmid = (i + 0.5) * myrange / RBins + RStart
    return binmid

def genRandomBurstEff(reff, size):
    randnrs = np.random.random(size)
    nacc = randnrs < reff.sum()
    effval = nacc / size
    return effval

def genEffIndex(effval,EffBins):
    effndx = int(EffBins * effval)
    if effndx >= EffBins:
        effndx -= 1
    return effndx

def getRange(rrange):
    return rrange[1]-rrange[0]

def modifyR0(R0,newkappa):
    return R0*pow((newkappa/(2./3)),1./6)


class TransferMatrix(object):
    def __init__(self,RBins,EffBins,BurstCount,burstGenerator, R0, RRange):
        self.tm=np.zeros((RBins,EffBins),np.float64)
        self.RBins=RBins
        self.EffBins=EffBins
        self.BurstCount=BurstCount
        self.burstGenerator =  burstGenerator
        self.RRange=RRange
        self.myrange=getRange(self.RRange)
        self.R0=R0

        self.matrixGenerated=False
            
    def generateMatrix(self):
        raise NotImplementedError("Please implement this function!")
    
    def getMatrix(self):
        if not self.matrixGenerated:
            self.generateMatrix()
        return self.tm


class GlobalAVGKappaTransferMatrix(TransferMatrix):
       
    def populateMatrixWithBursts(self,rbinindex,reff,bursts,weight):      
        for size in bursts:
            effval = genRandomBurstEff(reff, size)
            effval = efficiencyDeltaFix(effval)
            effndx = genEffIndex(effval,self.EffBins)
            self.tm[rbinindex,effndx]+=weight
    
    def populateMatrix(self,rbinindex):
        reff,bursts = self.getBinEfficiencies(i)
        weight=1./len(reff)
        for eff in zip(reff,bursts):
            self.populateMatrixWithBursts(rbinindex, reff, bursts, weight)

    def getR0(self,binnr):
        return self.R0

    def getBinEfficiencies(self, binnr):
        binmid = generateBinMid(self.myrange, binnr, self.RBins, self.RRange[0])
        myR0 = self.getR0(binnr)
        reff = rToEff(binmid, myR0)
        return [reff]*self.BurstCount,getBurstSizes(self.BurstCount,self.burstGenerator)

    def generateMatrix(self):
        myrange=getRange(self.RRange)
        for i in range(self.RBins):
            self.populateMatrix(i)
       
        self.matrixGenerated=True

class DistanceAVGKappaTransferMatrix(GlobalAVGKappaTransferMatrix):
    def __init__(self,RBins,EffBins,BurstCount,burstGenerator,R0,RSamples,KappaSamples,SampleWeights,RRange=None):        
        if not RRange:
            GlobalAVGKappaTransferMatrix.__init__(self, RBins, EffBins, BurstCount, burstGenerator, (RSamples.min(),RSamples.max()))
        else:
            GlobalAVGKappaTransferMatrix.__init__(self, RBins, EffBins, BurstCount, burstGenerator, R0, RRange)        
        self.RSamples=RSamples
        self.KSamples=KappaSamples
        self.SampleWeights=SampleWeights
        
        self.kappaavg=np.zeros((RBins),np.float64)
        self.kappaBinned=False
    

    def addToBin(self, kappaavgnum, K, prb, Rind):
        self.kappaavg[Rind] += K * prb
        kappaavgnum[Rind] += prb


    def genRIndex(self, R):
        Rind = int((R - self.RRange[0]) / (self.myrange / self.RBins))
        if Rind == self.RBins:
            Rind -= 1
        return Rind

    def binKappa(self):
        kappaavgnum=np.zeros((self.RBins),np.float64)
        for R,K,prb in zip(self.RSamples,self.KappaSamples,self.SampleWeights):
            Rind = self.genRIndex(R)
            if Rind < 0 or Rind > self.RBins:
                continue 
            self.addToBin(kappaavgnum, K, prb, Rind)
        
        self.kappaavg/=kappaavgnum
        self.kappaBinned=True
     
    def getR0(self,binnr):
        return modifyR0(self.R0, self.kappaavg[binnr])
        
    def getKappaAVG(self):
        if not self.kappaBinned:
            self.binKappa()
        return self.kappaavg


class DistanceKappaTransferMatrix(istanceAVGKappaTransferMatrix):    
    def __init__(self,RBins,EffBins,BurstCount,burstGenerator,R0,RSamples,KappaSamples,SampleWeights,RRange=None):
         DistanceKappaTransferMatrix.__init__(self, RBins, EffBins, BurstCount, burstGenerator, R0, RSamples, KappaSamples, SampleWeights, RRange)   
         self.rkappaBinned = [ [] for i in range(RBins) ]

    def addToBin(self, kappaavgnum, K, prb, Rind):
        super(DistanceKappaTransferMatrix,self).addToBin(kappaavgnum, K, prb, Rind)
        self.rkappaBinned[Rind].append((R,K,prb))
        
    def getBinEfficiencies(self, binnr):
        bursts=getBurstSizes(self.BurstCount,self.burstGenerator)
        beffs=[]
        rkprb=np.arr(self.rkappaBinned[binnr])
        Rarr=rkprb[:0]
        Kappaarr=rkprb[:1]
        Prbarr=rkprb[:2]
        R0_mod=modifyR0(self.R0,Kappaarr)

        cumulprb=Prbarr.cumsum()
        cumulprb/=cumulprb[-1]
        for burst in bursts:
            randnrs=numpy.random.random(len(burst))
            ndxchoice=cumulprb.searchsorted(randnrs)
            effs=rToEff(Rarr[ndxchoice],R0_mod[ndxchoice])
            beffs.append(effs.mean())

        return beffs,bursts
                
   
