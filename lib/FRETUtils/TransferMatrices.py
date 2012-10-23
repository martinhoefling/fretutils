'''
Created on 29.10.2010

@author: mhoefli
'''

from __future__ import division

from FRETUtils.Bursts import getBurstSizes
    
import numpy as np
from FRETUtils.Helpers import rToEff
from FRETUtils.Efficiencies import efficiencyDeltaFix
from matplotlib import pyplot as plt
from matplotlib import cm as cm


def generateBinMid(myrange, ind,RBins, RStart):
    binmid = (ind + 0.5) * myrange / RBins + RStart
    return binmid

def genRandomBurstEffs(reff, size):
    randnrs = np.random.random(size)
    nacc = randnrs < reff
    effval = nacc.sum() / size
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
    
    def plot(self):
        xdist = np.linspace(self.RRange[0]+self.myrange/self.RBins/2,self.RRange[1]-self.myrange/self.RBins/2,self.RBins)
        xeff = np.linspace(0., 1., self.EffBins, endpoint=False) + 1. / self.EffBins / 2
        
        xfifth=self.RBins/5
        yfifth=self.EffBins/5
        maxscale = self.tm[xfifth:-xfifth,yfifth:-yfifth].max()
        plt.imshow(self.tm.T, interpolation='nearest', aspect="auto", origin="lower", cmap=cm.jet, vmin=self.tm.min(), vmax=maxscale)
        plt.xlabel("Distance")
        plt.ylabel("Efficiency")
        nrlabels = 4
        nrx = len(xdist)
        nre = len(xeff)
        stepx = nrx / nrlabels
        stepe = nre / nrlabels
        ndxx = (np.arange(nrx) % stepx == 0).nonzero()
        ndxe = (np.arange(nre) % stepe == 0).nonzero()
        plt.xticks(np.arange(len(xdist))[ndxx], xdist[ndxx])
        plt.yticks((xeff * len(xeff))[ndxe], xeff[ndxe])



class GlobalAVGKappaTransferMatrix(TransferMatrix):
    def __init__(self,RBins,EffBins,BurstCount,burstGenerator, R0, RRange):
        TransferMatrix.__init__(self, RBins, EffBins, BurstCount, burstGenerator, R0, RRange)
       
    def populateMatrixWithBurst(self,rbinindex,reff,burst,weight):      
        effval = genRandomBurstEffs(reff, burst)
        effval = efficiencyDeltaFix((effval,))
        effndx = genEffIndex(effval,self.EffBins)
        self.tm[rbinindex,effndx]+=weight
    
    def populateMatrix(self,rbinindex):
        reff,bursts = self.getBinEfficiencies(rbinindex)
        weight=1./len(reff)
        for eff,burst in zip(reff,bursts):
            self.populateMatrixWithBurst(rbinindex, eff, burst, weight)

    def getR0(self,binnr):
        return self.R0

    def getBinEfficiencies(self, binnr):
        binmid = generateBinMid(self.myrange, binnr, self.RBins, self.RRange[0])
        myR0 = self.getR0(binnr)
        reff = rToEff(binmid, myR0)
        return [reff]*self.BurstCount,getBurstSizes(self.BurstCount,self.burstGenerator)

    def generateMatrix(self):
        for i in range(self.RBins):
            self.populateMatrix(i)
       
        self.matrixGenerated=True

class DistanceAVGKappaTransferMatrix(GlobalAVGKappaTransferMatrix):
    def __init__(self,RBins,EffBins,BurstCount,burstGenerator,R0,RSamples,KappaSamples,SampleWeights,RRange=None):        
        if not RRange:
            GlobalAVGKappaTransferMatrix.__init__(self, RBins, EffBins, BurstCount, burstGenerator, R0, (RSamples.min(),RSamples.max()))
        else:
            GlobalAVGKappaTransferMatrix.__init__(self, RBins, EffBins, BurstCount, burstGenerator, R0, RRange)        
        self.RSamples=RSamples
        self.KappaSamples=KappaSamples
        self.SampleWeights=SampleWeights
        
        self.kappaavg=np.zeros((RBins),np.float64)
        self.kappaBinned=False
    

    def addToBin(self, kappaavgnum, R, K, prb, Rind):
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
            self.addToBin(kappaavgnum, R, K, prb, Rind)
        
        self.binKappaSanityCheck(kappaavgnum)
        self.kappaavg/=kappaavgnum
        self.kappaBinned=True
    
    def binKappaSanityCheck(self,kappaavgnum):
        if not (kappaavgnum>0).all():
            raise ValueError("Not all kappa bins have at least one sample. Adjust number of distance-bins or distance range")
        print "There are at least %d kappa samples in each r-bin."%kappaavgnum.min()
     
    def getR0(self,binnr):
        if not self.kappaBinned:
            self.binKappa()
        return modifyR0(self.R0, self.kappaavg[binnr])
        
    def getKappaAVG(self):
        if not self.kappaBinned:
            self.binKappa()
        return self.kappaavg


class DistanceKappaTransferMatrix(DistanceAVGKappaTransferMatrix):    
    def __init__(self,RBins,EffBins,BurstCount,burstGenerator,R0,RSamples,KappaSamples,SampleWeights,RRange=None):
        DistanceAVGKappaTransferMatrix.__init__(self, RBins, EffBins, BurstCount, burstGenerator, R0, RSamples, KappaSamples, SampleWeights, RRange)   
        self.rkappaBinned = [ [] for _i in range(RBins) ]

    def addToBin(self, kappaavgnum, R, K, prb, Rind):
        super(DistanceKappaTransferMatrix,self).addToBin(kappaavgnum, R, K, prb, Rind)
        self.rkappaBinned[Rind].append((R,K,prb))
        
    def getBinEfficiencies(self, binnr):
        if not self.kappaBinned:
            self.binKappa()
            
        if len(self.rkappaBinned[binnr])==0:
            raise ValueError("Bin #%d is empty. This is a problem. Try setting a smaller range and fewer bins in R-direction."%binnr)
            
        bursts=getBurstSizes(self.BurstCount,self.burstGenerator)
        beffs=[]
        
        rkprb=np.array(self.rkappaBinned[binnr])
        Rarr=rkprb[:,0]
        Kappaarr=rkprb[:,1]
        Prbarr=rkprb[:,2]
        R0_mod=modifyR0(self.R0,Kappaarr)

        cumulprb=Prbarr.cumsum()
        cumulprb/=cumulprb[-1]
        for burst in bursts:
            randnrs=np.random.random(burst)
            ndxchoice=cumulprb.searchsorted(randnrs)
            effs=rToEff(Rarr[ndxchoice],R0_mod[ndxchoice])
            beffs.append(effs.mean())

        return beffs,bursts
                
   
