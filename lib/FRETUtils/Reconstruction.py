'''
Created on 15.10.2012

@author: mhoefli
'''
from FRETUtils.Bursts import readBurstSizes, genPowerlawTable,\
    getAcceptRejectBurst
import random
import numpy

from TransferMatrices import GlobalAVGKappaTransferMatrix, DistanceAVGKappaTransferMatrix, DistanceKappaTransferMatrix
from TransferMatrixInversion import GaussianRegularizationDistanceReconstruction

transferMatrices=("global","local","none")

def getExpBurstgen(bsdfile):
    print "Using experimental BSD from",bsdfile 
    bsizes=readBurstSizes(bsdfile)
    burstGenerator = lambda : random.choice(bsizes)
    return burstGenerator

def getFitBurstgen(burstsmin,burstsmax,exp=-2.3):
    print "Using fitted BSD from",burstsmin,"to",burstsmax
    bstable=genPowerlawTable(burstsmin,burstsmax,exp)
    burstGenerator = lambda : getAcceptRejectBurst(bstable,burstsmin,burstsmax)
    return burstGenerator

def readEfficiencies(efficiencyfile, efficiencybins):
    effs=numpy.loadtxt(efficiencyfile)
    effhist = numpy.histogram(effs, efficiencybins, (0.,1.), normed=True)
    return effhist

def readRKPrbSamples(rkappafile):
    arr = numpy.loadtxt(rkappafile,comments="#")
    #discarding time information
    return arr[:,1],arr[:,2],arr[:,3]

def constructGlobalTM(options,burstGenerator):
    TM = GlobalAVGKappaTransferMatrix(options.distancebins,options.efficiencybins,options.burstcount,burstGenerator, options.R0, (options.distancestart,options.distanceend))
    return TM.getMatrix()

def constructNonGlobalTM(TMType,options,burstGenerator):
    if options.distancestart and options.distanceend:
        myrange=(options.distancestart,options.distanceend)
    elif not options.distancestart and not options.distanceend:
        myrange=None
    else:
        raise ValueError("Specify distance range by using both, -rs and -re or none of both (autodetect from -r file).")
    RSamples,KappaSamples,SampleWeights = readRKPrbSamples(options.rkappafile)
    TM = TMType(options.distancebins,options.efficiencybins,options.burstcount,burstGenerator, options.R0,RSamples,KappaSamples,SampleWeights,RRange=myrange)
    return TM.getMatrix()

def constructLocalTM(options,burstGenerator):
    return constructNonGlobalTM(DistanceAVGKappaTransferMatrix,options,burstGenerator)
    

def constructNoAVGTM(options,burstGenerator):
    return constructNonGlobalTM(DistanceKappaTransferMatrix,options,burstGenerator)

def getBurstGenerator(options):
    if options.expbfile:
        return getExpBurstgen(options.expbfile)
    else:
        return getFitBurstgen(options.burstMinsize, options.burstMaxsize, options.burstLambda)

def constructTM(options):
    burstgen = getBurstGenerator(options)
    if options.transferMatrix=="global":
        return constructGlobalTM(options,burstgen)
    elif options.transferMatrix=="local":
        return constructLocalTM(options,burstgen)
    elif options.transferMatrix=="none":
        return constructNoAVGTM(options,burstgen)
    else:
        raise ValueError("Invalid transfer matrix type selected") 

def resolveDistances(options, TM, effhist):
    return GaussianRegularizationDistanceReconstruction(options, TM, effhist)

def writeDistances(distances,options):
    numpy.savetxt(options.outdistfile,distances)
