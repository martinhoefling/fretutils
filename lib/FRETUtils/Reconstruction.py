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
    effhist,_bins = numpy.histogram(effs, bins=efficiencybins, range=(0.,1.), normed=True)
    return effhist

def readRKPrbSamples(rkappafile):
    print "Reading rkappa file %s."%rkappafile
    arr = numpy.loadtxt(rkappafile,comments="#")
    print "Done!"
    #discarding time information
    return arr[:,1],arr[:,2],arr[:,3]

def constructGlobalTM(options,config,burstGenerator):
    myrange = getRange(config)
    if options.rkappafile:
        RSamples,KappaSamples,SampleWeights = readRKPrbSamples(options.rkappafile)    
        if not myrange:
            config.set("Transfer Matrix","from distance",RSamples.min())
            config.set("Transfer Matrix","to distance",RSamples.max())
            myrange = getRange(config)


    if not myrange:
        raise ValueError("Range must be specified in config file for TM type global. Alternatively read min and max values from rkappaprb file.")
    
    dbins = config.get("Transfer Matrix","distance bins")
    ebins = config.get("Transfer Matrix","efficiency bins")
    bcount= config.get("Transfer Matrix","bursts per bin")
    R0 = config.get("Transfer Matrix","R0")
    TM = GlobalAVGKappaTransferMatrix(dbins,ebins,bcount,burstGenerator, R0, myrange)
    return TM

def getRange(config):
    if config.get("Transfer Matrix", "from distance") > 0 and config.get("Transfer Matrix", "to distance") > 0:
        myrange = (config.get("Transfer Matrix", "from distance"), config.get("Transfer Matrix", "to distance"))
    elif config.get("Transfer Matrix", "from distance") < 0 and config.get("Transfer Matrix", "to distance"):
        myrange = None
    else:
        raise ValueError("Specify distance range by two positive values, or two negative for autodetect.")
    return myrange

def constructNonGlobalTM(TMType,options,config,burstGenerator):
    myrange = getRange(config)

    RSamples,KappaSamples,SampleWeights = readRKPrbSamples(options.rkappafile)
    
    if not myrange:
        config.set("Transfer Matrix","from distance",RSamples.min())
        config.set("Transfer Matrix","to distance",RSamples.max())
    
    dbins = config.get("Transfer Matrix","distance bins")
    ebins = config.get("Transfer Matrix","efficiency bins")
    bcount= config.get("Transfer Matrix","bursts per bin")
    R0 = config.get("Transfer Matrix","R0")    
    TM = TMType(dbins,ebins,bcount,burstGenerator, R0,RSamples,KappaSamples,SampleWeights,RRange=myrange)
    return TM

def constructLocalTM(options,config,burstGenerator):
    return constructNonGlobalTM(DistanceAVGKappaTransferMatrix,options,config,burstGenerator)
    

def constructNoAVGTM(options,config,burstGenerator):
    return constructNonGlobalTM(DistanceKappaTransferMatrix,options,config,burstGenerator)

def getBurstGenerator(config):
    if config.get("Burst Size Distribution","method")=="analytical":
        llimit=config.get("Burst Size Distribution","llimit")
        ulimit=config.get("Burst Size Distribution","ulimit")
        lmb=config.get("Burst Size Distribution","lambda")
        return getFitBurstgen(llimit, ulimit,lmb)
    else:
        return getExpBurstgen(config.get("Burst Size Distribution","bsdfile"))    

def constructTM(options,config):
    burstgen = getBurstGenerator(config)
    tm = config.get("Transfer Matrix","type")
    if tm=="global":
        return constructGlobalTM(options,config,burstgen)
    elif tm=="local":
        return constructLocalTM(options,config,burstgen)
    elif tm=="none":
        return constructNoAVGTM(options,config,burstgen)

def resolveDistances(config, TM, effhist):
    return GaussianRegularizationDistanceReconstruction(config, TM, effhist)

def writeDistances(xvrange,distances,options):
    arr=numpy.array((xvrange,distances))
    numpy.savetxt(options.outdistfile,arr.T)
