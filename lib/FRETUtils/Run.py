'''
Created on 04.10.2012

@author: mhoefli
'''

from FRETUtils.Efficiencies import calculateBursts,calcKineticRatesFromConfig
from FRETUtils.Ensemble import readProbabilities,assignTrajProbabilityClasses,cleanProbabilities
from FRETUtils.Photons import setPhotonGenerator
from FRETUtils.Trajectories import createTrajectoryList,readTrajs,calcFRETRates,writeRKProbTraj,floodTrajsWithPhotons
from FRETUtils.Config import FRETConfigParser,ReconstructionConfigParser,BurstDistAVGConfigParser
from FRETUtils.Reconstruction import readEfficiencies, constructTM,\
    resolveDistances, writeDistances
from FRETUtils.Distances import getDistanceBursts


import random
import sys
import cPickle
import multiprocessing
from numpy import array,savetxt
import os

def getSecureConfig(conffile,parser):
    config = parser()
    if not os.path.exists(conffile):
        with open(conffile,"w") as configfile:
            config.write(configfile)
    config.readfp(open(conffile))
    return config

def getFRETConfig(conffile):
    return getSecureConfig(conffile,FRETConfigParser)

def getReconstructionConfig(conffile):
    return getSecureConfig(conffile,ReconstructionConfigParser)

def getDistAVGConfig(conffile):
    return getSecureConfig(conffile,BurstDistAVGConfigParser)

def doMultiprocessRun(options, config, trajectories, eprobabilities):
    ncpu = config.get("System", "ncpu")

    print "Preparing %d clients for burst generation." % ncpu
    pool = multiprocessing.Pool(ncpu)
    
    results=[]
         
    nbursts = config.get("Monte Carlo", "nbursts")
    print "Total bursts to calculate:",nbursts
 
    try:
        blocksize = config.get("System", "blocksize")
    except:
        blocksize = 100
        print "Using default blocksize of", blocksize,"for each job."
 
    blockcount = nbursts // blocksize
    remainder = nbursts % blocksize
    blocks=[blocksize]*blockcount
    blocks.append(remainder)
    print "Setting up %d jobs to generate %d bursts" % (blockcount,nbursts)
    config.set("System","verbose",0)
        
    for block in blocks:
        options.rseed = random.randint(0, sys.maxint)
        config.makeReadonly()
        res= pool.apply_async(calculateBursts, (trajectories, eprobabilities, config,block,random.randint(0,sys.maxint)))
        results.append(res)
              
    bursts=[]
    
    for ene,res in enumerate(results,start=1):
        bursts+=res.get()
        print "\r%6d of %6d jobs processed."%(ene,blockcount),        
    
    return bursts

def efficienciesFromBursts(config,bursts):
    QD=config.get("Dye Constants","QD")
    QA=config.get("Dye Constants","QA")

    effs=[]
    for burst in bursts:
        effs.append(burst.getEfficiency(QD,QA))
    return array(effs)
    
def sizesFromBursts(bursts,corrected=False,QD=0.,QA=0.):
    sizes=[]
    for burst in bursts:
        if not corrected:
            sizes.append(burst.bsize)
        else:
            sizes.append(burst.bsizeCorr(QD,QA))
               
    return array(sizes)   

def writeOutputFiles(options, config, bursts):
    print 
    print "================================ Calculation complete ========================="
    print "Preparing output."
    if options.binaryofile:
        print "Binary output requested, this may take a while..."
        fh = open(options.binaryofile, "w")
        cPickle.dump(bursts, fh)
        fh.close()
        print "Binary output (pickled data) written to ", options.binaryofile
    if options.efficiencyofile:
        print "Burst efficiency output requested."
        fh = open(options.efficiencyofile, "w")
        savetxt(fh, efficienciesFromBursts(config, bursts))
        fh.close()
        print "Burst efficiencies written to ", options.efficiencyofile
    if options.burstsizeofile:
        print "Burst size output requested."
        fh = open(options.burstsizeofile, "w")
        norm = sizesFromBursts(bursts)
        corr = sizesFromBursts(bursts, corrected=True, QD=config.get("Dye Constants", "QD"), QA=config.get("Dye Constants", "QA"))
        savetxt(fh, array((norm, corr)).T)
        fh.close()
        print "Burst size written to ", options.burstsizeofile
    if options.burstcompofile:
        print "Burstcomposition output requested."
        fh = open(options.burstcompofile, "w")
        for myburst in bursts:
            fh.write("%d %d %d %d\n" % (myburst.donorphot, myburst.acceptorphot, myburst.donortherm, myburst.acceptortherm))
        
        fh.close()
        print "Burstcomposition written to ", options.burstcompofile
    if options.endtimeofile:
        print "Endtime output requested."
        fh = open(options.endtimeofile, "w")
        for burst in bursts:
            for phot in burst.photons:
                fh.write("%f " % phot.endtime)
            
            fh.write("\n")
        
        fh.close()
        print "Endtimes written to ", options.endtimeofile
    if options.decaytimeofile:
        print "Decaytime output requested."
        fh = open(options.decaytimeofile, "w")
        for burst in bursts:
            for phot in burst.photons:
                fh.write("%f " % phot.duration)
            
            fh.write("\n")
        
        fh.close()
        print "Decaytimes written to ", options.decaytimeofile
    print "Finished!"


def readTrajAndClasses(options):
    print "\nReading trajectories recursively from directory \"%s\"." % options.trajdirectory
    trajectories = createTrajectoryList(options.trajdirectory, options.trajformat)
    if len(trajectories) == 0:
        raise ValueError("No trajectories found. Did you specify the correct format with the -r switch? Exiting.")
    readTrajs(trajectories, options.trajformat)
    print "\nReading ensemble probabilities"
    eprobabilities = readProbabilities(options.pbfile)
    assignTrajProbabilityClasses(trajectories, eprobabilities)
    print "Removing empty classes"
    eprobabilities = cleanProbabilities(trajectories, eprobabilities)
    return trajectories, eprobabilities


def readConfigAndAssignFRETRate(options, trajectories):
    print "Reading configuration file \"%s\"." % options.configfilename
    config = getFRETConfig(options.configfilename)
    calcKineticRatesFromConfig(config)
    config.set("System", "verbose", 0)
    if options.rseed:
        print "Setting up RNG seed to %d" % options.rseed
        random.seed(options.rseed)
    setPhotonGenerator(config)
    calcFRETRates(trajectories, config)
    return config

def runMCFRET(options):
    trajectories, eprobabilities = readTrajAndClasses(options)   
    config = readConfigAndAssignFRETRate(options, trajectories)
    config.sethidden("Burst Size Distribution", "bsdfile", options.expbfile,str)
    print 
    print "================================ Input prepared ========================="
    print 
    print "Starting efficiency calculation"
    ncpu = config.get("System", "ncpu")
    if ncpu==-1:
        print "Determining number of available cpu's..."
        print "-> %d cpus detected" % multiprocessing.cpu_count()
        ncpu=multiprocessing.cpu_count()
        config.set("System", "ncpu","%d"%ncpu)
    
    if options.prffile:
        print "THIS IS A PROFILING RUN! - will write to logfile and run with only process", options.prffile
        import cProfile
        cProfile.runctx('bursts = calculateBursts(trajectories,eprobabilities,config,%d,%d)'%(config.get("Monte Carlo", "nbursts"),random.randint(0,sys.maxint)), globals(), locals(), options.prffile)
        print "Profiling runs write not output..."
        
    elif ncpu > 1:
        bursts = doMultiprocessRun(options, config, trajectories, eprobabilities)
    
    else:
        print "Doing single process run."
        config.set("System","verbose",1)
        print "Will calculate efficiencies from",config.get("Monte Carlo", "nbursts"),"bursts."
        print "Setting up burst generator."
        bursts = calculateBursts(trajectories, eprobabilities, config,config.get("Monte Carlo", "nbursts"),random.randint(0,sys.maxint))
    
    if not options.prffile:
        writeOutputFiles(options, config, bursts)

def runMultiprocessPhotonFlooding(trajectories,config):
    ncpu = config.get("System", "ncpu")

    print "Preparing %d clients for burst generation." % ncpu
    pool = multiprocessing.Pool(ncpu)

    results=[]    
    print "Setting up jobs for %d trajectories." % (len(trajectories))
    config.set("System","verbose",0)
    for traj in trajectories:
        trajs={}
        trajs[traj]=trajectories[traj]
        config.makeReadonly()
        res= pool.apply_async(floodTrajsWithPhotons, (trajs, config,random.randint(0,sys.maxint)))
        results.append(res)
              
    bursts=[]
    
    for ene,res in enumerate(results,start=1):
        restrajs = res.get()
        for restraj in restrajs:
            trajectories[restraj]=restrajs[restraj]
        print "\r%6d of %6d jobs processed."%(ene,len(trajectories)),        
    
    return bursts    

def runTrajPhotonFlooding(trajectories,config):
    print 
    print "================================ Input prepared ========================="
    print 
    print "Starting trajectory processing"
    ncpu = config.get("System", "ncpu")
    if ncpu==-1:
        print "Determining number of available cpu's..."
        print "-> %d cpus detected" % multiprocessing.cpu_count()
        ncpu=multiprocessing.cpu_count()
        config.set("System", "ncpu",ncpu)
           
    if ncpu > 1:
        runMultiprocessPhotonFlooding(trajectories,config)
    
    else:
        config.set("System","verbose",1)
        print "Doing single process run."
        trajectories=floodTrajsWithPhotons(trajectories,config,random.randint(0,sys.maxint))    

def runTrajPrbAdd(options):
    trajectories, eprobabilities = readTrajAndClasses(options)
    if options.configfilename:
        config = readConfigAndAssignFRETRate(options, trajectories)
        runTrajPhotonFlooding(trajectories,config)
    else:
        print "No config file specified, assigning class probabilities only."
        config = None
        
    print "Writing to",options.outtrajfile
    try:
        fh=open(options.outtrajfile,"w")
        writeRKProbTraj(fh,trajectories,eprobabilities,config)
    finally:
        fh.close()

def validateOptions(options):
    if not options.efficiencyfile:
        print "No efficiency file (-e) specified."
        sys.exit(1)
    if not os.path.exists(options.efficiencyfile):
        print "Efficiency file %s does not exist."%options.efficiencyfile
        sys.exit(1)

def runReconstruction(options):
    config = getReconstructionConfig(options.configfile)
    config.sethidden("Burst Size Distribution", "bsdfile", options.expbfile,str)

    if options.rseed:
        print "Setting up RNG seed to %d" % options.rseed
        random.seed(options.rseed)
    effhist=readEfficiencies(options.efficiencyfile,config.get("Transfer Matrix","efficiency bins"))
    TM = constructTM(options,config)
    r_prdist,xrange,e_fitprdist,fitvals = resolveDistances(config,TM,effhist)
    writeDistances(xrange,r_prdist,options)

def runBurstDistAVGs(options):
    if options.rseed:
        print "Setting up RNG seed to %d" % options.rseed
        random.seed(options.rseed)
    trajectories, eprobabilities = readTrajAndClasses(options)   
    config.sethidden("Burst Size Distribution", "bsdfile", options.expbfile,str)
    config = getDistAVGConfig(options.configfilename)
    distbursts = getDistanceBursts(trajectories,eprobabilities,config)
    if options.distoutfile:
        print "Burst efficiency output requested."
        fh = open(options.distoutfile, "w")
        savetxt(fh, distbursts)
        fh.close()
    
    