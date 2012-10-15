'''
Created on 24.06.2010

@author: mhoefli
'''

from numpy import load as npload
from numpy import zeros
#from FRETUtils.Ensemble import getTrajClassProbability
from os.path import join as pathjoin
import os
import random
import numpy

def createTrajectoryList(rootdir,fformat):
    """creates a dictionary of all npz files in a given directory"""
    trajs={}
    for root,dirs,files in os.walk(rootdir,followlinks=True):
        for filen in files:
            if filen.endswith(".%s"%fformat):
                trajs[pathjoin(root,filen)]=""
    return trajs

def loadSingleTraj(fname,fformat):
    """loads a single trajectory depending on the format chosen"""
    if fformat=="npz":
        ar=npload(fname)
        return ar['arr_0']
    elif fformat=="dat":
        return numpy.loadtxt(fname)
    else:
        raise ValueError("Invalid file format %s"%fformat)
    
def readTrajs(trajs,fformat):
    """populates dictionary by reading in trajectory files"""
    keys=trajs.keys()
    keys.sort()
    for key in keys:
        print "Reading trajectory \"%s\"."%key
        trajdict={}
        fullarr=loadSingleTraj(key,fformat)
        trajdict["t"]=fullarr[:,0]
        trajdict["R"]=fullarr[:,1]
        trajdict["k2"]=fullarr[:,2]
        if fullarr.shape[1]>3:
            trajdict["l"]=fullarr[:,3]
        trajdict["length"]=trajdict["t"].shape[0]
        print "-> %d samples read."%trajdict["length"]
        trajs[key]=trajdict
    print

def calcFRETRates(trajs,config):
    """calculate fret rates of each trajectory according to given configuration"""
    deltat=config.getfloat("Monte Carlo","deltat")
    keys=trajs.keys()
    keys.sort()
    R0=config.getfloat("FRET Constants","R0")/pow(config.getfloat("FRET Constants","kappa"),1./6)
    kDtot=config.getfloat("Dye Constants","kD")+config.getfloat("Dye Constants","kDi")
    print "Converting R/Kappa^2 trajectories to donor decay probability trajectories..."
    print "Constants for the donor decay calculation:"
    print "-> R_0 reduced by (kappa^2)^1/6 is %6.2f"%R0
    print "-> kDtot corresponding to the lifetime of %6.2e ps is %6.2e ps^-1."%(1/kDtot,kDtot)
    print "-> time increment is %f ps"%deltat
    for key in keys:
        print "Calculating donor decay rates of trajectory \"%s\"."%key
        trajs[key]["kRET"]=kDtot*trajs[key]["k2"]*(R0/trajs[key]["R"])**6
        trajs[key]["pRETpDtot"]=(trajs[key]["kRET"]+kDtot)*deltat
        
def getRandomTrajectory(trajs,species):
    """returns random trajectory of a given class and from given trajectory dictionary. It also takes the length of trajectories as weighting into account."""
    allowedkeys=[]
    maxsamples=0
    for key in trajs:
        if trajs[key]["species"]==species:
            if trajs[key]["length"] > maxsamples:
                maxsamples = trajs[key]["length"]
            allowedkeys.append((key,trajs[key]["length"]))
    
    while True:
        testkey = random.choice(allowedkeys)
        if testkey >= random.random()*maxsamples:
            return trajs[testkey[0]]  

#def getTrajWeight(trajs,key,prb):
#    """Calculates the trajectory weight in the full ensemble using the probability of the class of the trajectory and the trajectory class sample count"""
#    return getTrajClassProbability(trajs[key],prb)*(float(trajs[key]["length"])/getClassTrajSamples(trajs[key]["species"],trajs))
#    
#def getTrajWeightArray(traj,key,prb):
#    """Creates an array with trajectory ensemble weight as constant value"""
#    return zeros(len(traj[key]["R"])) + getTrajWeight(traj,key,prb) 