'''
Created on 13.11.2012

@author: mhoefli
'''

from FRETUtils.Efficiencies import getBurstSizeGenerator
from FRETUtils.Trajectories import getRandomTrajectory
from FRETUtils.Ensemble import pickFromEnsemble
from FRETUtils.Bursts import getBurstSizes
import numpy,random


def generateBurst(trajs,eprob,conf,burstsize):
    """calculates efficiency according to a single given burstsize length using trajectories, configuration and probabilities and given method."""
    #select method
    
    #use all trajectories per burstsize
    if conf.get("Burst Accumulation","method") == "all":
        return generateBurstFromAllTraj(eprob,trajs,conf,burstsize)
                    
    #select species once but allow all trajectories
    elif conf.get("Burst Accumulation","method") == "same-species":
        return generateBurstFromSameClassTraj(eprob,trajs,conf,burstsize)
                               
    #select species and trajectory once
    elif conf.get("Burst Accumulation","method") == "trajectory":
        return generateBurstFromSingleTraj(eprob,trajs,conf,burstsize)


def getRandomDistance(traj):
    return random.choice(traj["R"])


def generateBurstFromAllTraj(eprob,trajs,conf,burstsize):
    distances=[]
    while True:
        species=pickFromEnsemble(eprob)
        traj=getRandomTrajectory(trajs,species)
        distances.append(getRandomDistance(traj))
            
        if len(distances)==burstsize:
            return numpy.array(distances).mean()
          

def generateBurstFromSameClassTraj(eprob,trajs,conf,burstsize):
    species=pickFromEnsemble(eprob)
    distances=[]
    while True:
        traj=getRandomTrajectory(trajs,species)
        distances.append(getRandomDistance(traj))        
        if len(distances)==burstsize:
            return numpy.array(distances).mean()


def generateBurstFromSingleTraj(eprob,trajs,conf,burstsize):
    species=pickFromEnsemble(eprob)
    traj=getRandomTrajectory(trajs,species)
    distances=[]
    while True:
        distances.append(getRandomDistance(traj))
        if len(distances)==burstsize:
            return numpy.array(distances).mean()

def getDistanceBursts(trajectories,probabilities,config):
    burstGenerator = getBurstSizeGenerator(config,1)
    nbursts = config.get("Burst Size Distribution", "nbursts")
    burstsizelist=getBurstSizes(nbursts,burstGenerator)
    burstdists=[]
    for bs in burstsizelist:
        burstdists.append(generateBurst(trajectories,probabilities,config,bs))
    
    return numpy.array(burstdists)
    