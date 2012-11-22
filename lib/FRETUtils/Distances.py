'''
Created on 13.11.2012

@author: mhoefli
'''

from FRETUtils.Bursts import getBurstSizes
from FRETUtils.Efficiencies import getBurstSizeGenerator
from FRETUtils.Ensemble import pickFromEnsemble, getTrajClassProbability
from FRETUtils.Trajectories import getRandomTrajectory
import numpy
import random

def normalizeTrajProbForSpecies(trajs, probs):
    for trajkey in trajs:
        traj = trajs[trajkey]
        if not traj.has_key("l"):
            continue
        print "Correcting trajectory", trajkey, "for class probability in 4th column."
        prb = getTrajClassProbability(traj, probs)
        traj["l"] *= prb


def generateBurst(trajs, eprob, conf, burstsize):
    """calculates efficiency according to a single given burstsize length using trajectories, configuration and probabilities and given method."""
    # select method

    # use all trajectories per burstsize
    if conf.get("Burst Accumulation", "method") == "all":
        return generateBurstFromAllTraj(eprob, trajs, conf, burstsize)

    # select species once but allow all trajectories
    elif conf.get("Burst Accumulation", "method") == "same-species":
        return generateBurstFromSameClassTraj(eprob, trajs, conf, burstsize)

    # select species and trajectory once
    elif conf.get("Burst Accumulation", "method") == "trajectory":
        return generateBurstFromSingleTraj(eprob, trajs, conf, burstsize)


def getRandomDistanceBasedOnOccupancy(traj):
    if not traj.has_key("cumsum"):
        cs = traj["l"].cumsum()
        if cs[-1] == 0:
            raise ValueError("Trajectory %s has invalid probabilities in column 4" % (traj["fname"]))
        cs /= cs[-1]
        traj["cumsum"] = cs
    ndx = traj["cumsum"].searchsorted(random.random())
    return traj["R"][ndx]

def getRandomDistance(traj):
    if traj.has_key("l"):
        return getRandomDistanceBasedOnOccupancy(traj)
    else:
        return random.choice(traj["R"])

def generateBurstFromAllTraj(eprob, trajs, conf, burstsize):
    distances = []
    while True:
        species = pickFromEnsemble(eprob)
        traj = getRandomTrajectory(trajs, species)
        distances.append(getRandomDistance(traj))

        if len(distances) == burstsize:
            return numpy.array(distances).mean()


def generateBurstFromSameClassTraj(eprob, trajs, conf, burstsize):
    species = pickFromEnsemble(eprob)
    distances = []
    while True:
        traj = getRandomTrajectory(trajs, species)
        distances.append(getRandomDistance(traj))
        if len(distances) == burstsize:
            return numpy.array(distances).mean()


def generateBurstFromSingleTraj(eprob, trajs, conf, burstsize):
    species = pickFromEnsemble(eprob)
    traj = getRandomTrajectory(trajs, species)
    distances = []
    while True:
        distances.append(getRandomDistance(traj))
        if len(distances) == burstsize:
            return numpy.array(distances).mean()

def getDistanceBursts(trajectories, probabilities, config):
    print "Generating Bursts"
    burstGenerator = getBurstSizeGenerator(config, 1)
    nbursts = config.get("Burst Size Distribution", "nbursts")
    burstsizelist = getBurstSizes(nbursts, burstGenerator)
    burstdists = []
    for bs in burstsizelist:
        burstdists.append(generateBurst(trajectories, probabilities, config, bs))
        if (len(burstdists) % 100 == 0):
            print "\r %d of %d completed." % (len(burstdists), len(burstsizelist)),

    return numpy.array(burstdists)
