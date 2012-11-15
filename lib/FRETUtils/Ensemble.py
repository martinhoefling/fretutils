'''
Created on 24.06.2010

@author: mhoefli
'''

import random
import re

def readProbabilities(pbfile):
    """Reads in probabilities from a specific probability file, 1st colum is the regexp to match, 2nd the name and 3rd the probability in the ensemble"""
    with open(pbfile) as pfh:
        probabilities = []
        print
        for line in pfh:
            spl = line.split()
            if len(spl) == 0:
                continue

            if len(spl) != 3:
                print "Line split is :", spl
                raise ValueError("Line in probability file has not 3 entries.")
            probabilities.append((spl[1], re.compile(spl[0]), float(spl[2])))
            print "Found ensemble class", spl[1], "with probability %6.4f." % (float(spl[2]))
        print
    return probabilities

def assignTrajProbabilityClasses(trajs, probabilities):
    """applies classes to a trajectory dictionary with read in and compiled probabilities"""
    keys = trajs.keys()
    keys.sort()
    for key in keys:
        for pclass in probabilities:
            myclass = None
            if pclass[1].search(key):
                myclass = pclass[0]
                break
        if not myclass:
            raise ValueError("Cannot assign a probability class to simulation %s" % key)
        trajs[key]["species"] = myclass
        print "Assigned species", myclass, "to trajectory \"%s\"." % key

def getClassTrajCount(myclass, trajs):
    """counts trajectories of a distinct class"""
    counter = 0
    for trajk in trajs.keys():
        if trajs[trajk]["species"] == myclass:
            counter += 1
    return counter

def getClassTrajSamples(myclass,trajs):
    """counts the samples from all trajectories of a distinct class"""
    counter=0
    for trajk in trajs.keys():
        if trajs[trajk]["species"] == myclass:
            counter+=trajs[trajk]["length"]
    return counter

def getClassTrajWeight(key,trajs):
    myclass = trajs[key]["species"]
    classsamples = getClassTrajSamples(myclass,trajs)
    return float(trajs[key]["length"])/classsamples

def cleanProbabilities(trajs, probs):
    """removes all probability classes from probability dictionary which do not have any trajectory"""
    newprobs = []
    for prob in probs:
        if getClassTrajCount(prob[0], trajs) > 0:
            newprobs.append(prob)
    return newprobs


def pickFromEnsemble(eprob):
    """returns a class from ensemble according to its probability, this function (re-) normalizes the probability if necessary"""
    epsum = 0.
    for pr in eprob:
        epsum += pr[2]

    rnd = random.random() * epsum
    epsum = 0
    for pr in eprob:
        epsum += pr[2]
        if rnd < epsum:
            return pr[0]

def getTrajClassProbability(traj, probabilities):
    """returns the class probability of a distinct trajectory"""
    for pr in probabilities:
        if pr[0] == traj["species"]:
            return pr[2]


