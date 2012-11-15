'''
Created on 24.06.2010

@author: mhoefli
'''

from numpy import load as npload
from numpy import zeros
from FRETUtils.Ensemble import getTrajClassProbability, getClassTrajCount
from FRETUtils.Photons import getPhoton
from os.path import join as pathjoin
import os, sys
import random
import numpy

def writeRKProbTraj(fh, trajs, probabilities, config):
    """writes time, distance, kappa and the total probability in an opened file handle"""
    if config:
        startclip = config.get("Photon Flooding", "startclip")
        endclip = config.get("Photon Flooding", "endclip")
    else:
        startclip = 0
        endclip = 0

    keys = trajs.keys()
    keys.sort()
    for key in keys:
        fh.write("# starting traj %s\n" % key)

        classprob = getTrajClassProbability(trajs[key], probabilities)
        nrtrajs = getClassTrajCount(trajs[key]["species"], trajs)
        probprefact = classprob / nrtrajs

        if not trajs[key].has_key("photons"):
            trajs[key]["photons"] = numpy.ones(len(trajs[key]["t"]))
        trajs[key]["photons"] = trajs[key]["photons"] * probprefact
        arr = numpy.array((trajs[key]["t"], trajs[key]["R"], trajs[key]["k2"], trajs[key]["photons"]))
        if endclip == 0:
            numpy.savetxt(fh, arr.T[startclip:, :])
        else:
            numpy.savetxt(fh, arr.T[startclip:-endclip, :])

        print key, "written with traj %s probability with length %d, clipping %d and %d frames at beginning and end." % (probprefact, trajs[key]["length"], startclip, endclip)

def createTrajectoryList(rootdir, fformat):
    """creates a dictionary of all npz files in a given directory"""
    trajs = {}
    for root, _dirs, files in os.walk(rootdir, followlinks = True):
        for filen in files:
            if filen.endswith(".%s" % fformat):
                trajs[pathjoin(root, filen)] = ""
    return trajs

def loadSingleTraj(fname, fformat):
    """loads a single trajectory depending on the format chosen"""
    if fformat == "npz":
        ar = npload(fname)
        return ar['arr_0']
    elif fformat == "dat":
        return numpy.loadtxt(fname)
    else:
        raise ValueError("Invalid file format %s" % fformat)

def readTrajs(trajs, fformat):
    """populates dictionary by reading in trajectory files"""
    keys = trajs.keys()
    keys.sort()
    for key in keys:
        print "Reading trajectory \"%s\"." % key
        trajdict = {}
        fullarr = loadSingleTraj(key, fformat)
        trajdict["t"] = fullarr[:, 0]
        trajdict["R"] = fullarr[:, 1]
        trajdict["k2"] = fullarr[:, 2]
        if fullarr.shape[1] > 3:
            trajdict["l"] = fullarr[:, 3]
        trajdict["length"] = trajdict["t"].shape[0]
        print "-> %d samples read." % trajdict["length"]
        trajs[key] = trajdict
    print

def calcFRETRates(trajs, config):
    """calculate fret rates of each trajectory according to given configuration"""
    deltat = config.get("Monte Carlo", "deltat")
    keys = trajs.keys()
    keys.sort()
    R0 = config.get("FRET Constants", "R0") / pow(config.get("FRET Constants", "kappa"), 1. / 6)
    kDtot = config.get("Dye Constants", "kD") + config.get("Dye Constants", "kDi")
    print "Converting R/Kappa^2 trajectories to donor decay probability trajectories..."
    print "Constants for the donor decay calculation:"
    print "-> R_0 reduced by (kappa^2)^1/6 is %6.2f" % R0
    print "-> kDtot corresponding to the lifetime of %6.2e ps is %6.2e ps^-1." % (1 / kDtot, kDtot)
    print "-> time increment is %f ps" % deltat
    for key in keys:
        print "Calculating donor decay rates of trajectory \"%s\"." % key
        trajs[key]["kRET"] = kDtot * trajs[key]["k2"] * (R0 / trajs[key]["R"]) ** 6
        trajs[key]["pRETpDtot"] = (trajs[key]["kRET"] + kDtot) * deltat

def getRandomTrajectory(trajs, species):
    """returns random trajectory of a given class and from given trajectory dictionary. It also takes the length of trajectories as weighting into account."""
    allowedkeys = []
    maxsamples = 0
    for key in trajs:
        if trajs[key]["species"] == species:
            if trajs[key]["length"] > maxsamples:
                maxsamples = trajs[key]["length"]
            allowedkeys.append((key, trajs[key]["length"]))

    while True:
        testkey = random.choice(allowedkeys)
        if testkey[1] >= random.random() * maxsamples:
            return trajs[testkey[0]]

def floodTrajsWithPhotons(trajs, config, randseed):
    random.seed(randseed)
    numpy.random.seed(random.randint(0, sys.maxint))
    photcount = config.get("Photon Flooding", "photoncount")

    deltat = config.get("Monte Carlo", "deltat")
    verbose = config.get("System", "verbose")


    for key in trajs:
        traj = trajs[key]
        if verbose:
            print "Processing trajectory", key
        traj["photons"] = zeros(traj["length"])
        for ndx in range(traj["length"]):
            if verbose and ndx % 10 == 0:
                print "%d/%d\r" % (ndx, traj["length"]),

            for _ in range(photcount):
                try:
                    photon = getPhoton(traj, config, ndx)
                    endndx = int(photon.endtime / deltat)
                    if endndx != traj["length"]:
                        traj["photons"][endndx] += 1
                except ValueError:
                    pass


    return trajs
