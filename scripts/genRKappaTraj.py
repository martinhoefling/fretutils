#!/usr/bin/env python

from optparse import OptionParser, OptionGroup, SUPPRESS_HELP
import sys, os, numpy
from math import *
import random

program_version = "1.0"


def getCmdlineOptions():
    """Commandline parser, see -- long description for each option"""
    parser = OptionParser()

    group1 = OptionGroup(parser, "Required input options")
    group1.add_option("-o", "--outputfile", dest = "outputfile",
                      help = "R-Kappa output file", default = "rkappagen.dat")
    group1.add_option("-t", "--timestep", dest = "timestep",
                      help = "timestep in ps", default = 10, type = "float")
    group1.add_option("-s", "--samples", dest = "samples",
                      help = "number of samples", default = 10000, type = "int")
    group1.add_option("-r", "--distance", dest = "distexpr",
                      help = "expression to calculate distance, use t for time", default = "5.4;R", type = "string")
    group1.add_option("-c", "--kappacorr", dest = "kappacorr",
                      help = "correlation of kappa^2 and distance, e.g. kappa2*(sin)", default = "kappa2", type = "string")


    group2 = OptionGroup(parser, "Optional input")
    group2.add_option("-d", "--decay", dest = "decay",
                      help = "orientation decay per timestep", default = 1, type = "float")

    parser.add_option_group(group1)
    parser.add_option_group(group2)

    return parser.parse_args()

def genRandomVec():
    gauss = numpy.random.normal(size = 3)
    return gauss / sqrt((gauss ** 2).sum())

def getKappa(donor, acceptor, distance):
    return (donor * acceptor).sum() - 3 * (donor * distance).sum() * (acceptor * distance).sum()

def alterVec(options, vec):
    vec = (1 - options.decay) * vec + options.decay * genRandomVec()
    vec /= sqrt((vec ** 2).sum())
    return vec

def runGenTraj(options):
    donor = genRandomVec()
    acceptor = genRandomVec()
    distance = genRandomVec()
    with open(options.outputfile, "w") as fh:
        R = float(options.distexpr.split(";")[0])
        print "Starting distance is %f" % R
        for time in numpy.arange(0, options.timestep * options.samples, options.timestep):
            donor = alterVec(options, donor)
            acceptor = alterVec(options, acceptor)
            distance = alterVec(options, distance)
            kappa2 = getKappa(donor, acceptor, distance) ** 2
            t = float(time)
            R = eval(options.distexpr.split(";")[1])
            kappa2 = eval(options.kappacorr)
            fh.write("%f %f %f\n" % (time, R, kappa2))
    print "Trajectory generated."


def main():
    (options, args) = getCmdlineOptions()
    runGenTraj(options)
    print """
# %s - %s
# (C) 2012 Martin Hoefling and Helmut Grubmueller
#
# Please cite the usage as:
# Hoefling M, Lima N, Haenni D, Seidel CAM, Schuler B, Grubmueller H,  
# Structural Heterogeneity and Quantitative FRET Efficiency Distributions 
# of Polyprolines through a Hybrid Atomistic Simulation and Monte Carlo 
# Approach. (2011) PLoS ONE 6(5): e19791. doi:10.1371/journal.pone.0019791
    """ % (os.path.basename(sys.argv[0]), program_version)

if __name__ == "__main__":
    main()
