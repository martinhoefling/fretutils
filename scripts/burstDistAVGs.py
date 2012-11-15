#!/usr/bin/env python


from FRETUtils.Run import runBurstDistAVGs
from optparse import OptionParser, OptionGroup, SUPPRESS_HELP
import sys, os

program_version = "1.0"


def getCmdlineOptions():
    """Commandline parser, see -- long description for each option"""
    parser = OptionParser()

    group1 = OptionGroup(parser, "Required input options")
    group1.add_option("-d", "--directory", dest = "trajdirectory",
                      help = "directory with R-Kappa trajectories", default = ".", metavar = "RKDIR")
    group1.add_option("-p", "--probabilites", dest = "pbfile",
                      help = "definition and probability of trajectory classes", default = "probabilities.dat", metavar = "probabilities.dat")
    group1.add_option("-c", "--configfile", dest = "configfilename",
                      help = "configuration filename, default file will be written if file does not exist", default = "distAVG.conf")


    group2 = OptionGroup(parser, "Optional input")
    group2.add_option("-k", "--expbursts", dest = "expbfile",
                      help = "experimental bursts size distribution file", default = None, metavar = "exp.dat")
    group2.add_option("-s", "--seed", dest = "rseed",
                      help = "random number generator seed", default = None, type = "int" , metavar = "python_default")
    group2.add_option("-r", "--trajformat", default = "npz", dest = "trajformat",
                      help = "trajectory format: npz (numpy), dat (plaintext)")

    group3 = OptionGroup(parser, "Output options")
    group3.add_option("-o", "--distanceburst-out", dest = "distoutfile",
                      help = "distance burst output file", default = "distout.txt", metavar = "distout.txt")

    parser.add_option_group(group1)
    parser.add_option_group(group2)
    parser.add_option_group(group3)

    return parser.parse_args()

def main():
    (options, args) = getCmdlineOptions()
    runBurstDistAVGs(options)
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
