from FRETUtils.Run import runReconstruction
from optparse import OptionParser, OptionGroup, SUPPRESS_HELP
import sys,os

program_version="1.0"
    

def getCmdlineOptions():
    """Commandline parser, see -- long description for each option"""
    parser = OptionParser()
    
    group1 = OptionGroup(parser, "Required input options")
    group1.add_option("-e", "--efficiencies", dest="efficiencyfile",
                      help="file with burst efficiencies", default="effs.txt")
    group1.add_option("-t", "--transferMatrix", dest="transferMatrix",
                      help="transfer Matrix type from global kappa^2 average (global), distance dependent kappa^2 average (local) or no averaging (none)", default="global")
    group1.add_option("--R0", dest="R0",
                      help="Foerster Radius including kappa=2/3", type="float", default=None, metavar=5.4)
    group1.add_option("-b", "--burst-count", dest="burstCount",
                      help="Burst count per distance bin", type="int", default=40)
    
    group2 = OptionGroup(parser, "Optional input")
    group2.add_option("--bs", "--burst-minsize", dest="burstMinsize",
                      help="Minimal burstsize (when no experimental burstsize file is provided)", type="int", default=20)
    group2.add_option("--be", "--burst-maxsize", dest="burstMaxsize",
                      help="Maximal burstsize (when no experimental burstsize file is provided)", type="int", default=100)
    group2.add_option("--ba", "--burst-exponent", dest="burstLambda",
                      help="Analytical burstsize function exponent", type="float", default=-2.3)    
    group2.add_option("-k", "--expbursts",dest="expbfile", 
                      help="experimental bursts size distribution file", default = None,metavar="exp.dat")    
    group2.add_option("-r", "--rkappa", dest="rkappafile",
                      help="file with r kappa (and probabilities)", default="rkappaprb.txt")    
    group2.add_option("--eb", "--efficiency-bins",dest="efficiencybins", 
                      help="efficiency bin count", default = 50, type="int")
    group2.add_option("--rb", "--distance-bins",dest="distancebins", 
                      help="distance bin count", default = 100, type="int")
    group2.add_option("--rs", "--distance-start",dest="distancestart", 
                      help="distance range start", default = None, type="float")
    group2.add_option("--re", "--distance-end",dest="distanceend", 
                      help="distance range end", default = None, type="float")

    
    group2.add_option("-s", "--seed",dest="rseed", 
                      help="random number generator seed", default = None, type="int" ,metavar="python_default")

    group3 = OptionGroup(parser, "Output options")
    group3.add_option("-o", "--output-distances",dest="outdistfile", 
                      help="distance histogram output file", default = "distout.txt")
   
    parser.add_option_group(group1)
    parser.add_option_group(group2)
    parser.add_option_group(group3)

    return parser.parse_args()

def main():
    (options, args) = getCmdlineOptions() 
    runReconstruction(options)
    print """
# %s - %s
# (C) 2011 Martin Hoefling and Helmut Grubmueller
#
# Please cite the usage as:
# Hoefling M, Lima N, Haenni D, Seidel CAM, Schuler B, Grubmueller H,  
# Structural Heterogeneity and Quantitative FRET Efficiency Distributions 
# of Polyprolines through a Hybrid Atomistic Simulation and Monte Carlo 
# Approach. (2011) PLoS ONE 6(5): e19791. doi:10.1371/journal.pone.0019791
    """%(os.path.basename(sys.argv[0]),program_version)

if __name__ == "__main__":
    main()