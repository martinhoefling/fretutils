'''
Created on 04.10.2012

@author: mhoefli
'''

import numpy
import ConfigParser


class FakeOptions():
    def __init__(self):
        self.trajdirectory = "."
        self.configfilename = "fret.conf"
        self.pbfile = "probabilities.dat"
        
        self.expbfile = "exp.dat"
        self.rseed = 123
        self.trajformat = "npz"
        
        self.binaryofile = None
        self.efficiencyofile= None
        self.burstcompofile = None
        self.burstsizeofile = None
        self.endtimeofile = None
        self.decaytimeofile = None
        self.prffile = None

class FakeOptionsRK():
    def __init__(self):
        self.trajdirectory = "."
        self.configfilename = None
        self.pbfile = "probabilities.dat"
        self.rseed = 123
        self.trajformat = "npz"
        
        self.outtrajfile = "trajout.txt"

def createConstantDummyRKTraj(filename,nrsamples,dt,Rval,kappaval,fformat="numpy"):
    time = numpy.linspace(0,nrsamples*dt,nrsamples)
    R = numpy.ones(nrsamples)*Rval
    Kappa = numpy.ones(nrsamples)*kappaval
     
    comb = numpy.array((time,R,Kappa)).T
    if fformat is "numpy":
        numpy.savez(filename,comb)       
    elif fformat is "plain":
        numpy.savetxt(filename,comb)   
    else:
        raise ValueError("Invalid file format.",fformat)

def createStepFunctionDummyRKTrajWithL(filename,nrsamples,dt,Rval1,Rval2,kappaval,fformat="numpy"):
    time = numpy.linspace(0,nrsamples*dt,nrsamples)
    R = numpy.ones(nrsamples)*Rval1
    hlength=R.shape[0]/2
    R[hlength:]=Rval2
    Kappa = numpy.ones(nrsamples)*kappaval
     
    comb = numpy.array((time,R,Kappa,R)).T
    if fformat is "numpy":
        numpy.savez(filename,comb)       
    elif fformat is "plain":
        numpy.savetxt(filename,comb)   
    else:
        raise ValueError("Invalid file format.",fformat)

def createProbabilityClassFile(fname,classes,prbs):
    try:
        f=open(fname,"w")
        for i in range(len(classes)):
            f.write(".*%s.* %s %6.4f\n"%(classes[i],classes[i],prbs[i]))
    finally:
        f.close()     

def writeBurstFile():
    try:
        f=open("expbursts.bst","w")
        f.write(bursts)
    finally:
        f.close()   

def writeInvalidDummyFile():
    try:
        f=open("invalid.inv","w")
        f.write("invalid")
    finally:
        f.close()   


def writeInvalidBurstFile():
    try:
        f=open("expburstsinvalid.bst","w")
        f.write(bursts)
        f.write("1,2,3,4,5\n")
    finally:
        f.close()   
        
def writeConfigFiles():
    try:
        fh=open("standard.conf","w")
        fh.write(configfile)
    finally:
        fh.close()
    
    config = ConfigParser.RawConfigParser()
    config.read('standard.conf')
    
    config.set('Burst Accumulation','method','same-species')
    try:
        fh=open("same-species.conf","w")
        config.write(fh)
    finally:
        fh.close()    

    config.set('Burst Accumulation','method','invalid')
    try:
        fh=open("invalidburstacc.conf","w")
        config.write(fh)
    finally:
        fh.close()    


    config.set('Burst Accumulation','method','all')
    try:
        fh=open("all.conf","w")
        config.write(fh)
    finally:
        fh.close()    

    config.set('System','ncpu','2')
    try:
        fh=open("dual.conf","w")
        config.write(fh)
    finally:
        fh.close()            

    config.set('System','ncpu','-1')
    config.set('System','blocksize',"2")
    try:
        fh=open("multi.conf","w")
        config.write(fh)
    finally:
        fh.close()    

    config.set('Dye Constants','QD','0.5')
    config.set('Dye Constants','QA','0.5')
    config.set('System','ncpu','1')
    try:
        fh=open("thermal.conf","w")
        config.write(fh)
    finally:
        fh.close()    



    config.set('Burst Size Distribution','apply','invalid')
    try:
        fh=open("invalidbscutoff.conf","w")
        config.write(fh)
    finally:
        fh.close()             


    config.set('Burst Size Distribution','apply','corrected')
    try:
        fh=open("corrected.conf","w")
        config.write(fh)
    finally:
        fh.close()             

    config.set('Burst Size Distribution','method','invalidbsd')
    try:
        fh=open("invalidbsd.conf","w")
        config.write(fh)
    finally:
        fh.close() 

    config.set('Burst Size Distribution','method','file')
    try:
        fh=open("expbstcorr.conf","w")
        config.write(fh)
    finally:
        fh.close()  
                  
    config.set('Burst Size Distribution','method','file')
    config.set('Burst Size Distribution','apply','true-photon')
    try:
        fh=open("expbst.conf","w")
        config.write(fh)
    finally:
        fh.close()  
    
    config.set('Monte Carlo','minstarttraj','110000')
    config.set('Monte Carlo','maxstarttraj','110000')
    try:
        fh=open("wrongminmaxstart.conf","w")
        config.write(fh)
    finally:
        fh.close()
        
    config.set('Monte Carlo','minstarttraj','-100000')
    config.set('Monte Carlo','maxstarttraj','100000')
    try:
        fh=open("wrongminmaxstart2.conf","w")
        config.write(fh)
    finally:
        fh.close()   

    config.set('Monte Carlo','minstarttraj','0')
    config.set('Monte Carlo','maxstarttraj','1000')
    config.set('Monte Carlo','photrejectdist','10')
    config.set('Monte Carlo','globalrejectretry','10')

    for i in ("trajectory","same-species","all"):
        config.set('Burst Accumulation','method',i)
        try:
            fh=open("rejecttest-%s.conf"%i,"w")
            config.write(fh)
        finally:
            fh.close()  



    config.set('Monte Carlo','globalrejectretry','100000')
    config.set('Monte Carlo','photrejectdist','2')
    config.set("Monte Carlo","nbursts","5")

    try:
        fh=open("stepreject.conf","w")
        config.write(fh)
    finally:
        fh.close()  

    config.set("Monte Carlo","nbursts","50")
    config.set('System','photongenerator','cython')
    try:
        fh=open("cython.conf","w")
        config.write(fh)
    finally:
        fh.close()  

    config.set('System','photongenerator','cextension')
    try:
        fh=open("cextension.conf","w")
        config.write(fh)
    finally:
        fh.close()  

    config.set('System','photongenerator','invalid')
    try:
        fh=open("invalid_photongenerator.conf","w")
        config.write(fh)
    finally:
        fh.close()  

def writeRKConfigFiles():
    try:
        fh=open("standard.conf","w")
        fh.write(configfile)
    finally:
        fh.close()
    
    config = ConfigParser.RawConfigParser()
    config.read('standard.conf')
    config.set('Monte Carlo','minstarttraj','0')
    config.set('Monte Carlo','maxstarttraj','0')
    try:
        fh=open("standardRK.conf","w")
        config.write(fh)
    finally:
        fh.close()  
    
    config.set('System','ncpu','-1')
    config.set('Photon Flooding','photoncount','20')

    try:
        fh=open("standardRKmulti.conf","w")
        config.write(fh)
    finally:
        fh.close() 
    
    config.set('System','ncpu','1')
    config.set('Photon Flooding','photoncount','5')
    config.set('Photon Flooding','startclip','20')
    config.set('Photon Flooding','endclip','40')
    try:
        fh=open("standardRKclip.conf","w")
        config.write(fh)
    finally:
        fh.close()  


configfile="""[Dye Constants]
tauD=4000
tauA=3900
QD=1.00
QA=1.00

[FRET Constants]
R0 = 4.0
kappa = 0.66666666

[Burst Size Distribution]
method = analytical
llimit = 20
ulimit = 80
lambda = -2.3
apply = true-photon

[Burst Accumulation]
method = trajectory
#method = same-species
#method = all

[Monte Carlo]
minstarttraj = 0
maxstarttraj = 1000
deltat = 10
nbursts = 50
photrejectdist = 0.0
rejectretry = 10

[System]
photongenerator = python
ncpu =  1

[Photon Flooding]
"""

bursts="""6,34,1.395339993,38.22816933
12,119,-0.4387392944,144.7120697
15,8,14.17343894,8.605908493
26,24,21.46747034,23.48822066
75,53,67.4385116,59.3010967
12,50,6.09756327,59.69610025
144,105,134.0717494,125.476435
125,77,118.5791254,91.65645018
104,44,100.2808946,51.92793861
43,9,41.46921249,7.985643681
71,79,59.85444839,89.39337971
16,21,12.09963608,21.06266482
9,55,2.21971254,64.62174591
14,79,5.816911471,96.95722543
60,45,52.02342248,44.98337885
10,16,6.620163838,16.07181178

19,7,16.58802966,4.199380926
38,13,34.30130912,9.820283086
41,22,35.60566827,18.33920318
53,18,50.7423003,19.22946756
115,84,103.0088362,91.88876314
17,13,13.54572318,11.05131973
41,31,36.65179876,33.01846504
37,21,32.72479647,19.93417479

"""