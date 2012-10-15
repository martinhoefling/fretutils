# -*- coding: utf-8 -*-
'''
Created on 10.06.2010

@author: martin
'''

from numpy import random
import sys
import random as pyrand

class Photon:
    def __init__(self,donor,endtime,duration,thermal=False):
        self.donor=donor
        self.thermal=thermal
        self.endtime=endtime
        self.duration=duration
        
    def checkThermal(self,QD,QA):
        #we have donor photon roll if really a photon...
        if self.donor:
            if pyrand.random() > QD:
                self.thermal=True
            else:
                self.thermal=False
                
        #we have acceptor photon roll if really a photon...    
        else:
            if pyrand.random() > QA:
                self.thermal=True
            else:
                self.thermal=False  

#def fretPotential2Rate(potential,spectraloverlap):
#    """Converts the coupling potential (in J) and a spectral overlap to transition rate"""  
#    return 1./(Constants.h_bar**2*Constants.c_light)*Helpers.JouleInNm(potential)*spectraloverlap

def getPhoton(traj,conf,startndx=None):
    """Calculates a single photon with given configuration and trajectory, returns 1 for acceptor and 0 for donor photon"""

    if not startndx:
        try:
            startndx=random.randint(conf.getint("Monte Carlo","minstarttraj"),int(len(traj["t"])-conf.getint("Monte Carlo","maxstarttraj")))
        except ValueError:
            raise IndexError("Invalid range for startindex generation %d, %d. Adjust minstarttraj and maxstarttraj parameters."%(conf.getint("Monte Carlo","minstarttraj"),len(traj["t"])-conf.getint("Monte Carlo","maxstarttraj")))
    
    if not -1 < startndx < int(len(traj["t"])):
        raise IndexError("Start index outside of trajectory length. Adjust minstarttraj and maxstarttraj parameters.")
    
    pDtot=(conf.getfloat("Dye Constants","kD")+conf.getfloat("Dye Constants","kDi"))*conf.getfloat("Monte Carlo","deltat")
    pRETpDtot=traj["pRETpDtot"]
    photon=-1
    while photon == -1:
        rejectcount=0
        while rejectcount <= conf.getint("Monte Carlo","rejectretry"):
            photon,endtime=tryGetPhoton(pDtot,pRETpDtot,startndx,random.randint(0,sys.maxint))

            if not rejectPhoton(startndx,endtime,traj,conf):
                duration=endtime-startndx

                if photon==0:
                    donor=True
                else:
                    donor=False
                return Photon(donor,endtime*conf.getfloat("Monte Carlo","deltat"),duration*conf.getfloat("Monte Carlo","deltat"))
            else:
                rejectcount+=1  
        raise ValueError("Reject count reached")
           

def rejectPhoton(starttime,endtime,traj,conf):
    """tests if the photon has to be rejected, based on given starttime, endtime, trajectory and configuration"""
    if not traj.has_key("R_min"):
        traj["R_min"]=traj["R"].min()
    if conf.getfloat("Monte Carlo","photrejectdist") < traj["R_min"]:
        return False #no rejection possible for this trajecotry, since all distances are larger 
    if traj["R"][starttime:endtime-+1].min() > conf.getfloat("Monte Carlo","photrejectdist"):
        return False
    return True

#def getPhotonPBC(traj,conf,startndx=None):
#    """calculates photons  and restarts from time 0 if at trajectory end, no photon rejection implemented here so far"""
#    if not startndx:
#        startndx=random.randint(0,int(len(traj["t"])-1))
#    pDtot=(conf.getfloat("Dye Constants","kD")+conf.getfloat("Dye Constants","kDi"))*conf.getfloat("Monte Carlo","deltat")
#    pRETpDtot=traj["pRETpDtot"]
#    photon=-1
#    while photon == -1:
#        photon,endtime=tryGetPhoton(pDtot,pRETpDtot,startndx,random.randint(0,sys.maxint))
#        if photon ==-1:
#            startndx=0
#    return photon,endtime

def tryPlainGetPhoton(p0,pvar,start,seed):
    """Try to generate a photon and returns 0 and 1 for donor and acceptor and -1 for failure"""
    for curndx in range(start,len(pvar)):
        rnd=random.random()
        if rnd < p0:
            return 0,curndx
        if rnd < pvar[curndx]:
            return 1,curndx
        
    return -1,len(pvar)

def setPhotonGenerator(config):
    """sets the photon generation routine (python or c-extension)"""
    global tryGetPhoton
    print "Determining MC routine generation routine..."
    if config.get("System", "photongenerator")=="python":
        print "-> Using photon generator in plain python"
        tryGetPhoton=tryPlainGetPhoton    
    elif config.get("System", "photongenerator")=="cextension":
        print "-> Using c-extension photon generator"
        from FRETUtils.fretnumpyext import tryGetCPhoton as tryGetPhoton
    elif config.get("System", "photongenerator")=="cython":
        print "-> Using cython photon generator"
        from FRETUtils.PhotonGenerator import tryGetCythonPhoton as tryGetPhoton        
    else:
        raise ValueError("Invalid photon generator %s"%config.get("System", "photongenerator"))

tryGetPhoton=tryPlainGetPhoton