'''
Created on 26.10.2012

@author: mhoefli
'''

from ConfigParser import ConfigParser

class SecureConfigParser(object,ConfigParser):
    def __init__(self):
        ConfigParser.__init__(self)
        self.allowed={}
        self.errormsg={}
#        self.optionxform = str
        
    def setdefault(self,section,option,default,thetype,allowed,valerrormsg):
        if not self.has_section(section):
            self.add_section(section)
        if self.has_option(section, option.lower()):
            raise ValueError("Default for section %s, option %s already set."%(section,option))
        self.allowed[(section,option.lower())]=(thetype,allowed)
        self.errormsg[(section,option.lower())]=valerrormsg
        self.set(section,option.lower(),default)
        
    def sethidden(self,section, option, default, thetype):
        self.setdefault(section, option, default, thetype, lambda x:True, "invalid value")

    def check(self, section, option, value):
        try:
            allowed = self.allowed[(section,option.lower())]
        except:
            raise ValueError("No default option for section %s, option %s set."%(section,option))
        try:
            allowed[0](value)
        except ValueError as e:
            print e
            raise ValueError("Option %s in section %s must be %s" % (option, section, str(allowed[0])))
        
        if not allowed[1](value):
            raise ValueError("Option %s in section %s %s, not" % (option, section, self.errormsg[(section, option)]),value)

    def set(self,section,option,value):
        self.check(section, option.lower(), value)
        super(SecureConfigParser,self).set(section,option.lower(),str(value))
    
    def get(self,section,option):
        if not self.allowed.has_key((section,option.lower())):
            raise ValueError("Invalid config option requested.")
        return self.allowed[(section,option.lower())][0](super(SecureConfigParser,self).get(section,option.lower()))
        
    def checkall(self):
        for section in self.sections():
            for option in self.options(section):
                print section,option
                self.check(section,option.lower(),self.get(section, option.lower()))
                
    def read(self,fname):
        super(SecureConfigParser,self).read(fname)
        self.checkall()
            
    def readfp(self,fp):
        super(SecureConfigParser,self).readfp(fp)
        self.checkall()
    
    def makeReadonly(self):
        for key in self.allowed.keys():
            tp = self.allowed[key][0]
            self.allowed[key] = (tp,None)
    
        
class FRETConfigParser(SecureConfigParser):
    def __init__(self):
        SecureConfigParser.__init__(self)
        self.setdefault("Dye Constants","tauD",4000,float,lambda x: x>0,"must be positive")
        self.setdefault("Dye Constants","tauA",3900,float,lambda x: x>0,"must be positive")
        self.setdefault("Dye Constants","QD",0.9,float,lambda x: x>0,"must be positive")
        self.setdefault("Dye Constants","QA",0.8,float,lambda x: x>0,"must be positive")

        self.setdefault("FRET Constants","R0",5.4,float,lambda x: x>0,"must be positive")
        self.setdefault("FRET Constants","kappa",0.6666666,float,lambda x: x>0,"must be positive")

        self.setdefault("Burst Size Distribution", "method", "analytical", str, lambda x: x in ["analytical","file"], "must be analytical or file")
        self.setdefault("Burst Size Distribution", "llimit", 20, int, lambda x: x>0,"must be positive" )
        self.setdefault("Burst Size Distribution", "ulimit", 80, int, lambda x: x>0,"must be positive" )
        self.setdefault("Burst Size Distribution", "lambda", -2.3, float, lambda x: x<0,"must be negative" )
        self.setdefault("Burst Size Distribution", "apply", "true-photon", str, lambda x: x in ["true-photon","corrected"] , "must be true-photon or corrected")
        
        self.setdefault("Burst Accumulation", "method", "trajectory", str, lambda x: x in ["trajectory","same-species","all"] ,"must be one of trajectory, same-species or all" )
        
        self.setdefault("Monte Carlo", "minstarttraj", 0, int, lambda x: x>=0,"must be positive or zero" )
        self.setdefault("Monte Carlo", "maxstarttraj", 1000, int, lambda x: x>=0,"must be positive or zero" )
        self.setdefault("Monte Carlo", "deltat", 10., float, lambda x: x>0,"must be positive" )
        self.setdefault("Monte Carlo", "nbursts", 50, int, lambda x: x>0,"must be positive" )
        self.setdefault("Monte Carlo", "photrejectdist", 0.0, float, lambda x: x>=0,"must be positive or zero" )
        self.setdefault("Monte Carlo", "rejectretry", 10, int, lambda x: x>0,"must be positive" )
        self.setdefault("Monte Carlo", "globalrejectretry",100000,int, lambda x:x>=0, "must be larger or equal zero")
        
        self.setdefault("System", "photongenerator", "cextension", str, lambda x: x in ["cython","python","cextension"] ,"must be one of cython, python or cextension" )
        self.setdefault("System", "ncpu", -1, int , lambda x: x>0 or x==-1 ,"must be larger than zero or -1 for auto-detect" )
        self.setdefault("System", "verbose", 1, int , lambda x: x==0 or x==1 ,"must be 0 or 1" )
        self.setdefault("System", "blocksize", 100, int , lambda x: x>0 ,"must be larger than zero" )

        self.setdefault('Photon Flooding','photoncount',10, int,lambda x: x>0,"must be positive"  )
        self.setdefault('Photon Flooding','startclip',0, int, lambda x: x>=0,"must be positive or zero"  )
        self.setdefault('Photon Flooding','endclip',0, int, lambda x: x>=0,"must be positive or zero"  )
        
class ReconstructionConfigParser(SecureConfigParser):
    def __init__(self):
        SecureConfigParser.__init__(self)
        self.setdefault("Transfer Matrix", "type", "global", str, lambda x:x in ("global","local","none"), "must be one of global local or none")
        self.setdefault("Transfer Matrix", "R0", 5.4, float, lambda x: x>0,"must be positive" )
        self.setdefault("Transfer Matrix", "distance bins", 80, int,lambda x: x>0,"must be positive" )
        self.setdefault("Transfer Matrix", "efficiency bins", 50, int, lambda x: x>0,"must be positive" )
        self.setdefault("Transfer Matrix", "bursts per bin", 800, int, lambda x: x>0,"must be positive" )
        self.setdefault("Transfer Matrix", "from distance", -1, float, lambda x: x!=0 , "must be non-zero")
        self.setdefault("Transfer Matrix", "to distance", -1, float, lambda x: x!=0 , "must be non-zero")
        
        self.setdefault("Burst Size Distribution", "method", "analytical", str, lambda x: x in ["analytical","file"], "must be analytical or file")
        self.setdefault("Burst Size Distribution", "llimit", 20, int, lambda x: x>0,"must be positive" )
        self.setdefault("Burst Size Distribution", "ulimit", 80, int, lambda x: x>0,"must be positive" )
        self.setdefault("Burst Size Distribution", "lambda", -2.3, float, lambda x: x<0,"must be negative" )
        
        self.setdefault("Reverse Model Fit","maxruntime",100,int,lambda x: x>0, "must be positive")
        self.setdefault("Reverse Model Fit","nr gaussian",2,int,lambda x: x>0, "must be positive")
        self.setdefault("Reverse Model Fit","sigma min",0.1,float,lambda x: x>0, "must be positive")
        self.setdefault("Reverse Model Fit","sigma max",4,float,lambda x: x>0 , "must be positive")
        self.setdefault("Reverse Model Fit","x0 min",-1,float,lambda x: x!=0 , "must be non-zero")
        self.setdefault("Reverse Model Fit","x0 max",-1,float,lambda x: x!=0 , "must be non-zero")
        self.setdefault("Reverse Model Fit","prefact min",0.0,float,lambda x: x>=0, "must be positive or zero")
        self.setdefault("Reverse Model Fit","prefact max",1.0,float,lambda x: x>0, "must be positive")
        self.setdefault("Reverse Model Fit","penaltyfact",0.0,float,lambda x: x>=0, "must be zero or positive")
     
        