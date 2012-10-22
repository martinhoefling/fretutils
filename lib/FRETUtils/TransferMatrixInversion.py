'''
Created on 15.10.2012

@author: mhoefli
'''

import numpy
from openopt import GLP

def GaussianRegularizationDistanceReconstruction(options, TM, effhist):
    grmin=4.
    grmax=6.
    gsigmin=0.1
    gsigmax=2.
    gamin=0.
    gamax=1.
    ngauss = 2
    Rmin=0.
    Rmax=10.
    maxtime=72
    lbounds = [gamin] * ngauss + [grmin] * ngauss + [gsigmin] * ngauss
    ubounds = [gamax] * ngauss + [grmax] * ngauss + [gsigmax] * ngauss
    
    r_prdist,e_fitprdist,fitvals = fittingOpenopt(effhist,TM,Rmin,Rmax,lbounds,ubounds,maxtime)
    return r_prdist,e_fitprdist,fitvals


def gaussSQDiff(argvec,TM,targeteff,xxarr):
    nrgauss = len(argvec)/3
    a_vals=argvec[0:nrgauss]
    r_vals=argvec[nrgauss:2*nrgauss]
    sig_vals=argvec[2*nrgauss:]
    
    gaussians = (a_vals*numpy.exp(-(xxarr.T-r_vals)**2/(2.0*sig_vals**2)))
    r_prdist = gaussians.sum(axis=1)  
    
    e_prdist=numpy.dot(r_prdist,TM)  
    e_prdist=e_prdist/e_prdist.mean()
    
    devnew=((targeteff-e_prdist)**2).mean()
 
    return devnew


def fittingOpenopt(pearr,tmatrix,minR,maxR,lbounds,ubounds,gmaxtime):    
      
    rvecbins=tmatrix.shape[0]
    myrange=maxR-minR
    xarr = numpy.linspace(minR+myrange/rvecbins/2,maxR-myrange/rvecbins/2,rvecbins)
    xxarr = numpy.array((xarr,xarr,xarr))

    minfuncwrap = lambda  x: gaussSQDiff(x,tmatrix,pearr,xxarr)

    print "Starting openopt ##########################"
    prob = GLP(minfuncwrap,lb=lbounds,ub=ubounds,maxFunEvals=1e15,maxNonSuccess=200,maxIter=1e5,maxTime=gmaxtime)#,maxIter=1e5,maxFunEvals=1e7,maxTime=3,maxCPUTime=3)
    result=prob.solve('de',population=100*len(lbounds))
    xopt=result.xf
    print "Minimum function chisq",result.ff 
        
    nrgauss = len(xopt) / 3
          
    a_final=xopt[0:nrgauss]
    r_final=xopt[nrgauss:2*nrgauss]
    sig_final=xopt[2*nrgauss:]

    gaussians = (a_final*numpy.exp(-(xxarr.T-r_final)**2/(2.0*sig_final**2)))
    r_prdist = gaussians.sum(axis=1)  
    e_fitprdist=numpy.dot(r_prdist,tmatrix)    

    r_prdist/=r_prdist.sum()
    e_fitprdist/=e_fitprdist.sum()

    return r_prdist,e_fitprdist,(a_final,r_final,sig_final)
