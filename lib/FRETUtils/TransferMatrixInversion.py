'''
Created on 15.10.2012

@author: mhoefli


'''

import numpy
from openopt import GLP,NLP
import matplotlib as mpl
import matplotlib.pyplot as plt


setBackend="TkAgg"

if mpl.rcParams["backend"]!=setBackend:
    print "Switching matplotlib backend from",mpl.rcParams["backend"],"to",setBackend 
    mpl.rcParams["backend"]=setBackend
    plt.switch_backend(setBackend)

mpl.rcdefaults()  

def GaussianRegularizationDistanceReconstruction(config, TM, effhist):

    Rmin=config.get("Transfer Matrix", "from distance")
    Rmax=config.get("Transfer Matrix", "to distance")
    if config.get("Reverse Model Fit","x0 min")<0:
        config.set("Reverse Model Fit","x0 min",Rmin)
    if config.get("Reverse Model Fit","x0 max")<0:
        config.set("Reverse Model Fit","x0 max",Rmax)

    grmin=config.get("Reverse Model Fit","x0 min")
    grmax=config.get("Reverse Model Fit","x0 max")
    gsigmin=config.get("Reverse Model Fit","sigma min")
    gsigmax=config.get("Reverse Model Fit","sigma max")
    gamin=config.get("Reverse Model Fit","prefact min")
    gamax=config.get("Reverse Model Fit","prefact max")
    ngauss = config.get("Reverse Model Fit","nr gaussian")
    maxtime= config.get("Reverse Model Fit","maxruntime")
    
    lbounds = [gamin] * ngauss + [grmin] * ngauss + [gsigmin] * ngauss
    ubounds = [gamax] * ngauss + [grmax] * ngauss + [gsigmax] * ngauss
    
    r_prdist,xrange,e_fitprdist,fitvals = fittingOpenopt(effhist,TM,Rmin,Rmax,lbounds,ubounds,maxtime)
    return r_prdist,xrange,e_fitprdist,fitvals


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

def penalizeCloseGauss(argvec,TM,targeteff,xxarr):
    stddev = gaussSQDiff(argvec,TM,targeteff,xxarr)
    nrgauss = len(argvec)/3
    r_vals=argvec[nrgauss:2*nrgauss]
    sig_vals=argvec[2*nrgauss:]
    dists = numpy.subtract.outer(r_vals,r_vals.T)
    ssums = numpy.add.outer(sig_vals,sig_vals.T)
    absdist = numpy.sqrt(dists*dists)
    distsigdiff = absdist-0.5*ssums
    distsigdiffnotr = distsigdiff - distsigdiff * numpy.eye(distsigdiff.shape[0])
    diff = (distsigdiffnotr < 0).sum()/2
    return stddev * 10 ** diff

def plotCallback(p,lines_dist,lines_eff,lines_g,xxarr,TM):
#    print p.__dict__
    argvec = p.xk
    nrgauss = len(argvec)/3
    a_vals=argvec[0:nrgauss]
    r_vals=argvec[nrgauss:2*nrgauss]
    sig_vals=argvec[2*nrgauss:] 
         
    gaussians = (a_vals*numpy.exp(-(xxarr.T-r_vals)**2/(2.0*sig_vals**2)))
    r_prdist = gaussians.sum(axis=1)
#    r_prdist /= r_prdist.max()
    e_prdist=numpy.dot(r_prdist,TM.getMatrix())  
    e_prdist=e_prdist/e_prdist.mean()    
    
    for gauss in range(nrgauss):
        lines_g[gauss].set_ydata(gaussians[:,gauss])#/r_prdist.max()      
    
    lines_dist.set_ydata(r_prdist)
    lines_eff.set_ydata(e_prdist)
    
    plt.draw()
    
    return False

def createLivePlot(nrgauss,pearr,tmatrix,xarr,lbounds,ubounds):
    plt.figure(figsize=(10, 8))
    plt.ion()
    plt.subplot(221)
    g_lines = []
    for gauss in range(nrgauss):
        ln, = plt.plot(xarr, numpy.ones(len(xarr)), label="G%d" % gauss)
        g_lines.append(ln)
    
    maxgauss=numpy.array(ubounds)[nrgauss:2*nrgauss].sum()
    
    lines_distance, = plt.plot(xarr, numpy.ones(len(xarr))*maxgauss, label="sum", linewidth=2,linestyle=":")
    plt.yticks(())
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
    plt.xlabel("Reconstructed distance")
    plt.ylim(0., 1.1)
    plt.subplot(222)
    xeff = numpy.linspace(0., 1., len(pearr), endpoint=False) + 1. / len(pearr) / 2
    plt.plot(xeff, pearr, label="reference")
    lines_efficiency, = plt.plot(xeff, numpy.ones(len(pearr)), label="fit")
    plt.xlabel("FRET Efficiency")
    plt.yticks(())
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
    plt.subplot(223)
    tmatrix.plot()
    plt.show()
    return lines_distance,lines_efficiency,g_lines

def fittingOpenopt(pearr,tmatrix,minR,maxR,lbounds,ubounds,gmaxtime):    
    nrgauss = len(lbounds) / 3      
    rvecbins=tmatrix.getMatrix().shape[0]
    myrange=maxR-minR
    xarr = numpy.linspace(minR+myrange/rvecbins/2,maxR-myrange/rvecbins/2,rvecbins)
    xxarr = numpy.array([xarr] * nrgauss)

    #minfuncwrap = lambda  x: gaussSQDiff(x,tmatrix.getMatrix(),pearr,xxarr)
    minfuncwrap = lambda  x: penalizeCloseGauss(x,tmatrix.getMatrix(),pearr,xxarr)
    
    lines_distance,lines_efficiency,g_lines = createLivePlot(nrgauss,pearr,tmatrix,xarr,lbounds,ubounds)
    
    mycallback =  lambda p: plotCallback(p,lines_distance,lines_efficiency,g_lines,xxarr,tmatrix)

    print "Starting openopt ##########################"
    prob = GLP(minfuncwrap,lb=lbounds,ub=ubounds,callback=mycallback,maxFunEvals=1e15,maxNonSuccess=200,maxIter=1e5,maxTime=gmaxtime)
    result=prob.solve('de',population=1000*len(lbounds))
    #result=prob.solve('asa')
    #result=prob.solve('galileo') # not good
    #result=prob.solve('pswarm')
    #prob = GLP(minfuncwrap,lb=lbounds,ub=ubounds,callback=mycallback,maxNonSuccess=200,maxIter=1e5,maxTime=gmaxtime)
    #result=prob.solve('isres',population=100*len(lbounds))
    #prob = NLP(minfuncwrap,lb=lbounds,ub=ubounds,callback=mycallback,maxNonSuccess=200,maxIter=1e5,maxTime=gmaxtime)
    #result=prob.solve('scipy_lbfgsb')
    #result=prob.solve('scipy_tnc')
    #result=prob.solve('bobyqa')
    #result=prob.solve('ptn')
    #result=prob.solve('slmvm1')
    #result=prob.solve('slmvm2')
    #result=prob.solve('ralg')
    #result=prob.solve('scipy_cobyla') #good!!
    #result=prob.solve('mma')
    #result=prob.solve('auglag')
    #result=prob.solve('gsubg')

    
    xopt=result.xf
    print "Minimum function chisq",result.ff 
            
    a_final=xopt[0:nrgauss]
    r_final=xopt[nrgauss:2*nrgauss]
    sig_final=xopt[2*nrgauss:]

    gaussians = (a_final*numpy.exp(-(xxarr.T-r_final)**2/(2.0*sig_final**2)))
    r_prdist = gaussians.sum(axis=1)  
    e_fitprdist=numpy.dot(r_prdist,tmatrix.getMatrix())    

    r_prdist/=r_prdist.sum()
    e_fitprdist/=e_fitprdist.sum()

    return r_prdist,xarr,e_fitprdist,(a_final,r_final,sig_final)
