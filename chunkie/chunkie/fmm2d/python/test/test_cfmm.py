#!/usr/bin/env python

import fmm2dpy as fmm
import numpy as np
import numpy.linalg as la

def main():
    test_cfmm()

def test_cfmm():
    ntests = 36
    testres = np.zeros(ntests)
    #
    #  This is a testing code for making sure all the 
    #  fmm routines are accessible through fmm2d.py
    #

    n = 2000
    ntest = 100
    zk = 1.1 + 1j*0
    sources = np.random.uniform(0,1,(2,n))
    stmp = sources[:,0:ntest]

    nt = 1880
    targ = np.random.uniform(0,1,(2,nt))
    ttmp = targ[:,0:ntest]
    eps = 10**(-5)

    charges = np.random.uniform(0,1,n) + 0j
    dipstr = np.random.uniform(0,1,n)+1j*np.random.uniform(0,1,n)

    outex=fmm.Output()

    itest = 0
    out=fmm.cfmm2d(eps=eps,sources=sources,charges=charges,pg=1)
    out2 = fmm.c2ddir(sources=sources,targets=stmp,charges=charges,pgt=1)
    out2.pot = out2.pottarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=1,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges, pot")

    itest = itest+1

    out=fmm.cfmm2d(eps=eps,sources=sources,dipstr=dipstr,pg=1)
    out2 = fmm.c2ddir(sources=sources,targets=stmp,dipstr=dipstr,pgt=1)
    out2.pot = out2.pottarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=1,cauchy=1)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, dipoles, pot")


    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,charges=charges, \
        dipstr=dipstr,pg=1)
    out2 = fmm.c2ddir(sources=sources,targets=stmp,charges=charges, \
        dipstr=dipstr,pgt=1)
    out2.pot = out2.pottarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=1,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges and dipoles, pot")


    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,charges=charges,pg=2)
    out2 = fmm.c2ddir(sources=sources,targets=stmp,charges=charges, \
        pgt=2)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=2,cauchy=1)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges, pot and grad")
        

    itest = itest+1

    out=fmm.cfmm2d(eps=eps,sources=sources,dipstr=dipstr,pg=2)
    out2 = fmm.c2ddir(sources=sources,targets=stmp,dipstr=dipstr, \
        pgt=2)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=2,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, dipoles, pot and grad")

    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,charges=charges, \
        dipstr=dipstr,pg=2)
    out2 = fmm.c2ddir(sources=sources,targets=stmp,charges=charges,
       dipstr=dipstr,pgt=2)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=2,cauchy=1)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges and dipoles, pot and grad")

    itest=itest+1
    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=1)
    out2=fmm.c2ddir(sources=sources,targets=ttmp,charges=charges,pgt=1)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=1,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges, pot")

    itest = itest+1

    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,\
        dipstr=dipstr,pgt=1)
      
    out2=fmm.c2ddir(sources=sources,targets=ttmp,\
        dipstr=dipstr,pgt=1)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=1,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, dipoles, pot")

    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ, \
        charges=charges, \
        dipstr=dipstr,pgt=1)

    out2=fmm.c2ddir(sources=sources,targets=ttmp, \
        charges=charges, \
        dipstr=dipstr,pgt=1)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=1,cauchy=1)
    
    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges and dipoles, pot")

    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=2)
    out2=fmm.c2ddir(sources=sources,targets=ttmp,charges=charges,pgt=2)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=2,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges, pot and grad")

    itest = itest+1

    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,\
    dipstr=dipstr,pgt=2)
    out2=fmm.c2ddir(sources=sources,targets=ttmp,\
    dipstr=dipstr,pgt=2)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=2,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, dipoles, pot and grad")

    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,charges=charges,\
        dipstr=dipstr,pgt=2)
    out2 =fmm.c2ddir(sources=sources,targets=ttmp,charges=charges,\
        dipstr=dipstr,pgt=2)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=2,cauchy=1)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges and dipoles, pot and grad")

    itest = itest+1
    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=1,pg=1)
    out2=fmm.c2ddir(sources=sources,targets=stmp,charges=charges,pgt=1)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.c2ddir(sources=sources,targets=ttmp,charges=charges,pgt=1)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=1,pgt=1,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges, pot")
    itest = itest+1

    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,\
        dipstr=dipstr,pgt=1,pg=1)
    out2=fmm.c2ddir(sources=sources,targets=stmp, \
        dipstr=dipstr,pgt=1)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.c2ddir(sources=sources,targets=ttmp, \
          dipstr=dipstr,pgt=1)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=1,pgt=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, dipoles, pot")

    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ, \
        charges=charges, \
        dipstr=dipstr,pgt=1,pg=1)
    out2=fmm.c2ddir(sources=sources,targets=stmp,charges=charges, \
        dipstr=dipstr,pgt=1)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.c2ddir(sources=sources,targets=ttmp,charges=charges, \
          dipstr=dipstr,pgt=1)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=1,pgt=1,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges and dipoles, pot")

    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=2,pg=2)
    out2=fmm.c2ddir(sources=sources,targets=stmp,charges=charges, \
        pgt=2)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.c2ddir(sources=sources,targets=ttmp,charges=charges, \
          pgt=2)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=2,pgt=2,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges, pot and grad")
    itest = itest+1

    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,\
    dipstr=dipstr,pgt=2,pg=2)
    out2=fmm.c2ddir(sources=sources,targets=stmp, \
        dipstr=dipstr,pgt=2)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.c2ddir(sources=sources,targets=ttmp,dipstr=dipstr, \
          pgt=2)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=2,pgt=2,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, dipoles, pot and grad")

    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,charges=charges,\
        dipstr=dipstr,pgt=2,pg=2)
    out2=fmm.c2ddir(sources=sources,targets=stmp,charges=charges, \
        dipstr=dipstr,pgt=2)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.c2ddir(sources=sources,targets=ttmp,charges=charges, \
          dipstr=dipstr,pgt=2)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=2,pgt=2,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges and dipoles, pot and grad")



    nd = 2
    charges = np.random.uniform(0,1,(nd,n))
    dipstr = np.random.uniform(0,1,(nd,n))

    itest = itest+1
    out=fmm.cfmm2d(eps=eps,sources=sources,charges=charges,pg=1,nd=nd)
    out2 = fmm.c2ddir(sources=sources,targets=stmp,charges=charges,pgt=1,nd=nd)
    out2.pot = out2.pottarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=1,nd=nd,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges, pot, vectorized")

    itest = itest+1

    out=fmm.cfmm2d(eps=eps,sources=sources,dipstr=dipstr,pg=1,nd=nd)
    out2 = fmm.c2ddir(sources=sources,targets=stmp,dipstr=dipstr,pgt=1,nd=nd)
    out2.pot = out2.pottarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=1,nd=nd,cauchy=1)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, dipoles, pot, vectorized")


    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,charges=charges, \
        dipstr=dipstr,pg=1,nd=nd)
    out2 = fmm.c2ddir(sources=sources,targets=stmp,charges=charges, \
        dipstr=dipstr,pgt=1,nd=nd)
    out2.pot = out2.pottarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=1,nd=nd,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges and dipoles, pot, vectorized")


    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,charges=charges,pg=2,nd=nd)
    out2 = fmm.c2ddir(sources=sources,targets=stmp,charges=charges, \
        pgt=2,nd=nd)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=2,nd=nd,cauchy=1)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges, pot and grad, vectorized")
        

    itest = itest+1

    out=fmm.cfmm2d(eps=eps,sources=sources,dipstr=dipstr,pg=2,nd=nd)
    out2 = fmm.c2ddir(sources=sources,targets=stmp,dipstr=dipstr, \
        pgt=2,nd=nd)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=2,nd=nd,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, dipoles, pot and grad, vectorized")

    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,charges=charges, \
        dipstr=dipstr,pg=2,nd=nd)
    out2 = fmm.c2ddir(sources=sources,targets=stmp,charges=charges,
       dipstr=dipstr,pgt=2,nd=nd)
    out2.pot = out2.pottarg
    out2.grad = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pg=2,nd=nd,cauchy=1)


    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources, charges and dipoles, pot and grad, vectorized")
    

    itest=itest+1
    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=1,nd=nd)
    out2=fmm.c2ddir(sources=sources,targets=ttmp,charges=charges,pgt=1,nd=nd)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=1,nd=nd,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges, pot, vectorized")

    itest = itest+1

    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,\
        dipstr=dipstr,pgt=1,nd=nd)
      
    out2=fmm.c2ddir(sources=sources,targets=ttmp,\
        dipstr=dipstr,pgt=1,nd=nd)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=1,nd=nd,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, dipoles, pot, vectorized")

    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ, \
        charges=charges, \
        dipstr=dipstr,pgt=1,nd=nd)

    out2=fmm.c2ddir(sources=sources,targets=ttmp, \
        charges=charges, \
        dipstr=dipstr,pgt=1,nd=nd)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=1,nd=nd,cauchy=1)
    
    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges and dipoles, pot, vectorized")

    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=2,nd=nd)
    out2=fmm.c2ddir(sources=sources,targets=ttmp,charges=charges,pgt=2,nd=nd)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=2,nd=nd,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges, pot and grad, vectorized")

    itest = itest+1

    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,\
    dipstr=dipstr,pgt=2,nd=nd)
    out2=fmm.c2ddir(sources=sources,targets=ttmp,\
    dipstr=dipstr,pgt=2,nd=nd)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=2,nd=nd,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, dipoles, pot and grad, vectorized")

    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,charges=charges,\
        dipstr=dipstr,pgt=2,nd=nd)
    out2 =fmm.c2ddir(sources=sources,targets=ttmp,charges=charges,\
        dipstr=dipstr,pgt=2,nd=nd)
    err = fmm.comperr(ntest=ntest,out=out,outex=out2,pgt=2,nd=nd,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to targets, charges and dipoles, pot and grad, vectorized")

    itest = itest+1
    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=1,pg=1,nd=nd)
    out2=fmm.c2ddir(sources=sources,targets=stmp,charges=charges,pgt=1,nd=nd)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.c2ddir(sources=sources,targets=ttmp,charges=charges,pgt=1,nd=nd)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=1,pgt=1,nd=nd,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges, pot, vectorized")
    itest = itest+1

    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,\
        dipstr=dipstr,pgt=1,pg=1,nd=nd)
    out2=fmm.c2ddir(sources=sources,targets=stmp, \
        dipstr=dipstr,pgt=1,nd=nd)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.c2ddir(sources=sources,targets=ttmp, \
          dipstr=dipstr,pgt=1,nd=nd)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=1,pgt=1,nd=nd,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, dipoles, pot, vectorized")

    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ, \
        charges=charges, \
        dipstr=dipstr,pgt=1,pg=1,nd=nd)
    out2=fmm.c2ddir(sources=sources,targets=stmp,charges=charges, \
        dipstr=dipstr,pgt=1,nd=nd)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.c2ddir(sources=sources,targets=ttmp,charges=charges, \
          dipstr=dipstr,pgt=1,nd=nd)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=1,pgt=1,nd=nd,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges and dipoles, pot, vectorized")

    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,charges=charges,pgt=2,pg=2,nd=nd)
    out2=fmm.c2ddir(sources=sources,targets=stmp,charges=charges, \
        pgt=2,nd=nd)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.c2ddir(sources=sources,targets=ttmp,charges=charges, \
          pgt=2,nd=nd)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=2,pgt=2,nd=nd,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges, pot and grad, vectorized")
    itest = itest+1

    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,\
    dipstr=dipstr,pgt=2,pg=2,nd=nd)
    out2=fmm.c2ddir(sources=sources,targets=stmp, \
        dipstr=dipstr,pgt=2,nd=nd)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.c2ddir(sources=sources,targets=ttmp,dipstr=dipstr, \
          pgt=2,nd=nd)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=2,pgt=2,nd=nd,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, dipoles, pot and grad, vectorized")

    itest = itest + 1
    out=fmm.cfmm2d(eps=eps,sources=sources,targets=targ,charges=charges,\
        dipstr=dipstr,pgt=2,pg=2,nd=nd)
    out2=fmm.c2ddir(sources=sources,targets=stmp,charges=charges, \
        dipstr=dipstr,pgt=2,nd=nd)
    outex.pot = out2.pottarg
    outex.grad = out2.gradtarg
    out2=fmm.c2ddir(sources=sources,targets=ttmp,charges=charges, \
          dipstr=dipstr,pgt=2,nd=nd)
    outex.pottarg = out2.pottarg
    outex.gradtarg = out2.gradtarg
    err = fmm.comperr(ntest=ntest,out=out,outex=outex,pg=2,pgt=2,nd=nd,cauchy=1)

    if(err<eps):
        testres[itest] = 1
    else:
        print("Failed sources to sources and targets, charges and dipoles, pot and grad, vectorized")

    if(sum(testres)==ntests):
        print("all cfmm tests succeeded")

if __name__ == "__main__":
    main()
