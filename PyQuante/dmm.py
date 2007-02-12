#!/usr/bin/env python
"""\
 dmm.py - Implementation of density matrix methods, including 
  Niklasson, Tymczak and Challacombe's density matrix minimizer
  JCP 118, 8611 (2003)

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

# Status:
#  PM,TCP,TRS work for h2/* and h2o/*. Appears noisy for
#     h2o/631gss. Also, when lots of iterations are done in the
#     dm convergence, the program gives unreliable results. This
#     normally kicks in around 70 or 80 iterations.
#  MCW works, provided we have a guess for efermi

from PyQuante.Ints import getbasis,getints,get2JmK
from PyQuante.Molecule import Molecule
from PyQuante.LA2 import mkdens,SymOrth,simx
from PyQuante.hartree_fock import get_energy
from NumWrap import diagonal,matrixmultiply,identity

def gershgorin_minmax(A):
    n,m = A.shape
    mins = []
    maxs = []
    for i in range(n):
        offsum = sum(abs(A[i,:]))-abs(A[i,i])
        mins.append(A[i,i]-offsum)
        maxs.append(A[i,i]+offsum)
    return min(mins),max(maxs)
        
def trace(A): return sum(diagonal(A))

def Dinit_mcw(F,Ne,tol=1e-7,maxit=100):
    #Solve for efermi and D0 using bisection:
    beta = 0.5
    emin,emax = gershgorin_minmax(F)
    I = identity(F.shape[0],'d')
    elow = emin
    ehigh = emax+20
    de = emax-elow
    alpha = beta/de
    nelow = trace(alpha*(elow*I-F) + beta*I)
    nehigh = trace(alpha*(ehigh*I-F) + beta*I)

    for i in range(100):
        efermi = 0.5*(elow+ehigh)
        nefermi = trace(alpha*(efermi*I-F)+ beta*I)
        print elow,ehigh,nelow,Ne,nehigh
        if abs(Ne-nefermi) < tol: break
        if nefermi < Ne:
            elow = efermi
            nelow = nefermi
        elif nefermi > Ne:
            ehigh = efermi
            nehigh = nefermi
    alpha = min(beta/(emax-efermi),(1-beta)/(efermi-emin))
    return alpha*(efermi*I-F)+beta*I

def DMP(F,S,Ne,Method=0,MaxIter=50,ErrorLimit=1e-12):
    # Density Matrix Purification Methods
    # 0 -> Trace correcting purification (default)
    # 1 -> Trace resetting
    # 2 -> McWeeny purification
    # 3 -> Canonical purification
    methods = ['TCP','TRS','MCW','PM']

    # Step 1: Orthogonalize the Fock matrix:
    X = SymOrth(S)
    F = simx(F,X)

    # Step 2: Initialize the density matrix:
    emin,emax = gershgorin_minmax(F)
    N = F.shape[0]
    I = identity(N,'d')
    if Method == 0 or Method == 1:
        D = (emax*I-F)/(emax-emin)
    elif Method == 2:
        D = Dinit_mcw(F,Ne)
    elif Method == 3:
        efermi = trace(F)/N
        beta = Ne/float(N)
        alpha = min(Ne/(emax-efermi),(N-Ne)/(efermi-emin))/float(N)
        D = alpha*(efermi*I-F) + beta*I
    else:
        raise "Unknown method %d" % Method

    Dsumold = sum(sum(D))
    # Step 3: Iterate on DM updates:
    for iter in range(MaxIter):
        Ne_curr = trace(D)
        D2 = matrixmultiply(D,D)
        if Method == 0:
            if Ne_curr < Ne:
                D = 2.0*D-D2
            else:
                D = D2
            if abs(Ne_curr-Ne) < ErrorLimit: break
        elif Method == 1:
            Df = matrixmultiply(D2,4*D-3*D2)
            trf = trace(Df)
            Dp = I-D
            Dp2 = matrixmultiply(Dp,Dp)
            Dg = matrixmultiply(D2,Dp2)
            trg = trace(Dg)
            gamma = (Ne-trf)/trg
            if gamma > 2:
                D = 2*D-D2
            elif gamma < 0:
                D = D2
            else:
                D = Df-gamma*Dg
            if abs(Ne_curr-Ne) < ErrorLimit: break
        elif Method == 2:
            D = 3*D2-2*matrixmultiply(D,D2)
            if abs(Ne_curr-Ne) < ErrorLimit: break
        elif Method == 3:
            D3 = matrixmultiply(D,D2)
            cn = trace(D2-D3)/trace(D-D2)
            if cn < 0.5:
                D = ((1.0-2.0*cn)*D+(1.0+cn)*D2-D3)/(1.0-cn)
            else:
                D = ((1+cn)*D2-D3)/cn
            Dsum = sum(sum(D))
            if Dsum-Dsumold  < ErrorLimit: break
            Dsumold = Dsum
    else: print "DMM: Warning MaxIters reached"
    #print "%s converged in %d iters" % (methods[Method],iter)
    D = simx(D,X,'T')
    return D
        
    
