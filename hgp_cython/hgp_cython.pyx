import math
import numpy as np # is the duplicate import really necessary?
cimport numpy as np
DTYPE = np.double

cdef double product_center_1D(double alpha, double ax, double beta, double bx):
    return (alpha*ax+beta*bx)/(alpha+beta)

cdef double dist2(double xa,double ya,double za,double xb,double yb,double zb):
    return math.pow(xa-xb,2)+math.pow(ya-yb,2)+math.pow(za-zb,2)

#cdef double dist(double xa,double ya,double za,double xb,double yb,double zb):
#    return sqrt(dist2(xa,ya,za,xb,yb,zb))

def vrr(xyza,double norma,lmna,double alphaa,
        xyzb,double normb,double alphab,
        xyzc,double normc,lmnc,double alphac,
        xyzd,double normd,double alphad,int M):

    from PyQuante.cints import Fgamma

    cdef double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd,px,py,pz,qx,qy,qz,zeta,eta
    cdef double wx,wy,wz,rab2,rcd2,rpq2,Kab,Kcd,T
    cdef int la,ma,na,lc,mc,nc,mtot,im,i,j,k,q,r,s

    xa = xyza[0]
    ya = xyza[1]
    za = xyza[2]
    xb = xyzb[0]
    yb = xyzb[1]
    zb = xyzb[2]
    xc = xyzc[0]
    yc = xyzc[1]
    zc = xyzc[2]
    xd = xyzd[0]
    yd = xyzd[1]
    zd = xyzd[2]
    la = lmna[0]
    ma = lmna[1]
    na = lmna[2]
    lc = lmnc[0]
    mc = lmnc[1]
    nc = lmnc[2]

    px = product_center_1D(alphaa,xa,alphab,xb)
    py = product_center_1D(alphaa,ya,alphab,yb)
    pz = product_center_1D(alphaa,za,alphab,zb)
    qx = product_center_1D(alphac,xc,alphad,xd)
    qy = product_center_1D(alphac,yc,alphad,yd)
    qz = product_center_1D(alphac,zc,alphad,zd)
    zeta = alphaa+alphab
    eta = alphac+alphad
    wx = product_center_1D(zeta,px,eta,qx)
    wy = product_center_1D(zeta,py,eta,qy)
    wz = product_center_1D(zeta,pz,eta,qz)

    rab2 = dist2(xa,ya,za,xb,yb,zb)
    Kab = np.sqrt(2)*math.pow(np.pi,1.25)/(alphaa+alphab)\
          *np.exp(-alphaa*alphab/(alphaa+alphab)*rab2)
    rcd2 = dist2(xc,yc,zc,xd,yd,zd)
    Kcd = np.sqrt(2)*math.pow(np.pi,1.25)/(alphac+alphad)\
          *np.exp(-alphac*alphad/(alphac+alphad)*rcd2)
    rpq2 = dist2(px,py,pz,qx,qy,qz)
    T = zeta*eta/(zeta+eta)*rpq2

    mtot = la+ma+na+lc+mc+nc+M

    Fgterms = np.ndarray(mtot+1,DTYPE)
    Fgterms[mtot] = Fgamma(mtot,T)
    for im in xrange(mtot-1,-1,-1):
        Fgterms[im]=(2.*T*Fgterms[im+1]+np.exp(-T))/(2.*im+1)

    # Store the vrr values as a 7 dimensional array
    vrr_terms = np.ndarray([la,ma,na,lc,mc,nc,mtot+1],DTYPE)
    #vrr_terms = {}
    for im in xrange(mtot+1):
        vrr_terms[0,0,0,0,0,0,im] = (
            norma*normb*normc*normd*Kab*Kcd/np.sqrt(zeta+eta)*Fgterms[im]
            )

    for i in xrange(la):
        for im in xrange(mtot-i):
            vrr_terms[i+1,0,0, 0,0,0, im] = (
                (px-xa)*vrr_terms[i,0,0, 0,0,0, im]
                + (wx-px)*vrr_terms[i,0,0, 0,0,0, im+1]
                )
            if i:
                vrr_terms[i+1,0,0, 0,0,0, im] += (
                    i/2./zeta*( vrr_terms[i-1,0,0, 0,0,0, im]
                               - eta/(zeta+eta)*vrr_terms[i-1,0,0, 0,0,0, im+1]
                               ))

    for j in xrange(ma):
        for i in xrange(la+1):
            for im in xrange(mtot-i-j):
                vrr_terms[i,j+1,0, 0,0,0, im] = (
                    (py-ya)*vrr_terms[i,j,0, 0,0,0, im]
                    + (wy-py)*vrr_terms[i,j,0, 0,0,0, im+1]
                    )
                if j:
                    vrr_terms[i,j+1,0, 0,0,0, im] += (
                        j/2./zeta*(vrr_terms[i,j-1,0, 0,0,0, im]
                                  - eta/(zeta+eta)
                                  *vrr_terms[i,j-1,0, 0,0,0, im+1]
                                  ))


    for k in xrange(na):
        for j in xrange(ma+1):
            for i in xrange(la+1):
                for im in xrange(mtot-i-j-k):
                    vrr_terms[i,j,k+1, 0,0,0, im] = (
                        (pz-za)*vrr_terms[i,j,k, 0,0,0, im]
                        + (wz-pz)*vrr_terms[i,j,k, 0,0,0, im+1]
                        )
                    if k:
                        vrr_terms[i,j,k+1, 0,0,0, im] += (
                            k/2./zeta*(vrr_terms[i,j,k-1, 0,0,0, im]
                                      - eta/(zeta+eta)
                                      *vrr_terms[i,j,k-1, 0,0,0, im+1]
                                      ))

    for q in xrange(lc):
        for k in xrange(na+1):
            for j in xrange(ma+1):
                for i in xrange(la+1):
                    for im in xrange(mtot-i-j-k-q):
                        vrr_terms[i,j,k, q+1,0,0, im] = (
                            (qx-xc)*vrr_terms[i,j,k, q,0,0, im]
                            + (wx-qx)*vrr_terms[i,j,k, q,0,0, im+1]
                            )
                        if q:
                            vrr_terms[i,j,k, q+1,0,0, im] += (
                                q/2./eta*(vrr_terms[i,j,k, q-1,0,0, im]
                                         - zeta/(zeta+eta)
                                         *vrr_terms[i,j,k, q-1,0,0, im+1]
                                         ))
                        if i:
                            vrr_terms[i,j,k, q+1,0,0, im] += (
                                i/2./(zeta+eta)*vrr_terms[i-1,j,k, q,0,0, im+1]
                                )

    for r in xrange(mc):
        for q in xrange(lc+1):
            for k in xrange(na+1):
                for j in xrange(ma+1):
                    for i in xrange(la+1):
                        for im in xrange(mtot-i-j-k-q-r):
                            vrr_terms[i,j,k, q,r+1,0, im] = (
                                (qy-yc)*vrr_terms[i,j,k, q,r,0, im]
                                + (wy-qy)*vrr_terms[i,j,k, q,r,0, im+1]
                                )
                            if r:
                                vrr_terms[i,j,k, q,r+1,0, im] += (
                                    r/2./eta*(vrr_terms[i,j,k, q,r-1,0, im]
                                             - zeta/(zeta+eta)
                                             *vrr_terms[i,j,k, q,r-1,0, im+1]
                                             ))
                            if j:
                                vrr_terms[i,j,k, q,r+1,0, im] += (
                                    j/2./(zeta+eta)*vrr_terms[i,j-1,k,q,r,0,im+1]
                                    )

    for s in xrange(nc):
        for r in xrange(mc+1):
            for q in xrange(lc+1):
                for k in xrange(na+1):
                    for j in xrange(ma+1):
                        for i in xrange(la+1):
                            for im in xrange(mtot-i-j-k-q-r-s):
                                vrr_terms[i,j,k,q,r,s+1,im] = (
                                    (qz-zc)*vrr_terms[i,j,k,q,r,s,im]
                                    + (wz-qz)*vrr_terms[i,j,k,q,r,s,im+1]
                                    )
                                if s:
                                    vrr_terms[i,j,k,q,r,s+1,im] += (
                                        s/2./eta*(vrr_terms[i,j,k,q,r,s-1,im]
                                                 - zeta/(zeta+eta)
                                                 *vrr_terms[i,j,k,q,r,s-1,im+1]
                                                 ))
                                if k:
                                    vrr_terms[i,j,k,q,r,s+1,im] += (
                                        k/2./(zeta+eta)*vrr_terms[i,j,k-1,q,r,s,im+1]
                                        )
    return vrr_terms[la,ma,na,lc,mc,nc,M]
