"""Linear-in-doubles triple-block corrections (H33_TWO/PHH_PHH/HPP_HPP)."""
import numpy as np
E=np.einsum

class TripleCorrection:
    def __init__(self,ints,t):
        g=ints.g('o','o','v','v');self.t=t;self.g=g
        r=.5*E('ikab,jkab->ij',g,t,optimize=True);self.ro=r+r.T
        r=.5*E('ijac,ijbc->ab',g,t,optimize=True);self.rv=r+r.T
        x=.5*E('ijab,klab->ijkl',g,t,optimize=True);self.xo=x+x.transpose(2,3,0,1)
        w=E('ikac,jkbc->iajb',g,t,optimize=True);self.w=w+w.transpose(2,3,0,1)

    def action(self,x,y):
        z=-.25*E('ijkl,kla->ija',self.xo,x,optimize=True)
        r=-.5*E('iakb,kjb->ija',self.w,x,optimize=True);z+=r-r.swapaxes(0,1)
        r=.5*E('ila,jl->ija',x,self.ro,optimize=True);z+=r-r.swapaxes(0,1)
        z+=.5*E('ijb,ab->ija',x,self.rv,optimize=True)
        # X_abcd factorized through occupied pairs; no V^4 intermediate.
        temp=E('klcd,icd->ikl',self.t,y,optimize=True)
        w=.125*E('klab,ikl->iab',self.g,temp,optimize=True)
        temp=E('klcd,icd->ikl',self.g,y,optimize=True)
        w+=.125*E('klab,ikl->iab',self.t,temp,optimize=True)
        r=.5*E('iajc,jcb->iab',self.w,y,optimize=True);w+=r-r.swapaxes(1,2)
        r=-.5*E('icb,ac->iab',y,self.rv,optimize=True);w+=r-r.swapaxes(1,2)
        w-=.5*E('jab,ij->iab',y,self.ro,optimize=True)
        return z,w

    def diagonal(self,ip,ea):
        dh=[];dp=[]
        for i,j,a in ip:
            dh.append(-.5*self.xo[i,j,i,j]-.5*(self.w[i,a,i,a]+self.w[j,a,j,a])
                      +.5*(self.ro[i,i]+self.ro[j,j]+self.rv[a,a]))
        for i,a,b in ea:
            xab=np.sum(self.g[:,:,a,b]*self.t[:,:,a,b])
            dp.append(.5*xab+.5*(self.w[i,a,i,a]+self.w[i,b,i,b])
                      -.5*(self.rv[a,a]+self.rv[b,b]+self.ro[i,i]))
        return np.array(dh),np.array(dp)
