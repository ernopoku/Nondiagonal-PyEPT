"""Spin-orbital tensor equations; unique pairs carry no sqrt(2) factors."""
import numpy as np
from itertools import combinations
from .methods import method_spec

E = np.einsum

def amplitudes(ints, denominator_tol=1e-10):
    o=ints.nocc; eo=ints.energy[:o]; ev=ints.energy[o:]
    den=eo[:,None,None,None]+eo[None,:,None,None]-ev[None,None,:,None]-ev[None,None,None,:]
    if np.min(np.abs(den))<denominator_tol: raise ValueError('Near-zero MP2 denominator; perturbative reference is singular.')
    t=ints.g('o','o','v','v')/den
    # Eq. (33), 2023: first-order doubles, second-order singles.
    numerator=E('jkab,bijk->ia',t,ints.g('v','o','o','o'),optimize=True)
    numerator+=E('ijbc,bcaj->ia',t,ints.g('v','v','v','o'),optimize=True)
    d=eo[:,None]-ev[None,:]
    if np.min(np.abs(d))<denominator_tol: raise ValueError('Near-zero occupied-virtual gap.')
    return t, numerator/(2*d)


def vertices(ints,t,need_ip=True,need_ea=True):
    """Full-weight second-order vertex pieces. New metric multiplies by 1/2."""
    o,v,n=ints.nocc,ints.nvir,ints.n
    ip1=ints.g('p','v','o','o').transpose(0,2,3,1)
    ea1=ints.g('p','o','v','v')
    ip2=np.zeros_like(ip1); ea2=np.zeros_like(ea1)
    for p in range(n):
        if need_ip:
            ip2[p]=.5*E('abc,ijbc->ija',ints.g([p],'v','v','v')[0],t,optimize=True)
            ring=E('kib,jkab->ija',ints.g([p],'o','o','v')[0],t,optimize=True)
            ip2[p]+=ring-ring.swapaxes(0,1)
        if need_ea:
            ea2[p]=.5*E('ikl,klab->iab',ints.g([p],'o','o','o')[0],t,optimize=True)
            ring=E('cak,ikbc->iab',ints.g([p],'v','v','o')[0],t,optimize=True)
            ea2[p]+=ring-ring.swapaxes(1,2)
    return ip1,ip2,ea1,ea2


def static_self_energy(ints,t,singles,kind):
    n,o=ints.n,ints.nocc
    rho=np.zeros((n,n))
    if kind=='none':return rho
    factor=.5 if kind=='linear' else 1.
    rho[:o,o:]=factor*singles;rho[o:,:o]=factor*singles.T
    if kind=='quadratic':
        rho[:o,:o]=-.5*E('ikab,jkab->ij',t,t,optimize=True)
        rho[o:,o:]=.5*E('ijac,ijbc->ab',t,t,optimize=True)
    sigma=np.empty((n,n))
    for p in range(n):
        sigma[p]=E('rqs,rs->q',ints.g([p],'p','p','p')[0],rho,optimize=True)
    if np.max(np.abs(sigma-sigma.T))>1e-9: raise ArithmeticError('Non-Hermitian static self-energy.')
    return sigma

class Hamiltonian:
    """Matrix-free Hermitian super-operator in one conserved spin sector.

    Layout: simple alpha/beta creators; i<j,a (2hp); i,a<b (2ph).
    Triple tensors are unpacked antisymmetrically for contractions.
    """
    def __init__(self,ints,method='NRL3',sector='ip',spin=0,*,doubles=None,
                 static_tol=1e-10,static_max_cycle=500):
        if spin not in (0,1):raise ValueError('spin must be 0 (alpha) or 1 (beta).')
        self.ints=ints;self.spec=method_spec(method,sector);self.sector=sector;self.spin=spin
        if self.spec.name == 'nD-NRL3':
            raise ValueError('Use EPT or StaticNRL3 to construct nD-NRL3.')
        o,v,n=ints.nocc,ints.nvir,ints.n
        self.simple=np.flatnonzero(ints.spin==spin)
        # Spin values +1 alpha, -1 beta. Operator spin change must match creator.
        sz=1-2*ints.spin;target=1-2*spin
        self.ip=np.array([(i,j,a) for i,j in combinations(range(o),2) for a in range(v)
                          if sz[i]+sz[j]-sz[o+a]==target],dtype=int).reshape(-1,3)
        self.ea=np.array([(i,a,b) for i in range(o) for a,b in combinations(range(v),2)
                          if sz[o+a]+sz[o+b]-sz[i]==target],dtype=int).reshape(-1,3)
        self.ns=len(self.simple);self.nh=len(self.ip);self.np=len(self.ea)
        self.shape=(self.ns+self.nh+self.np,)*2
        self.dtype=np.dtype(float)
        m=self.spec
        need_ip=bool(m.ip_vertex_h or m.ip_vertex_p)
        need_ea=bool(m.ea_vertex_h or m.ea_vertex_p)
        if m.triple_correction:
            if doubles is None:raise ValueError('BD-T1 requires converged Brueckner doubles, not MP2 amplitudes.')
            t=np.asarray(doubles);s=np.zeros((o,v))
            if t.shape!=(o,o,v,v):raise ValueError('Wrong doubles shape.')
        else:
            t,s=amplitudes(ints) if need_ip or need_ea or m.static!='none' else (None,None)
        b,c,d,f=vertices(ints,t,need_ip,need_ea)
        wi=np.r_[np.full(o,m.ip_vertex_h),np.full(v,m.ip_vertex_p)]
        we=np.r_[np.full(o,m.ea_vertex_h),np.full(v,m.ea_vertex_p)]
        b+=wi[:,None,None,None]*c;d+=we[:,None,None,None]*f
        self.bh=b[(self.simple[:,None],)+tuple(self.ip.T)] if self.nh else np.zeros((self.ns,0))
        self.bp=d[(self.simple[:,None],)+tuple(self.ea.T)] if self.np else np.zeros((self.ns,0))
        sigma=static_self_energy(ints,t,s,'none' if m.static=='dem' else m.static)
        self.a=np.diag(ints.energy[self.simple])+sigma[np.ix_(self.simple,self.simple)]
        eo=ints.energy[:o];ev=ints.energy[o:]
        if m.triple_correction:
            from .renormalization import TripleCorrection
            self.correction=TripleCorrection(ints,t)
            self.a=ints.fock[np.ix_(self.simple,self.simple)].copy()
        else:self.correction=None
        self.dh=eo[:,None,None]+eo[None,:,None]-ev[None,None,:]
        self.dp=ev[None,:,None]+ev[None,None,:]-eo[:,None,None]
        self.goooo=ints.g('o','o','o','o') if m.ip_interaction else None
        self.gvoov=ints.g('v','o','v','o') if m.ip_interaction else None
        self.govov=ints.g('o','v','o','v') if m.ea_interaction else None
        self._diag=None
        if m.static=='dem':
            from .dyson_adc import dem_static
            sigma,self.static_diagnostics=dem_static(self,tol=static_tol,
                                                     max_cycle=static_max_cycle)
            self.a+=sigma

    def unpack(self,h,p):
        o,v=self.ints.nocc,self.ints.nvir
        x=np.zeros((o,o,v));y=np.zeros((o,v,v))
        if self.nh:
            i,j,a=self.ip.T;x[i,j,a]=h;x[j,i,a]=-h
        if self.np:
            i,a,b=self.ea.T;y[i,a,b]=p;y[i,b,a]=-p
        return x,y

    def triple_action(self,h,p):
        x,y=self.unpack(h,p);z=self.dh*x;w=self.dp*y
        if self.spec.ip_interaction:
            z-=.5*E('ijkl,kla->ija',self.goooo,x,optimize=True)
            ring=E('biak,kjb->ija',self.gvoov,x,optimize=True)
            z+=ring-ring.swapaxes(0,1)
        if self.spec.ea_interaction:
            # Contract in spatial spin blocks; never allocate spin V^4 ERIs.
            w+=self.ints.virtual_ladder(y)
            ring=E('jaic,jcb->iab',self.govov,y,optimize=True)
            w-=ring-ring.swapaxes(1,2)
        if self.correction is not None:
            zh,wp=self.correction.action(x,y);z+=zh;w+=wp
        return z[tuple(self.ip.T)],w[tuple(self.ea.T)]

    def matvec(self,vec):
        vec=np.asarray(vec).reshape(-1)
        if vec.size!=self.shape[0]:raise ValueError('Wrong vector dimension.')
        s=vec[:self.ns];h=vec[self.ns:self.ns+self.nh];p=vec[self.ns+self.nh:]
        zh,zp=self.triple_action(h,p)
        return np.r_[self.a@s+self.bh@h+self.bp@p, self.bh.T@s+zh,self.bp.T@s+zp]

    def diagonal(self):
        if self._diag is None:
            # Analytical first-order diagonal; distinct-pair normalization.
            dh=self.dh[tuple(self.ip.T)].copy();dp=self.dp[tuple(self.ea.T)].copy()
            if self.spec.ip_interaction:
                for k,(i,j,a) in enumerate(self.ip):
                    dh[k]+=-self.goooo[i,j,i,j]+self.gvoov[a,i,a,i]+self.gvoov[a,j,a,j]
            if self.spec.ea_interaction:
                for k,(i,a,b) in enumerate(self.ea):
                    g=self.ints.g([self.ints.nocc+a],[self.ints.nocc+b],[self.ints.nocc+a],[self.ints.nocc+b]).item()
                    dp[k]+=g-self.govov[i,a,i,a]-self.govov[i,b,i,b]
            if self.correction is not None:
                ch,cp=self.correction.diagonal(self.ip,self.ea);dh+=ch;dp+=cp
            self._diag=np.r_[np.diag(self.a),dh,dp]
        return self._diag.copy()

    def dense(self,max_dimension=2500):
        if self.shape[0]>max_dimension:raise MemoryError('Dense dimension guard exceeded; use kernel(targets=...).')
        eye=np.eye(self.shape[0]);return np.column_stack([self.matvec(x) for x in eye])
