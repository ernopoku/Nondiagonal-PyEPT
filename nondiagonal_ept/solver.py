"""Residual-controlled Davidson solver and spectral observables."""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator, minres

HARTREE_TO_EV=27.211386245981

class ConvergenceError(RuntimeError): pass

@dataclass
class Pole:
    energy: float
    strength: float
    dyson_mo: np.ndarray
    dyson_ao: np.ndarray
    residual: float
    iterations: int
    target: int | None
    sector: str
    vector: np.ndarray

    @property
    def energy_ev(self): return self.energy*HARTREE_TO_EV
    @property
    def binding_energy(self):
        """-omega: positive IP for a bound hole, positive EA for a bound attachment."""
        return -self.energy
    @property
    def binding_energy_ev(self): return self.binding_energy*HARTREE_TO_EV
    @property
    def normalized_dyson_mo(self):
        return self.dyson_mo/np.sqrt(self.strength) if self.strength>1e-15 else np.zeros_like(self.dyson_mo)


def davidson(ham,index,tol=1e-9,max_cycle=150,max_space=40):
    """Track the Ritz root with maximum overlap on the requested simple MO.

    This targets primary states; use dense_spectrum for every satellite.
    Failure is an exception, never a silently returned unconverged pole.
    """
    n=ham.shape[0]
    if not 0<=index<ham.ns:raise ValueError('Target must be in the simple-operator space.')
    if max_space<4:raise ValueError('max_space must be at least four.')
    if tol<=0 or max_cycle<1:raise ValueError('tol and max_cycle must be positive.')
    q=np.zeros(n)
    q[index]=1.;basis=[q];actions=[ham.matvec(q)];diag=ham.diagonal()
    for cycle in range(1,max_cycle+1):
        b=np.column_stack(basis);ab=np.column_stack(actions)
        projected=b.T@ab
        if np.max(np.abs(projected-projected.T))>1e-8:
            raise ArithmeticError('Hamiltonian failed Hermiticity check.')
        vals,u=eigh((projected+projected.T)/2)
        overlaps=np.abs(b[index]@u)**2
        k=np.argmax(overlaps);omega=vals[k];x=b@u[:,k];ax=ab@u[:,k]
        r=ax-omega*x;res=np.linalg.norm(r)
        if res<tol:return float(omega),x,float(res),cycle
        if len(basis)>=min(n,max_space):
            # Retain target and nearby Ritz vectors. Always reorthogonalize.
            selected=[k]+[j for j in np.argsort(np.abs(vals-omega)) if j!=k][:min(5,max_space//2)-1]
            basis=[b@u[:,j] for j in selected];actions=[ab@u[:,j] for j in selected]
        den=omega-diag
        den=np.where(np.abs(den)<1e-5,np.where(den<0,-1e-5,1e-5),den)
        q=r/den
        for _ in range(2):
            for v in basis:q-=v*np.dot(v,q)
        norm=np.linalg.norm(q)
        if norm<1e-12:
            q=r.copy()
            for _ in range(2):
                for v in basis:q-=v*np.dot(v,q)
            norm=np.linalg.norm(q)
        if norm<1e-13:raise ConvergenceError(f'Davidson stagnated at residual {res:.3e}.')
        q/=norm;basis.append(q);actions.append(ham.matvec(q))
    raise ConvergenceError(f'Davidson did not converge after {max_cycle} cycles; residual {res:.3e}.')


def make_pole(ham,omega,vector,residual,iterations,target=None):
    raw=vector[:ham.ns].copy();strength=float(raw@raw)
    if raw.size and raw[np.argmax(np.abs(raw))]<0: vector=-vector;raw=-raw
    mo=np.zeros(ham.ints.coefficients.shape[1]);mo[ham.ints.spatial[ham.simple]]=raw
    ao=ham.ints.coefficients@mo
    norm=ao@ham.ints.overlap@ao
    if not np.isclose(norm,strength,atol=1e-8):raise ArithmeticError('AO Dyson norm does not equal pole strength.')
    return Pole(float(omega),strength,mo,ao,float(residual),iterations,target,ham.sector,vector)


def self_energy(ham,energy,derivative=False,tol=1e-10):
    """Schur-complement self-energy, with dSigma/dE = -B R^2 B.T.

    Uses residual-checked MINRES rather than an explicit inverse. A singular or
    poorly conditioned triple resolvent raises ConvergenceError.
    """
    n=ham.nh+ham.np;b=np.hstack([ham.bh,ham.bp])
    def action(x):
        h,p=ham.triple_action(x[:ham.nh],x[ham.nh:]);return energy*x-np.r_[h,p]
    op=LinearOperator((n,n),matvec=action,dtype=float)
    sols=[]
    for rhs in b:
        if np.linalg.norm(rhs)==0:sols.append(np.zeros(n));continue
        x,info=minres(op,rhs,rtol=tol,maxiter=max(200,3*n))
        residual=np.linalg.norm(action(x)-rhs)
        if info!=0 or residual>max(1e-9,tol*100*np.linalg.norm(rhs)):
            raise ConvergenceError(f'Triple resolvent solve failed: info={info}, residual={residual:.3e}.')
        sols.append(x)
    x=np.array(sols)
    if derivative:return -x@x.T
    return ham.a-np.diag(ham.ints.energy[ham.simple])+b@x.T
