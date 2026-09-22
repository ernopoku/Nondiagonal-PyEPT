import numpy as np
import pytest
from nondiagonal_ept.integrals import Integrals
from nondiagonal_ept.blocks import Hamiltonian,amplitudes,vertices
from oracle import operator_matrix,physical_hamiltonian,reference_states

@pytest.fixture(params=[(3,1),(4,2)])
def ints(request):
    rng=np.random.default_rng(31);n,no=request.param
    z=rng.normal(size=(n,n,n,n))*.08
    z=(z+z.swapaxes(0,1)+z.swapaxes(2,3)+z.swapaxes(0,1).swapaxes(2,3))/4
    z=(z+z.transpose(2,3,0,1))/2
    return Integrals(np.repeat(np.r_[np.linspace(-1.4,-.8,no),np.linspace(.3,.9,n-no)],2),z,np.repeat(np.arange(n),2),np.tile([0,1],n),2*no,np.eye(n),np.arange(n),np.eye(n))

def test_vertices_against_fermion_algebra(ints):
    t,s=amplitudes(ints);b,c,d,e=vertices(ints,t)
    h=physical_hamiltonian(ints);ref,corr=reference_states(ints,t)
    o=ints.nocc;n=ints.n
    for p in range(ints.n):
        annihilator=operator_matrix(ints.n,((p,0),))
        comm=annihilator@h-h@annihilator
        for i,j,a in [(i,j,a) for i in range(o) for j in range(i+1,o) for a in range(o,n)]+[(a,b,i) for a in range(o,n) for b in range(a+1,n) for i in range(o)]:
            triple=operator_matrix(ints.n,((i,1),(j,1),(a,0)))
            anticom=comm@triple+triple@comm
            first=ref@anticom@ref
            second=ref@anticom@corr+corr@anticom@ref
            if i<o:
                np.testing.assert_allclose(first,b[p,i,j,a-o],atol=1e-13)
                np.testing.assert_allclose(second,c[p,i,j,a-o],atol=1e-13)
            else:
                np.testing.assert_allclose(first,d[p,a,i-o,j-o],atol=1e-13)
                np.testing.assert_allclose(second,e[p,a,i-o,j-o],atol=1e-13)

def test_triple_blocks_against_determinants(ints):
    ham=Hamiltonian(ints,'2ph-TDA');h=physical_hamiltonian(ints)
    ref=np.zeros(2**ints.n);ref[(1<<ints.nocc)-1]=1;eref=ref@h@ref
    ih=[];ip=[]
    for i,j,a in ham.ip:
        ih.append(operator_matrix(ints.n,((a+ints.nocc,1),(j,0),(i,0)))@ref)
    for i,a,b in ham.ea:
        ip.append(operator_matrix(ints.n,((a+ints.nocc,1),(b+ints.nocc,1),(i,0)))@ref)
    bh=np.array(ih).T;bp=np.array(ip).T
    dense=ham.dense()
    np.testing.assert_allclose(dense[ham.ns:ham.ns+ham.nh,ham.ns:ham.ns+ham.nh],eref*np.eye(ham.nh)-bh.T@h@bh,atol=1e-13)
    np.testing.assert_allclose(dense[-ham.np:,-ham.np:],bp.T@h@bp-eref*np.eye(ham.np),atol=1e-13)
    np.testing.assert_allclose(ham.diagonal(),np.diag(dense),atol=1e-13)

@pytest.mark.parametrize('method',['ND2','2ph-TDA','NR2','NRP3','NRQ3','NRL3','3+'])
def test_hermiticity_and_spin(ints,method):
    a=Hamiltonian(ints,method,spin=0).dense();b=Hamiltonian(ints,method,spin=1).dense()
    np.testing.assert_allclose(a,a.T,atol=1e-13)
    np.testing.assert_allclose(np.linalg.eigvalsh(a),np.linalg.eigvalsh(b),atol=1e-13)

def test_linear_triple_corrections_against_metric(ints):
    from nondiagonal_ept.renormalization import TripleCorrection
    t,_=amplitudes(ints);correction=TripleCorrection(ints,t);ham=Hamiltonian(ints,'NRL3')
    h=physical_hamiltonian(ints);ref,corr=reference_states(ints,t);o=ints.nocc
    operators=[operator_matrix(ints.n,((i,1),(j,1),(o+a,0))) for i,j,a in ham.ip]
    operators += [operator_matrix(ints.n,((o+a,1),(o+b,1),(i,0))) for i,a,b in ham.ea]
    raw=np.zeros((len(operators),)*2)
    for i,x in enumerate(operators):
        for j,y in enumerate(operators):
            comm=h@y-y@h;op=x.T@comm+comm@x.T
            raw[i,j]=ref@op@corr+corr@op@ref
    columns=[]
    for vec in np.eye(ham.nh+ham.np):
        x,y=ham.unpack(vec[:ham.nh],vec[ham.nh:]);z,w=correction.action(x,y)
        columns.append(np.r_[z[tuple(ham.ip.T)],w[tuple(ham.ea.T)]])
    actual=np.array(columns).T;expected=(raw+raw.T)/2
    # BD-T1 explicitly omits the cross-manifold block.
    expected[:ham.nh,ham.nh:]=0;expected[ham.nh:,:ham.nh]=0
    np.testing.assert_allclose(actual,expected,atol=1e-13)
    a,b=correction.diagonal(ham.ip,ham.ea)
    np.testing.assert_allclose(np.r_[a,b],np.diag(actual),atol=1e-13)

def test_static_density_against_reference_state(ints):
    from nondiagonal_ept.blocks import static_self_energy
    t,s=amplitudes(ints);ref,corr=reference_states(ints,t);n,o=ints.n,ints.nocc
    singles=np.zeros_like(ref)
    for i in range(o):
        for a in range(o,n):singles+=s[i,a-o]*operator_matrix(n,((a,1),(i,0)))@ref
    rho=np.zeros((n,n))
    for r in range(n):
        for q in range(n):
            op=operator_matrix(n,((r,1),(q,0)))
            rho[r,q]=corr@op@corr-(corr@corr)*(ref@op@ref)+ref@op@singles+singles@op@ref
    expected=np.einsum('prqs,rs->pq',ints.g('p','p','p','p'),rho)
    np.testing.assert_allclose(static_self_energy(ints,t,s,'quadratic'),expected,atol=1e-13)

def test_singles_against_rayleigh_schrodinger(ints):
    t,s=amplitudes(ints);ref,corr=reference_states(ints,t);h=physical_hamiltonian(ints)
    for i in range(ints.nocc):
        for a in range(ints.nocc,ints.n):
            single=operator_matrix(ints.n,((a,1),(i,0)))@ref
            expected=single@h@corr/(ints.energy[i]-ints.energy[a])
            np.testing.assert_allclose(s[i,a-ints.nocc],expected,atol=1e-13)

def test_strict_adc3_matches_fci_through_third_order(ints):
    from dataclasses import replace
    n,o=ints.n,ints.nocc
    neutral=[d for d in range(2**n) if d.bit_count()==o]
    ion=[d for d in range(2**n) if d.bit_count()==o-1]
    errors=[]
    for scale in (.5,.25):
        scaled=replace(ints,spatial_eri=ints.spatial_eri*scale)
        physical=physical_hamiltonian(scaled)
        exact=np.linalg.eigvalsh(physical[np.ix_(neutral,neutral)])[0]-np.linalg.eigvalsh(physical[np.ix_(ion,ion)])[0]
        ham=Hamiltonian(scaled,'3+');values,vectors=np.linalg.eigh(ham.dense())
        k=np.argmax(abs(vectors[o//2-1])**2)
        errors.append(abs(values[k]-exact))
    # Halving the fluctuation potential reduces the leading O(lambda^4)
    # error by approximately 16. This checks scientific order, not code shape.
    assert 12 < errors[0]/errors[1] < 20


def test_bdt1_vertices_against_half_weight_operator_expansion(ints):
    # Evaluate the metric in determinant space, independently of tensor vertices.
    from dataclasses import replace
    ints = replace(ints, fock=np.diag(ints.energy))
    t, _ = amplitudes(ints)
    ham = Hamiltonian(ints, 'BD-T1', doubles=t)
    h = physical_hamiltonian(ints)
    ref, corr = reference_states(ints, t)
    o = ints.nocc
    simple = [operator_matrix(ints.n, ((p, 1),)) for p in ham.simple]
    triples = [operator_matrix(ints.n, ((i,1),(j,1),(o+a,0))) for i,j,a in ham.ip]
    triples += [operator_matrix(ints.n, ((o+a,1),(o+b,1),(i,0))) for i,a,b in ham.ea]
    expected = np.empty((len(simple), len(triples)))
    for i, x in enumerate(simple):
        for j, y in enumerate(triples):
            comm = x.T@h-h@x.T
            op = comm@y+y@comm
            zeroth = ref@op@ref
            linear = ref@op@corr+corr@op@ref
            expected[i,j] = zeroth + .5*linear
    np.testing.assert_allclose(np.hstack([ham.bh,ham.bp]), expected, atol=1e-13)
