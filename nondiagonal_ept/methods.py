"""Table IV, JCP 159, 124109 (2023). Orders refer to whole blocks."""
from dataclasses import dataclass

@dataclass(frozen=True)
class Method:
    name: str
    static: str = 'none'
    ip_vertex_h: float = 0.
    ip_vertex_p: float = 0.
    ea_vertex_h: float = 0.
    ea_vertex_p: float = 0.
    ip_interaction: bool = False
    ea_interaction: bool = False
    triple_correction: bool = False

METHODS = {
    'BD-T1': Method('BD-T1','none',.5,.5,.5,.5,True,True,True),
    'ND2': Method('ND2'),
    '2PH-TDA': Method('2ph-TDA',ip_interaction=True,ea_interaction=True),
    'NR2': Method('NR2',ip_vertex_h=.5,ip_interaction=True),
    'NRP3': Method('NRP3',ip_vertex_h=.5,ip_vertex_p=.5,ip_interaction=True),
    'NRQ3': Method('NRQ3','linear',.5,.5,0.,0.,True,False),
    'NRL3': Method('NRL3','linear',.5,.5,.5,.5,True,True),
    '3+': Method('3+','quadratic',1.,1.,1.,1.,True,True),
}

def method_spec(name, sector='ip'):
    name=name.upper().replace(' ','')
    if name=='ADC(2)':name='ND2'
    if name in ('ADC(3)-STRICT','STRICT-ADC(3)'): name='3+'
    if name=='ADC(3)':
        raise ValueError('ADC(3) is ambiguous. Use ADC(3)-strict for 3+; renormalized static ADC(3) is not implemented.')
    if name not in METHODS: raise ValueError(f'Unknown/unsupported method {name!r}; choose {list(METHODS)}.')
    if sector not in ('ip','ea'): raise ValueError("sector must be 'ip' or 'ea'.")
    m=METHODS[name]
    if sector=='ea' and name in ('NR2','NRP3','NRQ3'):
        # Particle-hole counterpart of the detachment-oriented truncation.
        return Method(m.name,m.static,m.ea_vertex_p,m.ea_vertex_h,
                      m.ip_vertex_p,m.ip_vertex_h,m.ea_interaction,m.ip_interaction)
    return m
