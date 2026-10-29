"""Exact family and analytic probability results for validated inputs."""
import math
from decimal import Decimal, localcontext
from fractions import Fraction

import numpy as np

from .analytic import analytic
from .bounds import global_witness
from .conditioned import evaluate, solve


def exact_number(value):
    """Keep an exact rational, a non-underflowing decimal, and its log."""
    if value is None:
        return {'numerator':None,'denominator':None,'integer_base':None,'decimal':None,'log':None}
    if not isinstance(value,Fraction):
        value=Fraction(value)
    with localcontext() as ctx:
        ctx.prec=40
        decimal=str(Decimal(value.numerator)/Decimal(value.denominator))
    base=16 if max(value.numerator.bit_length(),value.denominator.bit_length())>2000 else 10
    encode=hex if base==16 else str
    return {'numerator':encode(value.numerator),'denominator':encode(value.denominator),'integer_base':base,'decimal':decimal,
            'log':math.log(value.numerator)-math.log(value.denominator) if value else None}


def analyze(case):
    sampler=case['sampler']
    context,records=solve(sampler)
    lower=sampler.counts(global_witness(sampler.neff,sampler.rna,False))
    must=bool((lower>sampler.maxima).any())
    can=context['possible']
    exact=analytic(sampler.neff,sampler.rna,sampler.maxima)
    if (exact['kind']=='exact_impossible') != (not can) or (exact['kind']=='exact_certain') != must:
        raise ValueError('analytic and family gate states disagree')
    for row in records:
        row['default_core']=bool(row['raw_core'] and must)
        row['default_union']=bool(row['raw_union'] and can)
    states={'gate_state':'impossible' if not can else 'certain' if must else 'possible_not_certain',
            'original_gate_pass':bool((case['original_counts']>sampler.maxima).any())
                                 if case['original_counts'] is not None else None,
            'original_gate_status':'evaluated' if case['original_counts'] is not None else 'not_evaluated_missing_DNA_rank',
            'raw_core':int(context['bound']['must'].sum()),'raw_union':int(context['bound']['possible'].sum()),
            'default_core':int(context['bound']['must'].sum()) if must else 0,
            'default_union':int(context['bound']['possible'].sum()) if can else 0,
            'passing_core':sum(bool(r['conditioned_core']) for r in records) if can else None,
            'passing_union':sum(bool(r['conditioned_union']) for r in records) if can else None,
            'passing_family_status':'exact' if can else 'undefined_G_impossible',
            'probability_method':exact['kind'],'gate_probability':exact_number(exact['gate_probability'])}
    return context,records,exact,states


def member_witness(context,cell,selected):
    sampler=context['s']
    positions=np.flatnonzero(sampler.ids==cell)
    if len(positions)!=1:raise ValueError('Requested member ID is not in the case')
    result=evaluate(context,int(positions[0]),selected)
    return {'cell':cell,'selected':selected,'status':result['status'],'reachable':result['reachable'],
            'reason':result['reason'],'gate_counts':result['counts'].tolist() if result['counts'] is not None else None,
            'rank':result['rank']}
