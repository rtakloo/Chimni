#!/usr/bin/env python3
"""Independent finite checks for chimni-TB2.tex (reconstructed 2026-11-09).

Only Python's standard library is required. No supplied verification script
was available. Lattice closure is derived from actual truncated-polynomial
multiplication and exact row division. Last-column counts solve the resulting
linear systems; displayed paper congruences are only comparison targets.
Finite agreement is not a proof of any unbounded statement.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from fractions import Fraction
import argparse
import math
import time

START = time.monotonic()

def report(*args):
    print(*args, flush=True)

def mul(u, v, dimension=None):
    """Coordinates are t, t^2,...,t^dimension; no constant term."""
    d = dimension if dimension is not None else len(u)
    out = [0] * d
    for i, a in enumerate(u):
        if not a:
            continue
        for j, b in enumerate(v):
            k = i + j + 1
            if k < d:
                out[k] += a * b
    return out

def row_division(v, H):
    rem = list(v)
    coeff = []
    for i, row in enumerate(H):
        q, r = divmod(rem[i], row[i])
        if r:
            return None
        coeff.append(q)
        for j in range(i, len(H)):
            rem[j] -= q * row[j]
    assert not any(rem)
    return coeff

def closed(H):
    return all(row_division(mul(H[i], H[j]), H) is not None
               for i in range(len(H)) for j in range(i, len(H)))

def leading_admissible(exps):
    # Leading terms of every basis-row product give these necessary conditions.
    return all(exps[i] + exps[j] >= exps[i+j+1]
               for i in range(len(exps)) for j in range(i, len(exps))
               if i+j+1 < len(exps))

def lift_system(H):
    """Derive A x = b from row products, with no paper closure formulas."""
    r = len(H)
    A, b = [], []
    for i in range(r):
        for j in range(i, r):
            v = mul(H[i], H[j], r + 1)
            coefficients = row_division(v[:r], H)
            assert coefficients is not None
            if any(coefficients) or v[r]:
                A.append(coefficients)
                b.append(v[r])
    return A, b

def solve_count(A, b, p, exponent, variables):
    """Count solutions modulo p^exponent using unit pivots and p-division."""
    if exponent == 0:
        return 1
    q = p ** exponent
    if not A:
        return q ** variables
    A = [[a % q for a in row] for row in A]
    b = [v % q for v in b]
    pivot = next(((i, j) for i, row in enumerate(A)
                  for j, a in enumerate(row) if a % p), None)
    if pivot is not None:
        r, c = pivot
        inv = pow(A[r][c], -1, q)
        remaining = [j for j in range(variables) if j != c]
        Ar, br = [], []
        for i in range(len(A)):
            if i == r:
                continue
            factor = A[i][c] * inv % q
            Ar.append([(A[i][j] - factor*A[r][j]) % q for j in remaining])
            br.append((b[i] - factor*b[r]) % q)
        return solve_count(Ar, br, p, exponent, variables-1)
    if any(v % p for v in b):
        return 0
    return p ** variables * solve_count([[a//p for a in row] for row in A],
                                         [v//p for v in b], p, exponent-1, variables)

def solve_vectors(A, b, p, exponent, variables):
    """Enumerate the same modular system by independent reconstruction."""
    if exponent == 0:
        yield (0,) * variables
        return
    q = p ** exponent
    if not A:
        yield from product(range(q), repeat=variables)
        return
    A = [[a % q for a in row] for row in A]
    b = [v % q for v in b]
    pivot = next(((i, j) for i, row in enumerate(A)
                  for j, a in enumerate(row) if a % p), None)
    if pivot is not None:
        r, c = pivot
        inv = pow(A[r][c], -1, q)
        remaining = [j for j in range(variables) if j != c]
        Ar, br = [], []
        for i in range(len(A)):
            if i == r:
                continue
            factor = A[i][c] * inv % q
            Ar.append([(A[i][j]-factor*A[r][j]) % q for j in remaining])
            br.append((b[i]-factor*b[r]) % q)
        for reduced in solve_vectors(Ar, br, p, exponent, variables-1):
            x = [0] * variables
            for j, value in zip(remaining, reduced):
                x[j] = value
            x[c] = ((b[r] - sum(A[r][j]*x[j] for j in remaining))*inv) % q
            yield tuple(x)
        return
    if any(v % p for v in b):
        return
    q0 = q // p
    for reduced in solve_vectors([[a//p for a in row] for row in A],
                                 [v//p for v in b], p, exponent-1, variables):
        for digits in product(range(p), repeat=variables):
            yield tuple(x+q0*t for x, t in zip(reduced, digits))

@lru_cache(maxsize=None)
def lattices(p, exps):
    if len(exps) == 1:
        return (((p**exps[0],),),)
    if not leading_admissible(exps):
        return ()
    d = len(exps)
    out = []
    for H in lattices(p, exps[:-1]):
        A, b = lift_system(H)
        for x in solve_vectors(A, b, p, exps[-1], d-1):
            out.append(tuple(tuple(row)+(value,) for row, value in zip(H, x))
                       + ((0,)*(d-1)+(p**exps[-1],),))
    return tuple(out)

def rank_mod(H, p):
    rows = [[v % p for v in row] for row in H]
    rank = 0
    for col in range(len(rows[0]) if rows else 0):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inv = pow(rows[rank][col], -1, p)
        rows[rank] = [v*inv % p for v in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][col]:
                factor = rows[i][col]
                rows[i] = [(a-factor*b) % p for a,b in zip(rows[i],rows[rank])]
        rank += 1
    return rank

def left_null_vector(H, p):
    r = len(H)
    # Used only for nullity one; search has at most p^3 cases here.
    return next(y for y in product(range(p), repeat=r)
                if any(y) and all(sum(y[i]*H[i][j] for i in range(r)) % p == 0
                                  for j in range(r)))

def diagonal_count(p, exps, cyclic=False):
    if not leading_admissible(exps):
        return 0
    if len(exps) == 1:
        return 1
    total = 0
    d = len(exps)
    a = exps[-1]
    for H in lattices(p, exps[:-1]):
        A, b = lift_system(H)
        count = solve_count(A, b, p, a, d-1)
        if not count:
            continue
        if cyclic:
            rank = rank_mod(H, p)
            if a == 0:
                count = count if rank >= d-2 else 0
            elif rank < d-2:
                count = 0
            elif rank == d-2:
                y = left_null_vector(H, p)
                bad = solve_count(A + [[v*p**(a-1) for v in y]], b+[0], p, a, d-1)
                count -= bad
        total += count
    return total

def compositions(total, parts):
    if parts == 1:
        yield (total,)
    else:
        for i in range(total+1):
            for rest in compositions(total-i, parts-1):
                yield (i,) + rest

def coefficients(n, p, bound):
    total, cyc = [], []
    for m in range(bound+1):
        types = list(compositions(m,n-1))
        total.append(sum(diagonal_count(p,a) for a in types))
        cyc.append(sum(diagonal_count(p,a,True) for a in types))
    return total, cyc

def rational_series(numerator, factors, bound):
    values = [0]*(bound+1)
    values[0] = 1
    for degree, coefficient in factors:
        for k in range(degree,bound+1):
            values[k] += coefficient*values[k-degree]
    return [sum(v*values[k-i] for i,v in numerator.items() if i <= k)
            for k in range(bound+1)]

def quartic_expected(p, bound, cyclic=False):
    if cyclic:
        numerator = {0:1,1:1,3:p**3-p**2,4:-p**3,5:-p**4}
        if p == 2:
            numerator={0:1,1:1,3:4,4:-16,5:-16,6:64,7:-64,9:512}
        return rational_series(numerator,[(2,p**2),(4,p**4)],bound)
    numerator={0:1,2:p,3:p**2,5:-p**4,6:-2*p**5,7:-2*p**6,
               8:-p**7,10:p**9,11:p**10,13:p**11}
    expected=rational_series(numerator,[(1,1),(2,p**2),(3,p**3),(4,p**4),(6,p**5)],bound)
    if p==2:
        correction=rational_series({6:64},[(2,4),(3,8),(4,16)],bound)
        expected=[a+b for a,b in zip(expected,correction)]
    return expected

def rho(p,j):
    q=p**j
    return sum((z*z+z-1)%q==0 for z in range(q))

def cc5_formula(p,m):
    if m==0:
        return 1
    first=p**(2*m-(m+3)//4-(m+1)//2)
    return first+sum((p-1)*rho(p,m-4*r)*p**(m+r-1)
                     for r in range((m+4)//5,(m-1)//4+1))

def cyclic_character_closed(values, modulus):
    """Check actual kernel generators t^j-b_j*t by multiplication."""
    d=len(values)
    rows=[]
    for j in range(1,d):
        row=[0]*d
        row[0]=-values[j]
        row[j]=1
        rows.append(row)
    return all(sum(a*b for a,b in zip(mul(rows[i],rows[j]),values))%modulus==0
               for i in range(d-1) for j in range(i,d-1))

def recurrence_pair(n,B,C,q):
    U=[0,1,B,C]
    for j in range(4,n+2):
        U.append((B*U[j-1]+(C-B*B)*U[j-2])%q)
    return U

def reduced_cc5(p,m):
    """Scaling t by a unit allows one representative B for each valuation."""
    q=p**m
    total=0
    B_values=[(0,1)]+[(p**r,(p-1)*p**(m-r-1)) for r in range(1,m)]
    for B,multiplicity in B_values:
        count=0
        for C in range(0,q,p):
            D=(2*B*C-B**3)%q
            count += cyclic_character_closed([1,B,C,D],q)
        total+=multiplicity*count
    return total

def det(M):
    if not M:
        return 1
    if len(M)==1:
        return M[0][0]
    return sum((-1)**j*M[0][j]*det([row[:j]+row[j+1:] for row in M[1:]])
               for j in range(len(M)))

def valuation(n,p):
    if n==0:
        return math.inf
    result=0
    while n%p==0:
        result+=1
        n//=p
    return result

def smith_weight(H,p):
    d=len(H)
    delta=[0]
    for k in range(1,d+1):
        delta.append(min(valuation(det([[H[i][j] for j in cols] for i in rows]),p)
                         for rows in combinations(range(d),k)
                         for cols in combinations(range(d),k)))
    return tuple(delta[i]-delta[i-1] for i in range(d,0,-1))

def cubic_expected(p,bound):
    series=Counter({(0,0):1})
    for a,b,c in [(1,1,1),(2,1,p),(3,0,p*p)]:
        previous=series
        series=Counter()
        for (x,y),v in previous.items():
            k=0
            while x+y+k*(a+b)<=bound:
                series[x+k*a,y+k*b]+=v*c**k
                k+=1
    numerator={(0,0):1,(1,0):1,(2,0):p,(3,1):-p,(4,1):-p*p,(5,1):-p*p}
    answer=Counter()
    for (x,y),v in series.items():
        for (a,b),c in numerator.items():
            if x+y+a+b<=bound:
                answer[x+a,y+b]+=v*c
    return +answer

def verify_solver():
    systems=0
    for p in (2,3):
        for a in (1,2):
            q=p**a
            for t in range(100):
                # Deterministic heterogeneous systems, including singular ones.
                A=[[(t*(i+2)+j*j+3*i*j)%q for j in range(3)] for i in range(3)]
                b=[(t*t+2*i*t+i)%q for i in range(3)]
                brute={x for x in product(range(q),repeat=3)
                       if all((sum(c*v for c,v in zip(row,x))-target)%q==0
                              for row,target in zip(A,b))}
                assert len(brute)==solve_count(A,b,p,a,3)
                assert brute==set(solve_vectors(A,b,p,a,3))
                systems+=1
    report('PASS modular linear solver versus brute force:',systems,'systems; p=2,3; exponents=1,2; 3 variables')

def displayed_closure(H,p,exps):
    """Comparison target only: the paper's displayed closure conditions."""
    if len(exps)==2:
        a,b=exps
        return b<=2*a
    if len(exps)==3:
        a,b,c=exps
        if b>2*a or c>a+b:return False
        u,v,x=H[0][1],H[0][2],H[1][2]
        return (2*p**a*u-p**(2*a-b)*x)%p**c==0
    a,b,c,d=exps
    if b>2*a or c>a+b or d>a+c or d>2*b:return False
    u,v,w,x,y,z=H[0][1],H[0][2],H[0][3],H[1][2],H[1][3],H[2][3]
    T=2*p**a*u-p**(2*a-b)*x
    if T%p**c:return False
    return ((p**a*x+p**b*u-p**(a+b-c)*z)%p**d==0 and
            (u*u+2*p**a*v-p**(2*a-b)*y-(T//p**c)*z)%p**d==0)

def displayed_smith(H,p,exps):
    """Use the explicitly listed minors as a target for generic minors."""
    generic=smith_weight(H,p)
    if len(exps)==2:return generic
    if len(exps)==3:
        a,b,c=exps
        u,v,x=H[0][1],H[0][2],H[1][2]
        d1=min(a,b,c,valuation(u,p),valuation(v,p),valuation(x,p))
        d2=min(a+b,a+valuation(x,p),a+c,valuation(u*x-p**b*v,p),
               c+valuation(u,p),b+c)
        return (a+b+c-d2,d2-d1,d1)
    a,b,c,d=exps
    u,v,w,x,y,z=H[0][1],H[0][2],H[0][3],H[1][2],H[1][3],H[2][3]
    terms=[p**(a+b+c),p**(a+b)*z,p**a*(x*z-p**c*y),
           u*x*z-u*p**c*y-p**b*v*z+p**(b+c)*w,p**(a+b+d),
           p**(a+d)*x,p**d*(u*x-p**b*v),p**(a+c+d),p**(c+d)*u,p**(b+c+d)]
    d3=min(valuation(v,p) for v in terms)
    assert d3==sum(generic[1:])
    return (sum(exps)-d3,)+generic[1:]

def verify_full_hnf():
    raw_checked=0
    closed_checked=0
    for n,p,bound in [(3,2,8),(3,3,5),(3,5,4),(4,2,5),(4,3,3),(5,2,3),(5,3,2)]:
        for m in range(bound+1):
            for exps in compositions(m,n-1):
                d=n-1
                positions=[(i,j) for j in range(d) for i in range(j)]
                observed=[]
                for entries in product(*(range(p**exps[j]) for i,j in positions)):
                    H=[[0]*d for _ in range(d)]
                    for i in range(d):H[i][i]=p**exps[i]
                    for (i,j),x in zip(positions,entries):H[i][j]=x
                    raw_checked+=1
                    is_closed=closed(H)
                    assert is_closed==displayed_closure(H,p,exps),(n,p,exps,H)
                    if is_closed:
                        assert smith_weight(H,p)==displayed_smith(H,p,exps)
                        observed.append(tuple(map(tuple,H)))
                independent=lattices(p,exps)
                assert set(observed)==set(independent),(n,p,exps)
                assert len(observed)==diagonal_count(p,exps)
                cyclic=sum(rank_mod(H,p)>=d-1 for H in observed)
                assert cyclic==diagonal_count(p,exps,True)
                closed_checked+=len(observed)
        report('PASS raw row-HNF enumeration: n,p,max_m=',n,p,bound)
    report('Raw matrices tested:',raw_checked,'; closed lattices checked:',closed_checked)
    report('PASS displayed cubic/quartic/quintic closure conditions on all raw matrices; explicit cotype weights versus generic integer minors on all closed matrices in these ranges')

def verify_cubic():
    for p,bound in [(2,10),(3,7),(5,5)]:
        observed=Counter()
        for m in range(bound+1):
            for exps in compositions(m,2):
                for H in lattices(p,exps):
                    assert closed(H)
                    observed[smith_weight(H,p)]+=1
        assert observed==cubic_expected(p,bound),(p,observed-cubic_expected(p,bound))
        report('PASS cubic full cotype via integer minors: p,max_m=',p,bound,'lattices=',sum(observed.values()))

def verify_coefficients():
    for p,bound in [(2,9),(3,7),(5,6)]:
        total,cyc=coefficients(4,p,bound)
        assert total==quartic_expected(p,bound),(p,total,quartic_expected(p,bound))
        assert cyc==quartic_expected(p,bound,True),(p,cyc,quartic_expected(p,bound,True))
        report('PASS quartic row-product linear lifting: p,max_m=',p,bound,'total=',total,'cyclic=',cyc)
    expected={2:[1,1,7,27,115,323,1411],3:[1,1,13,76,589,2371,21325],5:[1,1,31,306,5431,46681,815431]}
    for p in (2,3,5):
        total,cyc=coefficients(5,p,6)
        assert total==expected[p],(p,total,expected[p])
        assert cyc==[cc5_formula(p,m) for m in range(7)]
        report('PASS quintic row-product linear lifting: p,max_m=',p,6,'total=',total,'cyclic=',cyc)
        lattices.cache_clear()

def verify_exceptional():
    types=[(1,0,1,0),(1,1,0,1),(1,1,1,1),(1,1,1,2),(2,1,1,2),(2,2,0,2),(2,2,1,3),(2,2,2,4)]
    for p in (2,3,5,11):
        expected=[p**2,p**3,p**5,p**5*(1+(p-1)*rho(p,1)),p**8,p**7,p**10,
                  p**10*(p*p+(p*p-p)*rho(p,2))]
        # Large split-prime diagonals are counted through exact linear systems,
        # not by listing their (potentially trillion-sized) sets of lifts.
        tested_types = types if p != 11 else types[:-1]
        actual=[diagonal_count(p,exps) for exps in tested_types]
        assert actual==expected[:len(tested_types)],(p,actual,expected)
        report('PASS exceptional diagonals from row products: p=',p,'number=',len(tested_types),'counts=',actual)
        lattices.cache_clear()

def verify_family():
    for p,max_k in [(2,5),(3,3),(5,2),(11,1)]:
        values=[]
        for k in range(max_k+1):
            q=p**k
            N=sum((x*x-u*x-u*u)%q==0 for u in range(q) for x in range(q))
            expected=p**(2*(k//2))+sum((p-1)*p**(k-1)*rho(p,k-2*r) for r in range((k+1)//2)) if k else 1
            assert N==expected
            count=diagonal_count(p,(k,k,k,2*k))
            assert count==p**(5*k)*N,(p,k,count,N)
            values.append(count)
        report('PASS diagonal family via row products: p,k=0..',p,max_k,'counts=',values)
        lattices.cache_clear()
    for p,max_k in [(2,8),(3,5),(5,4),(11,3)]:
        for k in range(max_k+1):
            q=p**k
            N=sum((x*x-u*x-u*u)%q==0 for u in range(q) for x in range(q))
            if p==5:expected=p**k
            elif rho(p,1)==0:expected=p**(2*(k//2))
            else:expected=p**k+k*(p-1)*p**(k-1) if k else 1
            assert N==expected,(p,k,N,expected)
        report('PASS direct quadratic-pair enumeration and case formulas: p,k=0..',p,max_k)

def verify_universal():
    total_characters=0
    for p in (2,3,5):
        for n in range(4,9):
            max_m=3 if p==2 and n<=6 else 2
            for m in range(1,max_m+1):
                q=p**m
                actual=set()
                for tail in product(range(0,q,p),repeat=n-2):
                    total_characters+=1
                    if cyclic_character_closed([1]+list(tail),q):actual.add(tail)
                predicted=set()
                for B,C in product(range(0,q,p),repeat=2):
                    U=recurrence_pair(n,B,C,q)
                    if U[n]%q==0 and U[n+1]%q==0:predicted.add(tuple(U[2:n]))
                assert actual==predicted,(n,p,m,len(actual),len(predicted))
                if n==4:assert len(actual)==quartic_expected(p,m,True)[m]
                if n==5:assert len(actual)==cc5_formula(p,m)
        report('PASS full normalized-character product checks versus universal recurrence: p=',p,'n=4..8; m<=2 plus p=2,n<=6,m=3')
    report('Normalized characters tested directly:',total_characters)
    pair_count=0
    for p,max_m in [(2,7),(3,5),(5,4)]:
        for m in range(1,max_m+1):
            q=p**m
            for n in range(4,11):
                for B,C in product(range(0,q,p),repeat=2):
                    U=recurrence_pair(n,B,C,q)
                    actual=cyclic_character_closed(U[1:n],q)
                    predicted=U[n]%q==0 and U[n+1]%q==0
                    assert actual==predicted,(n,p,m,B,C)
                    pair_count+=1
        report('PASS recurrence-derived character product checks: p,max_m=',p,max_m,'n=4..10')
    report('Recurrence-pair checks:',pair_count)
    for p,max_m in [(2,9),(3,7),(5,6),(11,5)]:
        values=[]
        for m in range(1,max_m+1):
            actual=reduced_cc5(p,m)
            assert actual==cc5_formula(p,m),(p,m,actual,cc5_formula(p,m))
            values.append(actual)
        report('PASS quintic characters reduced by unit scaling: p,m=1..',p,max_m,'counts=',values)
    assert cc5_formula(11,5)==3382071
    assert cc5_formula(5,5)==15625
    assert rho(5,2)==0 and rho(11,1)==rho(11,2)==2

def exponent_bound(exps):
    a,b,c,d=exps
    e=min(a,2*a-b)
    h=min(d,e)
    delta=min(c+d,c+b,c+a,a+b,d+a,d+2*a-b,2*a,2*a+b-c,3*a-c)
    return min(b+c+d+delta+h,b+c+2*d+min(c+h,e))

def verify_comparison():
    exceptional=[]
    table=[]
    expected_counts=[1,3,6,11,17,27,35,50,64,87,108,136,165,202,242,290]
    for m in range(1,17):
        types=[t for t in compositions(m,4) if leading_admissible(t)]
        assert len(types)==expected_counts[m-1]
        E=[exponent_bound(t) for t in types]
        if m<=10:
            exceptional.extend(t for t,e in zip(types,E) if 2*e>=3*m-2)
            nonexceptional=[e for e in E if 2*e<3*m-2]
            table.append((m,len(types),len(E)-len(nonexceptional),str(Fraction(max(nonexceptional)+1,m))))
        else:
            assert max(E)==[15,16,18,19,21,22][m-11]
        if m==15:
            assert [t for t,e in zip(types,E) if e==21]==[(3,3,3,6)]
            assert max(e for t,e in zip(types,E) if t!=(3,3,3,6))==20
    assert exceptional==[(1,0,1,0),(1,1,0,1),(1,1,1,1),(1,1,1,2),(2,1,1,2),(2,2,0,2),(2,2,1,3),(2,2,2,4)]
    assert sum(expected_counts[:10])==301
    report('PASS finite E comparison m=1..16; m<=10 rows=',table)
    report('PASS m=11..16 counts=',expected_counts[10:],'max E=',[15,16,18,19,21,22])

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--only',choices=['solver','raw','cubic','coefficients','exceptional','family','universal','comparison'])
    args=parser.parse_args()
    suites=[('solver',verify_solver),('raw',verify_full_hnf),('cubic',verify_cubic),
            ('coefficients',verify_coefficients),('exceptional',verify_exceptional),
            ('family',verify_family),('universal',verify_universal),('comparison',verify_comparison)]
    report('Independent reconstruction; Python standard library; integer arithmetic only.')
    report('Membership engine: truncated-polynomial row products + exact triangular division.')
    report('Larger finite counts use exact linear lifting, independently derived from these products.')
    for name,fn in suites:
        if args.only is None or args.only==name:
            t=time.monotonic()
            fn()
            report('Suite',name,'completed in',round(time.monotonic()-t,3),'seconds')
    report('ALL REQUESTED CHECKS PASSED; elapsed seconds=',round(time.monotonic()-START,3))
    report('LIMIT: finite checks neither prove unbounded coefficient identities nor analytic convergence/asymptotics.')
