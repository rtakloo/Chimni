from collections import defaultdict

def mon(a,b,c=1): return {(a,b):c} if c else {}
def add(*ps):
 r=defaultdict(int)
 for p in ps:
  for m,c in p.items(): r[m]+=c
 return {m:c for m,c in r.items() if c}
def neg(p): return {m:-c for m,c in p.items()}
def mul(p,q):
 r=defaultdict(int)
 for (a,b),c in p.items():
  for (d,e),f in q.items(): r[a+d,b+e]+=c*f
 return {m:c for m,c in r.items() if c}
def factors(fs):
 r=mon(0,0)
 for f in fs:r=mul(r,add(mon(0,0),mon(*f,-1)))
 return r

def check(terms):
 fs=set(f for num,den in terms for f in den)
 out={}
 for num,den in terms:
  assert len(set(den))==len(den)
  out=add(out,mul(num,factors(fs-set(den))))
 assert not out,out
 print('PASS exact polynomial identity:',len(terms),'rational summands,',len(fs),'distinct binomial factors')
N=add(mon(0,0),mon(1,2),mon(2,3),mon(4,5,-1),mon(5,6,-2),mon(6,7,-2),mon(7,8,-1),mon(9,10),mon(10,11),mon(11,13))
D=[(0,1),(2,2),(3,3),(4,4),(5,6)]
F=[(mon(4,4),[(1,1),(2,2),(3,3)]),(mon(5,5,-1),[(1,1),(2,2),(4,4)]),
(mon(0,0),[(2,1),(0,1),(1,2)]),(mon(2,1,-1),[(2,1),(2,2),(3,3)]),
(mon(4,5),[(1,1),(3,3),(3,4)]),(mon(6,7,-1),[(1,1),(4,4),(5,6)]),
(mon(3,4),[(1,1),(3,3),(2,3)]),(mon(4,5,-1),[(1,1),(3,3),(3,4)]),
(mon(2,3),[(2,1),(1,2),(2,3)]),(mon(4,4,-1),[(2,1),(3,3),(2,3)])]
check(F+[(neg(N),D)])
M=add(mon(0,0),mon(0,1),mon(3,3),mon(2,3,-1),mon(3,4,-1),mon(4,5,-1))
A=[(add(mon(0,0),mon(0,1),mon(1,2)),[(2,3)]),
(add(mon(2,2),mon(1,2,-1)),[(0,1),(2,2)])]
pref=add(mon(3,2),mon(2,2,-2),mon(1,2))
for num,den in [(mon(0,2),[(0,1)]),(mon(4,4,-1),[(2,2)]),(mon(1,3,-1),[(1,2)]),(mon(5,5),[(3,3)])]:
 A.append((mul(pref,num),[(1,1),(2,1)]+den))
A.extend([(add(mon(5,6),mon(4,6,-2),mon(3,6)),[(1,2),(3,3),(2,3)]),(add(mon(3,3),mon(2,3,-2),mon(1,3)),[(1,2),(3,3)]),(add(mon(4,4),mon(3,4,-1)),[(3,3),(4,4)])])
check(A+[(neg(M),[(2,2),(4,4)])])

# Exact finite comparison as stated, testing both totals and maxima.
counts=[];maxima=[];exceptional=[]
for m in range(1,17):
 rows=[]
 for a in range(m+1):
  for b in range(min(2*a,m-a)+1):
   for c in range(min(a+b,m-a-b)+1):
    d=m-a-b-c
    if not 0<=d<=min(a+c,2*b):continue
    e=min(a,2*a-b);h=min(d,e)
    delta=min(c+d,c+b,c+a,a+b,d+a,d+2*a-b,2*a,2*a+b-c,3*a-c)
    E=min(b+c+d+delta+h,b+c+2*d+min(c+h,e))
    rows.append((E,(a,b,c,d)))
    if m<=10 and 2*E>=3*m-2:exceptional.append((a,b,c,d))
 counts.append(len(rows));maxima.append(max(E for E,t in rows))
 if m==15: assert [t for E,t in rows if E==21]==[(3,3,3,6)]
assert sum(counts[:10])==301
assert counts[10:]==[108,136,165,202,242,290]
assert maxima[10:]==[15,16,18,19,21,22]
assert exceptional==[(1,0,1,0),(1,1,0,1),(1,1,1,1),(1,1,1,2),(2,1,1,2),(2,2,0,2),(2,2,1,3),(2,2,2,4)]
print('PASS finite diagonal comparison:',counts,maxima)
