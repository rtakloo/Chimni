# Exact consecutive listings from Appendix A.
from fractions import Fraction

def diagonal_data(m):
    for a in range(m + 1):
        for b in range(min(2*a, m-a) + 1):
            for c in range(min(a+b, m-a-b) + 1):
                d = m-a-b-c
                if not 0 <= d <= min(a+c, 2*b):
                    continue
                e = min(a, 2*a-b)
                h = min(d, e)
                delta = min(c+d, c+b, c+a, a+b,
                            d+a, d+2*a-b, 2*a,
                            2*a+b-c, 3*a-c)
                E = min(b+c+d+delta+h,
                        b+c+2*d+min(c+h, e))
                yield (a,b,c,d), E

counts = [1, 3, 6, 11, 17, 27, 35, 50, 64, 87]
nexceptional = [0, 1, 1, 1, 1, 2, 0, 1, 0, 1]
thresholds = [(1,1), (1,1), (4,3), (5,4), (7,5),
              (4,3), (10,7), (11,8), (13,9), (7,5)]
exceptional = []
for m in range(1, 11):
    rows = list(diagonal_data(m))
    bad = [typ for typ, E in rows if 2*E >= 3*m-2]
    good = [E for typ, E in rows if 2*E < 3*m-2]
    assert len(rows) == counts[m-1]
    assert len(bad) == nexceptional[m-1]
    assert max(Fraction(E+1, m) for E in good) == \
           Fraction(*thresholds[m-1])
    exceptional.extend(bad)
assert sum(counts) == 301
assert exceptional == [
    (1,0,1,0), (1,1,0,1), (1,1,1,1), (1,1,1,2),
    (2,1,1,2), (2,2,0,2), (2,2,1,3), (2,2,2,4)]


expected_counts = [108, 136, 165, 202, 242, 290]
expected_maxima = [15, 16, 18, 19, 21, 22]
for m in range(11, 17):
    rows = list(diagonal_data(m))
    assert len(rows) == expected_counts[m-11]
    assert max(E for typ, E in rows) == expected_maxima[m-11]
    if m == 15:
        assert [typ for typ, E in rows if E == 21] == [(3,3,3,6)]
        assert max(E for typ, E in rows
                   if typ != (3,3,3,6)) == 20

print("PASS all printed finite comparisons")
