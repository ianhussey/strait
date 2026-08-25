# Window-level analytic certification: is there an achievable Q in [Qlo, Qhi]?
# hit(S, cap, n, dlo, dhi): exists achievable deficit in [dlo, dhi] measured from
# qmax(S,cap,n), via demotion decomposition (t parts at cap, rest <= cap-1). Exact.
import math, random, sys
sys.setrecursionlimit(100000)

def qmin_i(S, n):
    m, r = divmod(S, n)
    return (n-r)*m*m + r*(m+1)*(m+1)

def qmax_i(S, cap, n):
    if cap == 0: return 0 if S == 0 else None
    a, b = divmod(S, cap)
    if a > n or (a == n and b > 0): return None
    return a*cap*cap + b*b

def hit(S, cap, n, dlo, dhi):
    """exists partition of S, <=n parts, parts<=cap, with deficit from qmax in [dlo,dhi]?"""
    if S == 0: return dlo <= 0 <= dhi
    if cap == 0 or n <= 0: return False
    if cap > S: cap = S
    Qx = qmax_i(S, cap, n)
    if Qx is None: return False
    if dhi < 0: return False
    # max possible deficit
    if dlo > Qx - qmin_i(S, n): return False
    if dlo <= 0 <= dhi: return True   # the extremal config itself
    a, b = divmod(S, cap)
    # j = number of demoted cap-parts (t = a - j at cap); j=0 branch: re-partition b below cap
    for j in range(0, a + 1):
        S_j = j*cap + b
        n_j = n - (a - j)
        if n_j <= 0:
            continue
        qx = qmax_i(S_j, cap - 1, n_j)
        if qx is None:
            continue
        base = j*cap*cap + b*b     # deficit when sub-config q = qx is 'base - qx'
        d_off = base - qx
        if d_off > dhi:
            continue               # cheap probe; keep scanning j (d_off not proven monotone)
        if hit(S_j, cap - 1, n_j, dlo - d_off, dhi - d_off):
            return True
    return False

def window_possible(S, Qlo, Qhi, W, n):
    if 2*S > n*W:
        shift = n*W*W - 2*W*S
        S = n*W - S
        Qlo += shift; Qhi += shift
    Qx = qmax_i(S, W, n)
    if Qx is None: return False
    Qlo = max(Qlo, qmin_i(S, n)); Qhi = min(Qhi, Qx)
    if (Qlo - S) % 2: Qlo += 1
    if (Qhi - S) % 2: Qhi -= 1
    if Qlo > Qhi: return False
    return hit(S, W, n, Qx - Qhi, Qx - Qlo)

if __name__ == "__main__":
    # validation vs DP: (1) point windows == membership, (2) random windows
    def ach_sets(W, n):
        cur = {(0,0)}
        for _ in range(n):
            cur = {(S+v, Q+v*v) for (S,Q) in cur for v in range(W+1)}
        out = {}
        for (S,Q) in cur: out.setdefault(S, set()).add(Q)
        return out
    import time
    t0=time.time(); bad=tot=0
    random.seed(3)
    for W in range(1, 9):
        for n in range(2, 13):
            A = ach_sets(W, n)
            for S in range(0, n*W+1):
                qs = A[S]
                hi = max(qs); lo = min(qs)
                # point windows: all Q
                for Q in range(lo, hi+1):
                    tot += 1
                    if (Q in qs) != window_possible(S, Q, Q, W, n):
                        bad += 1; print("PT MM", W,n,S,Q)
                # random windows
                for _ in range(20):
                    x = random.randint(lo-3, hi+3); y = x + random.randint(0, 12)
                    truth = any(q in qs for q in range(max(x,0), y+1))
                    tot += 1
                    if truth != window_possible(S, x, y, W, n):
                        bad += 1; print("WIN MM", W,n,S,x,y)
    print(f"window validation: {tot} checks, {bad} mismatches, {time.time()-t0:.0f}s")
