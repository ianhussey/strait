# Reported-tuple level: find all tuples that pass GRIM + SD-bounds + GRIMMER(parity)
# yet are integer-unachievable ("phantoms"), and map their location on the umbrella.
import math
def ach_sets(W, n):
    cur = {(0,0)}
    for i in range(n):
        nxt = set()
        for (S,Q) in cur:
            for v in range(W+1):
                nxt.add((S+v, Q+v*v))
        cur = nxt
    out = {}
    for (S,Q) in cur:
        out.setdefault(S, set()).add(Q)
    return out
def qmin_i(S,n):
    m,r = divmod(S,n); return (n-r)*m*m + r*(m+1)*(m+1)
def qmax_i(S,W,n):
    a,b = divmod(S,W); return a*W*W + b*b

def phantom_map(l, u, n, d, sd_d=None):
    W = u - l
    sd_d = sd_d if sd_d is not None else d
    A = ach_sets(W, n)
    hm = 0.5*10**-d; hs = 0.5*10**-sd_d
    res = []
    m_grid = [round(l + i*10**-d, d) for i in range(int((u-l)*10**d)+1)]
    # sd grid up to max possible sd overall
    sdmax_glob = (W/2)*math.sqrt(n/(n-1))*1.001
    s_grid = [round(i*10**-sd_d, sd_d) for i in range(int(sdmax_glob*10**sd_d)+2)]
    for m in m_grid:
        my = m - l
        Slo = math.ceil((my-hm)*n); Shi = math.floor((my+hm)*n)
        # closed interval both ends (union of half rules): use <= on both
        Ss = [S for S in range(Slo, Shi+1) if 0 <= S <= n*W]
        if not Ss: continue  # GRIM fail
        for s in s_grid:
            # per S: Q window from sd in [s-hs, s+hs]
            grimmer_pass = False; achieved = False; bounds_pass=False
            for S in Ss:
                q_lo = (n-1)*max(0.0,(s-hs))**2 + S*S/n
                q_hi = (n-1)*(s+hs)**2 + S*S/n
                Qlo = math.ceil(q_lo - 1e-9); Qhi = math.floor(q_hi + 1e-9)
                Qlo = max(Qlo, qmin_i(S,n)); Qhi = min(Qhi, qmax_i(S,W,n))
                if Qlo > Qhi: continue
                bounds_pass = True
                # parity
                if (Qlo - S) % 2: Qlo += 1
                if Qlo > Qhi: continue
                grimmer_pass = True
                for Q in range(Qlo, Qhi+1, 2):
                    if Q in A[S]: achieved=True; break
                if achieved: break
            if grimmer_pass and not achieved:
                res.append((m, s))
    return res

for (l,u,n,d) in [(1,5,9,1),(1,7,9,1),(1,7,25,2),(0,10,15,2),(1,7,50,2)]:
    ph = phantom_map(l,u,n,d)
    print(f"=== scale {l}-{u}, n={n}, digits={d}: {len(ph)} phantom tuples ===")
    # summarize by mean
    from collections import defaultdict
    bym = defaultdict(list)
    for (m,s) in ph: bym[m].append(s)
    for m in sorted(bym):
        ss = bym[m]
        print(f"  mean={m}: sd in {sorted(ss)[:12]}{'...' if len(ss)>12 else ''}")
