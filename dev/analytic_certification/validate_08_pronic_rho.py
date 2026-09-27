# Test the pronic-weight reformulation of the corner problem.
#   I = Q - S = sum p_i(p_i - 1);  attainable at sum S  <=>  rho(I) <= S,
#   where rho(I) = min total weight sum p_i over representations of I
#   as a sum of pronic numbers p(p-1), p >= 2.  (1-parts are free.)
import math

IMAX = 26000
INF = 10**9
# DP: rho[I] = min weight
rho = [INF]*(IMAX+1); rho[0] = 0
pr = []
p = 2
while p*(p-1) <= IMAX:
    pr.append(p); p += 1
for I in range(2, IMAX+1, 2):
    best = INF
    for p in pr:
        c = p*(p-1)
        if c > I: break
        w = rho[I-c] + p
        if w < best: best = w
    rho[I] = best

# (1) exact equivalence with corner spectra: A(S) = { S+I : I even, rho(I) <= S }
def corner_spectra(N):
    # B[s] bitmask of achievable sum-of-squares over partitions of s
    B = [0]*(N+1); B[0] = 1
    cur = [0]*(N+1); cur[0] = 1
    for p in range(1, N+1):
        pp = p*p
        for s in range(p, N+1):
            cur[s] |= cur[s-p] << pp
    return cur
N = 158
cur = corner_spectra(N)
bad = 0; tot = 0
for S in range(1, N+1):
    m = cur[S]
    for Q in range(S, S*S+1):
        I = Q - S
        att = (m >> Q) & 1 == 1
        pred = (I % 2 == 0) and I <= IMAX and rho[I] <= S
        if I > IMAX: continue
        tot += 1
        if att != pred:
            bad += 1
            if bad < 6: print("RHO MISMATCH", S, Q, I, rho[I], att)
print(f"[1] rho-characterisation: {tot} checks, {bad} mismatches")

# (2) greedy telescoping upper bound on rho
def greedy_weight(I):
    w = 0
    while I > 0:
        u = int((1 + math.isqrt(1 + 4*I)) // 2)
        while u*(u-1) > I: u -= 1
        while (u+1)*u <= I: u += 1
        w += u
        I -= u*(u-1)
    return w
worst_excess = 0
for I in range(2, IMAX+1, 2):
    g = greedy_weight(I)
    assert g >= rho[I]
    u1 = int((1 + math.isqrt(1 + 4*I)) // 2)
    bound = u1 + math.isqrt(2*u1) + 2*int((2*u1)**0.25) + 8
    if g > bound:
        worst_excess += 1
        if worst_excess < 5: print("greedy exceeds claimed bound", I, g, bound)
print(f"[2] greedy vs claimed bound B(I): {worst_excess} violations; "
      f"max greedy-rho gap = {max(greedy_weight(I)-rho[I] for I in range(2, IMAX+1, 400))} (sampled)")

# (3) deepest-gap formula: D*(S) = S^2 - S - min{ even I : rho(I) > S }
print("[3]  S   D*_true   2sqrt2*S^1.5   ratio   lower-bound 2mS-2m^2+2m-4")
for S in range(40, N+1, 20):
    m = cur[S]
    # deepest missing = smallest even I with rho(I) > S  (within range)
    Imiss = None
    for I in range(0, IMAX+1, 2):
        if rho[I] > S: Imiss = I; break
    Dstar_rho = S*S - S - Imiss
    # from spectra:
    miss = [Q for Q in range(S, S*S+1) if (Q-S)%2==0 and not (m>>Q)&1]
    Dstar_true = S*S - miss[0] if miss else 0
    mm = math.isqrt(2*S)
    if mm*mm < 2*S: mm += 1
    lb = 2*mm*S - 2*mm*mm + 2*mm - 4
    print(f"   {S:4d} {Dstar_true:8d} {2*math.sqrt(2)*S**1.5:12.0f} "
          f"{Dstar_true/S**1.5:6.3f}   {lb:8d}  (rho-based D* = {Dstar_rho})")
