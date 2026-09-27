import math
# exact rho DP up to 40000
IMAX = 40000
INF = 10**9
rho = [INF]*(IMAX+1); rho[0] = 0
p = 2; pr = []
while p*(p-1) <= IMAX: pr.append(p); p += 1
for I in range(2, IMAX+1, 2):
    best = INF
    for q in pr:
        c = q*(q-1)
        if c > I: break
        w = rho[I-c] + q
        if w < best: best = w
    rho[I] = best

# (A) adversarial construction: for S, m = floor(sqrt(2S)), u = S-m,
#     I* = u(u-1) + (m-1)(m-2) + 2.  Claim rho(I*) > S  (so I* is missing at S),
#     giving D*(S) >= 2mS - 2m^2 + 2m - 4.
print("[A] adversarial lower bound: rho(I*) vs S")
bad = 0
for S in range(32, 190):
    m = math.isqrt(2*S)
    u = S - m
    Istar = u*(u-1) + (m-1)*(m-2) + 2
    if Istar > IMAX: continue
    ok = rho[Istar] > S
    if not ok:
        bad += 1
        print("  FAIL", S, m, Istar, rho[Istar])
print(f"  checked S in [32,190) (where I* <= {IMAX}): {bad} failures")

# (B) greedy bound at large I: B(I) = u1 + isqrt(2 u1) + 2*(2u1)^(1/4) + 8
def greedy_weight(I):
    w = 0
    while I > 0:
        u = (1 + math.isqrt(1 + 4*I)) // 2
        while u*(u-1) > I: u -= 1
        while (u+1)*u <= I: u += 1
        w += u; I -= u*(u-1)
    return w
import random
random.seed(4)
worst = 0
for _ in range(200000):
    I = 2*random.randint(1, 5*10**11)
    g = greedy_weight(I)
    u1 = (1 + math.isqrt(1 + 4*I)) // 2
    B = u1 + math.isqrt(2*u1) + 2*int((2*u1)**0.25) + 8
    if g > B: worst += 1; print("  greedy bound FAIL", I, g, B)
print(f"[B] greedy bound to I ~ 1e12: {worst} violations in 200k random trials")

# (C) explicit sufficient condition: every even I with
#     S^2 - S - I >= 2*sqrt(2)*S^1.5 + 5*S^1.25 + 24*S  attainable (rho <= S)?
print("[C] explicit upper-bound formula check")
bad = 0
for S in range(20, 190):
    D_needed = 2*math.sqrt(2)*S**1.5 + 5*S**1.25 + 24*S
    for I in range(0, min(IMAX, S*S - S) + 1, 2):
        if S*S - S - I >= D_needed and rho[I] > S:
            bad += 1
            if bad < 5: print("  FAIL", S, I, rho[I])
print(f"  {bad} violations for S < 190")
