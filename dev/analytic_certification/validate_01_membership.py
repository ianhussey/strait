# Validation 01: exact membership certifier vs brute-force DP.
#
# Exhaustive over all designs W <= 8, 2 <= n <= 12: every S, every Q.
# Plus randomized spot checks at W = 63 (0-63 inventory), n = 6.
#
# Expected output: 0 mismatches everywhere.
# Run: python3 validate_01_membership.py            (~1 min)
import math, random, sys, time
sys.setrecursionlimit(100000)

# ---- certifier under test (self-contained copy of certifier.py core) --------
def qmin_i(S, n):
    m, r = divmod(S, n)
    return (n - r) * m * m + r * (m + 1) * (m + 1)

def qmax_i(S, cap, n):
    if cap == 0:
        return 0 if S == 0 else None
    a, b = divmod(S, cap)
    if a > n or (a == n and b > 0):
        return None
    return a * cap * cap + b * b

def member(S, Q, cap, n, memo):
    if S == 0:
        return Q == 0
    if n <= 0 or cap <= 0:
        return False
    if cap > S:
        cap = S
    key = (S, Q, cap, n)
    v = memo.get(key)
    if v is not None:
        return v
    lo = qmin_i(S, n)
    hi = qmax_i(S, cap, n)
    if hi is None or Q < lo or Q > hi or (Q - S) % 2:
        memo[key] = False
        return False
    if Q == lo or Q == hi:
        memo[key] = True
        return True
    res = False
    for p1 in range(min(cap, math.isqrt(Q)), -(-S // n) - 1, -1):
        if member(S - p1, Q - p1 * p1, p1, n - 1, memo):
            res = True
            break
    memo[key] = res
    return res

def member_top(S, Q, W, n):
    if 2 * S > n * W:
        Q = Q + n * W * W - 2 * W * S
        S = n * W - S
    return member(S, Q, W, n, {})

# ---- ground truth ------------------------------------------------------------
def ach_sets(W, n):
    cur = {(0, 0)}
    for _ in range(n):
        cur = {(S + v, Q + v * v) for (S, Q) in cur for v in range(W + 1)}
    out = {}
    for (S, Q) in cur:
        out.setdefault(S, set()).add(Q)
    return out

def ach_masks(W, n):
    d = {0: 1}
    for _ in range(n):
        nd = {}
        for S, m in d.items():
            for v in range(W + 1):
                nd[S + v] = nd.get(S + v, 0) | (m << (v * v))
        d = nd
    return d

if __name__ == "__main__":
    t0 = time.time()
    bad = total = 0
    for W in range(1, 9):
        for n in range(2, 13):
            A = ach_sets(W, n)
            for S in range(0, n * W + 1):
                qs = A[S]
                for Q in range(0, W * W * n + 1):
                    total += 1
                    if (Q in qs) != member_top(S, Q, W, n):
                        bad += 1
                        print("MISMATCH", W, n, S, Q)
    print(f"[01a] exhaustive W<=8 n<=12: {total} checks, {bad} mismatches "
          f"({time.time()-t0:.0f}s)")

    t0 = time.time()
    W, n = 63, 6
    d = ach_masks(W, n)
    random.seed(1)
    bad2 = tot2 = 0
    for S in range(0, n * W + 1, 7):
        m = d[S]
        hi = qmax_i(S, W, n)
        lo = qmin_i(S, n)
        Qs = list(range(lo, hi + 1))
        random.shuffle(Qs)
        for Q in Qs[:200]:
            tot2 += 1
            if ((m >> Q) & 1 == 1) != member_top(S, Q, W, n):
                bad2 += 1
                print("MISMATCH63", S, Q)
    print(f"[01b] W=63 n=6 randomized: {tot2} checks, {bad2} mismatches "
          f"({time.time()-t0:.0f}s)")
    sys.exit(1 if bad + bad2 else 0)
