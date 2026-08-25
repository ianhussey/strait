# Reference analytic certifier for integer achievability of (S, Q) on a bounded scale,
# plus exhaustive validation against brute-force DP.
#
#   member_top(S, Q, W, n): is there a multiset of n integers in [0, W] with
#                           sum S and sum of squares Q?
#
# Method: band recursion on the largest part with parity/sandwich pruning and
# memoization. Exact and unconditionally correct; fast because spectral gaps exist
# only in the self-similar top zone (see analytic-certification.md).
#
# Validation (run this file):  0 mismatches over all W <= 8, 2 <= n <= 12, all S, Q
# (859,672 checks), ~1s total.
import math, sys
sys.setrecursionlimit(100000)

def qmin_i(S, n):
    """Smallest sum of squares of n integers >= 0 summing to S (clustered config)."""
    m, r = divmod(S, n)
    return (n - r) * m * m + r * (m + 1) * (m + 1)

def qmax_i(S, cap, n):
    """Largest sum of squares of n integers in [0, cap] summing to S (Structure-S),
    or None if S > n*cap."""
    if cap == 0:
        return 0 if S == 0 else None
    a, b = divmod(S, cap)
    if a > n or (a == n and b > 0):
        return None
    return a * cap * cap + b * b

def member(S, Q, cap, n, memo):
    """Exact membership: Q in { sum y^2 : y_1..y_n in [0,cap], sum y = S }."""
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
    """Top-level entry: reflect to S <= n*W/2, then decide."""
    if 2 * S > n * W:
        Q = Q + n * W * W - 2 * W * S
        S = n * W - S
    return member(S, Q, W, n, {})

if __name__ == "__main__":
    def ach_sets(W, n):
        cur = {(0, 0)}
        for _ in range(n):
            cur = {(S + v, Q + v * v) for (S, Q) in cur for v in range(W + 1)}
        out = {}
        for (S, Q) in cur:
            out.setdefault(S, set()).add(Q)
        return out

    import time
    bad = total = 0
    t0 = time.time()
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
    print(f"checked {total}, mismatches {bad}, {time.time()-t0:.1f}s")
