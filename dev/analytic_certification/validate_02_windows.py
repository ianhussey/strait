# Validation 02: window-level certifier (demotion decomposition, window.py)
# vs brute-force DP.
#
# For all designs W <= 8, 2 <= n <= 12: every point window [Q, Q] inside the
# sandwich, plus random windows of width 0-12 (including out-of-range edges).
#
# Expected output: 0 mismatches.
# Run: python3 validate_02_windows.py               (~1 min)
import random, sys, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from window import window_possible  # implementation under test

def ach_sets(W, n):
    cur = {(0, 0)}
    for _ in range(n):
        cur = {(S + v, Q + v * v) for (S, Q) in cur for v in range(W + 1)}
    out = {}
    for (S, Q) in cur:
        out.setdefault(S, set()).add(Q)
    return out

if __name__ == "__main__":
    t0 = time.time()
    random.seed(3)
    bad = tot = 0
    for W in range(1, 9):
        for n in range(2, 13):
            A = ach_sets(W, n)
            for S in range(0, n * W + 1):
                qs = A[S]
                hi, lo = max(qs), min(qs)
                for Q in range(lo, hi + 1):          # all point windows
                    tot += 1
                    if (Q in qs) != window_possible(S, Q, Q, W, n):
                        bad += 1
                        print("PT MISMATCH", W, n, S, Q)
                for _ in range(20):                   # random windows
                    x = random.randint(lo - 3, hi + 3)
                    y = x + random.randint(0, 12)
                    truth = any(q in qs for q in range(max(x, 0), y + 1))
                    tot += 1
                    if truth != window_possible(S, x, y, W, n):
                        bad += 1
                        print("WIN MISMATCH", W, n, S, x, y)
    print(f"[02] windows W<=8 n<=12: {tot} checks, {bad} mismatches "
          f"({time.time()-t0:.0f}s)")
    sys.exit(1 if bad else 0)
