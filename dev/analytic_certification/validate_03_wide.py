# Validation 03: wide scales and awkward n (sum-score-style designs).
#
# (a) Brute-force DP ground truth at W = 300, n = 4 (the largest exactly
#     tractable wide DP), three-way against both analytic algorithms.
# (b) Cross-validation of the two algorithmically independent certifiers
#     (largest-part recursion, certifier.py vs demotion decomposition,
#     window.py) at W in {173, 300, 517}, n in {29, 37, 61, 83, 217},
#     oversampling the top zone (where all gaps live) plus random controls,
#     point queries and windows.
#
# Expected output: 0 mismatches.
# Run: python3 -u validate_03_wide.py               (~10-30 min)
import random, sys, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from certifier import member_top, qmin_i, qmax_i
from window import window_possible

def dp_masks(W, n):
    d = {0: 1}
    for _ in range(n):
        nd = {}
        for S, m in d.items():
            for v in range(W + 1):
                nd[S + v] = nd.get(S + v, 0) | (m << (v * v))
        d = nd
    return d

if __name__ == "__main__":
    fail = 0
    # (a) DP ground truth
    t0 = time.time()
    random.seed(11)
    W, n = 300, 4
    d = dp_masks(W, n)
    bad = tot = 0
    for S in random.sample(range(0, n * W + 1), 120):
        m = d[S]
        lo, hi = qmin_i(S, n), qmax_i(S, W, n)
        qs = [hi - k for k in range(0, min(hi - lo, 12000), 2)][:400] + \
             [random.randint(lo, hi) for _ in range(200)] + [lo, lo + 2, (lo + hi) // 2]
        for Q in qs:
            if Q < lo or Q > hi:
                continue
            truth = (m >> Q) & 1 == 1
            p1 = member_top(S, Q, W, n)
            p2 = window_possible(S, Q, Q, W, n)
            tot += 1
            if not (truth == p1 == p2):
                bad += 1
                print("MM", W, n, S, Q, truth, p1, p2)
    fail += bad
    print(f"[03a] DP ground truth W={W} n={n}: {tot} checks, {bad} mismatches "
          f"({time.time()-t0:.0f}s)")

    # (b) cross-validation on wide scales, awkward n
    random.seed(12)
    for (W, n) in [(300, 37), (300, 83), (173, 61), (517, 29), (300, 217)]:
        bad = tot = 0
        t0 = time.time()
        Ss = random.sample(range(1, n * W // 2 + 1), 60) + \
             [random.randint(1, 3 * W) for _ in range(30)] + \
             [k * W + b for k in (1, 2, 3) for b in (0, 1, W // 2, W - 1)]
        for S in Ss:
            lo, hi = qmin_i(S, n), qmax_i(S, W, n)
            if hi is None:
                continue
            qs = [hi - k for k in range(0, min(hi - lo, 8000), 2)][:250] + \
                 [random.randint(lo, hi) for _ in range(80)]
            for Q in qs:
                if Q < lo or Q > hi:
                    continue
                p1 = member_top(S, Q, W, n)
                p2 = window_possible(S, Q, Q, W, n)
                tot += 1
                if p1 != p2:
                    bad += 1
                    print("XMM", W, n, S, Q, p1, p2)
            for _ in range(40):                       # windows near the top zone
                x = hi - random.randint(0, 6000)
                y = x + random.randint(0, 300)
                p2 = window_possible(S, x, min(y, hi), W, n)
                p1 = any(member_top(S, Q, W, n)
                         for Q in range(max(x, lo) + ((max(x, lo) - S) % 2),
                                        min(y, hi) + 1, 2))
                tot += 1
                if p1 != p2:
                    bad += 1
                    print("WMM", W, n, S, x, y, p1, p2)
        fail += bad
        print(f"[03b] cross-val W={W} n={n}: {tot} checks, {bad} mismatches "
              f"({time.time()-t0:.0f}s)")
    sys.exit(1 if fail else 0)
