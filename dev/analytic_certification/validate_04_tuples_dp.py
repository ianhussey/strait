# Validation 04: reported-tuple layer vs brute-force ground truth,
# at d = 1, 2 AND 3 reported decimals.
#
# For each design the DP enumerates the exact achievable (S, Q) states. A
# reported (mean, sd) is POSSIBLE iff some integer S in the mean's rounding
# window admits an integer Q in the sd's rounding window that is achievable
# (window semantics: closed interval on both sides = union of half_up and
# half_down, matching brimmest()'s default).
#
# The analytic pipeline (GRIM window -> window_possible) must reproduce this
# exactly. d = 3 is the hardest case: windows are so narrow they usually pin
# a single Q, stressing parity/endpoint semantics and the gap structure.
#
# Coverage per design/digits: the ENTIRE top-edge strip (every mean on the
# reporting grid, all sd within 0.35 of sd_max — where every phantom lives),
# plus random interior/bottom tuples as controls.
#
# Expected output: 0 mismatches.
# Run: python3 -u validate_04_tuples_dp.py          (~2-10 min)
import math, random, sys, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from window import window_possible, qmin_i, qmax_i

def ach_sets(W, n):
    cur = {(0, 0)}
    for _ in range(n):
        cur = {(S + v, Q + v * v) for (S, Q) in cur for v in range(W + 1)}
    out = {}
    for (S, Q) in cur:
        out.setdefault(S, set()).add(Q)
    return out

TOL = 1e-9

def s_window(mean_y, n, d):
    h = 0.5 * 10 ** -d
    return (math.ceil((mean_y - h) * n - TOL), math.floor((mean_y + h) * n + TOL))

def q_window(sd, S, n, d):
    h = 0.5 * 10 ** -d
    q_lo = (n - 1) * max(0.0, sd - h) ** 2 + S * S / n
    q_hi = (n - 1) * (sd + h) ** 2 + S * S / n
    return (math.ceil(q_lo - TOL), math.floor(q_hi + TOL))

def truth_possible(A, mean_y, sd, W, n, d):
    S_lo, S_hi = s_window(mean_y, n, d)
    for S in range(max(0, S_lo), min(n * W, S_hi) + 1):
        Q_lo, Q_hi = q_window(sd, S, n, d)
        if any(Q in A[S] for Q in range(max(Q_lo, 0), Q_hi + 1)):
            return True
    return False

def analytic_possible(mean_y, sd, W, n, d):
    S_lo, S_hi = s_window(mean_y, n, d)
    for S in range(max(0, S_lo), min(n * W, S_hi) + 1):
        Q_lo, Q_hi = q_window(sd, S, n, d)
        if window_possible(S, Q_lo, Q_hi, W, n):
            return True
    return False

def check_design(l, u, n, d, n_random=4000, seed=5):
    W = u - l
    A = ach_sets(W, n)
    random.seed(seed)
    step = 10 ** -d
    sd_cap = (W / 2) * math.sqrt(n / (n - 1))
    tuples = []
    # (i) full top-edge strip: every grid mean, sd within 0.35 of that mean's max
    n_means = round((u - l) / step)
    for i in range(n_means + 1):
        m = l + i * step
        S_lo, S_hi = s_window(m - l, n, d)
        if S_hi < 0 or S_lo > n * W or S_lo > S_hi:
            continue
        sm = max(math.sqrt((qmax_i(S, W, n) - S * S / n) / (n - 1))
                 for S in range(max(0, S_lo), min(n * W, S_hi) + 1))
        k_hi = math.floor(sm / step)
        for k in range(max(0, k_hi - round(0.35 / step)), k_hi + 2):
            tuples.append((round(m, d), round(k * step, d)))
    # (ii) random controls anywhere
    for _ in range(n_random):
        m = round(l + random.random() * (u - l), d)
        s = round(random.random() * sd_cap, d)
        tuples.append((m, s))
    bad = 0
    n_phantom = 0
    t0 = time.time()
    for (m, s) in tuples:
        tr = truth_possible(A, m - l, s, W, n, d)
        pr = analytic_possible(m - l, s, W, n, d)
        if tr != pr:
            bad += 1
            print("MISMATCH", l, u, n, d, m, s, tr, pr)
        if not tr:
            # count phantoms among the impossible: those passing the parity screen
            S_lo, S_hi = s_window(m - l, n, d)
            for S in range(max(0, S_lo), min(n * W, S_hi) + 1):
                Q_lo, Q_hi = q_window(s, S, n, d)
                Q_lo = max(Q_lo, qmin_i(S, n))
                Q_hi = min(Q_hi, qmax_i(S, W, n))
                if Q_lo > Q_hi:
                    continue
                if (Q_lo - S) % 2:
                    Q_lo += 1
                if Q_lo <= Q_hi:
                    n_phantom += 1
                    break
    print(f"[04] {l}-{u} n={n} d={d}: {len(tuples)} tuples, {bad} mismatches, "
          f"{n_phantom} phantoms ({time.time()-t0:.0f}s)")
    return bad

if __name__ == "__main__":
    fail = 0
    fail += check_design(1, 5, 9, 1)
    fail += check_design(1, 5, 9, 2)
    fail += check_design(1, 5, 9, 3)
    fail += check_design(1, 7, 25, 2)
    fail += check_design(1, 7, 25, 3)
    fail += check_design(0, 10, 15, 3)
    sys.exit(1 if fail else 0)
