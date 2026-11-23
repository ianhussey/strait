# Every multiset of n values drawn from `values`, one per row. A sample's mean
# and SD do not depend on the order of its values, so a brute-force check
# needs these choose(n + k - 1, n) rows rather than all k^n ordered samples:
# 1,716 instead of 823,543 for n = 7 on a 0-6 scale.
multisets <- function(values, n) {
  k <- length(values)
  bars <- utils::combn(n + k - 1, k - 1) # stars and bars
  counts <- diff(rbind(0, bars, n + k)) - 1 # one column of counts each
  t(apply(counts, 2, \(cnt) rep(values, cnt)))
}
