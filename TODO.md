* Check that all functions are in line with manuscript math, e.g., re parity and counter parity corrections (delta correction in muilwijk etc)
* Maybe add the muilwijk bound to strait's functions? Currently only does the sharp ones, but maybe folk want the better known ones. Note in docs that it's more often called bhatia-davis

- **Revisit the phrasing of the inconsistency decisions.** `consistent` is a single Boolean over tests with very different epistemic status, and the wording should probably reflect that:
    - `in_scale_range` and `bounds` failures are *arithmetic proofs of impossibility* — no dataset with those summary statistics exists.
    - `grim` and `grimmer` are deferred third-party verdicts with documented false-positive cases, which is why a GRIMMER-only failure already earns a caveat in `note`.
    - `feasibility` is now a residual, meaning "no sample exists for a reason none of the named tests accounts for" — accurate but uninformative to a reader.
    Calling all of these "inconsistent" flattens a proof and a flag into one word. Consider distinguishing *impossible* from *flagged*, deciding whether `failed_tests` tokens should be the test names (`brim`, `brimmer`, `grim`, `grimmer`) rather than the current mixture of column names and concepts (`in_scale_range`, `bounds`), and settling how the verdict should be phrased in a report a reader may act on. Worth fixing before CRAN pins the return shape.
