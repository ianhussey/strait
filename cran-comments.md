## Submission

This is a new submission.

The package was developed under the name `tides` and renamed to `strait`
before submission, to avoid a case-insensitive collision with the existing
CRAN package `Tides`. It has not previously been submitted to CRAN under
either name.

## Test environments

* local macOS 26.5 (aarch64), R 4.5.2

## R CMD check results

0 errors | 0 warnings | 2 notes.

* "New submission" — expected for a first submission. The same note lists
  `https://jamesheathers.medium.com/...` in the vignette as possibly invalid
  (status 403): Medium refuses automated requests, and the page opens in a
  browser.
* "Skipping checking HTML validation: 'tidy' doesn't look like recent enough
  HTML Tidy" — a property of the HTML Tidy binary on the local machine, not of
  the package.

## Notes for the reviewer

* The package contains no compiled code. It writes no files, and does not
  modify the user's options, graphics parameters or working directory.
* The granularity tests GRIM and GRIMMER are deferred to the `scrutiny`
  package rather than reimplemented, and reported values are unrounded with
  `roundwork`. The package imports ggplot2, roundwork, scrutiny, splines and
  tibble.
* All exported functions have running examples. Two examples in
  `sd_region_plot()` that enumerate the attainable (mean, SD) lattice are
  wrapped in `\donttest{}`; they complete in a few seconds and pass under
  `--run-donttest`.
* `Language` is `en-GB`. Domain terms, cited author surnames and possessive
  forms that the spell checker does not recognise are recorded in
  `inst/WORDLIST`, so `spelling::spell_check_package()` reports nothing.
* The reference in the Description field is given as `Authors (year) <doi:...>`.
