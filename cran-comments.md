## R CMD check results

Duration: ~2m 30s

0 errors | 0 warnings | 0 notes

Tested on:
- Windows 11, R 4.5.0 (local)

## Test results

277 tests, 0 failures.

## Spell check

No spelling errors found.

## Notes to CRAN reviewers

This is the first submission of the SingleArmMRCT package.

The package provides functions to calculate and visualise the Regional
Consistency Probability (RCP) for single-arm multi-regional clinical trials
(MRCTs) using the Effect Retention Approach (ERA). Six endpoint types are
supported: continuous, binary, count (negative binomial), time-to-event via
hazard ratio, milestone survival, and restricted mean survival time (RMST).

The package extends classical two-arm MRCT consistency evaluation methods
(Japanese MHLW Method 1 and Method 2) to the single-arm setting, addressing
a methodological gap that is of practical relevance for oncology and rare
disease drug development.

There are no reverse dependencies as this is an initial submission.
