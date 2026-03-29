## Resubmission

This is a resubmission (version 0.1.1). The following issues identified
in the initial submission have been addressed:

* Added `Language: en-GB` to DESCRIPTION to correctly handle British
  spellings (e.g., 'modelled', 'visualise').
* Expanded abbreviations in DESCRIPTION: RCP, MRCT, RMST, and MHLW are
  now defined in full at the end of the Description field.
* Fixed invalid URL in DESCRIPTION (BugReports and URL fields now point
  to the correct GitHub repository: gosukehommaEX/SingleArmMRCT).
* Wrapped examples for `plot_rcp1armRMST()` and
  `plot_rcp1armMilestoneSurvival()` in `\donttest{}` to avoid exceeding
  the 10-second example run time limit.

## R CMD check results

Checked as SingleArmMRCT 0.1.1.

Duration: ~2m 30s

0 errors | 0 warnings | 0 notes

Tested on:
- Windows 11, R 4.5.0 (local)

## Test results

277 tests, 0 failures.

## Spell check

No spelling errors found.
