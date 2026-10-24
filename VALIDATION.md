# LiDAR 0.1.1 validation

## Fix validation

* Version 0.1.1 source build and installation succeeded.
* The installed-package core tests also passed with two dnet-dependent skips.

* Core test suite: 33 expectations passed, 0 failures, 0 warnings, 2 optional
  dnet integration tests skipped because the full dnet package is unavailable.
* Separately downloaded and ran the upstream dRWR serial function from
  hfang-bristol/dnet against both connected and original graph fixtures.
* Reproduced zero cell activity for the original 3-neighbor graph (two components).
  Its normalized activity previously produced NaN. The 5-neighbor graph was connected,
  but the constant receptor labels tied candidate selection at zero and selected the
  first graph, matching the original implementation.
* Corrected zero-sum normalization preserves zeros and emits a warning.
* Positive-sum normalization is numerically unchanged in the regression test.
* The connected positive-signal fixture produced positive finite activity.
* This upstream-function check does not replace installation of the complete dnet
  package or the end-to-end Seurat/CelliD biological reference workflow.

## Historical baseline

Version 0.1.0 completed R CMD check with 0 errors, 0 warnings and 2 NOTEs
(optional dependencies unavailable; utils import declaration). Its RWR integration
test skipped locally and subsequently failed in the user's dnet-enabled environment.
The baseline check log is provided as historical evidence, not a full 0.1.1 check.

## Still pending

* Full R CMD check for 0.1.1.
* Full dnet-package integration tests in the user's environment.
* QS to Seurat to MCA/RWR to biological reference comparison.
* Real maintainer details, license selection and remote GitHub hosting.
