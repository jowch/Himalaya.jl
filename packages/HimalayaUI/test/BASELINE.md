# Test-suite baseline (pre-perf-refactor)

Captured: 2026-06-19, perf branch based on main `0a7d44d` (SQLite_jll pinned 3.51.2).

- Total Pass: **2145**   Fail: 0   Error: 0   Broken: 0
- Wall time (`real`): **268.65 s** (~4m29s); test execution proper: 4m18.6s
- Command: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'`

Every milestone gate must show Pass == 2145 (and Fail==Error==0), modulo the
assertions the M0.4 equivalence harness adds.

## ⚠️ Discrepancy with the investigation's ~1247s / ~20-min premise

The original investigation measured ~1247s. This uncontended baseline is **268s —
a 5× difference**. The investigation run was almost certainly **CPU-contended**:
during this session `ps` showed ~5 concurrent julia processes (PlutoMCP eval
servers, Malt workers, plus another agent running the full suite). Solo, the suite
is ~4.5 min, not ~20.

Implication for scope: M1 (in-process dispatch) remains clearly worthwhile, but the
ROI of M2 (template DB) and especially M3 (parallel process buckets — 5× Julia
startup+compile overhead) is much weaker at a 4.5-min baseline than at 20 min.
**Re-measure after M1 before committing to M2/M3.**
