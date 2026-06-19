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

## M1 result (in-process dispatch migration complete)

- Total Pass: **2170** (baseline 2145 + 25 equivalence-harness assertions), Fail 0, Error 0.
- Wall time (`real`): **186.90 s** (~3m07s); test execution 3m02.4s.
- vs baseline 268.65s → **~30% faster** (uncontended). 21 route-test files migrated from
  per-test Oxygen server boots to in-process `internalrequest` dispatch; ~145 boots eliminated.
- Beyond the raw 30%: in-process dispatch removes the server-boot/port-bind step entirely, so
  the suite no longer balloons under CPU contention the way the wire version did (the original
  ~1247s contended measurement). That robustness is arguably the bigger win.

Go/no-go threshold was <600s; we are at 187s. The <3min goal is essentially met by M1 alone.
M2/M3 ROI at this scale is marginal — see scope note above; deferred pending re-measure decision.

## M2 result (template-DB clone for fresh fixtures)

- Measured per-fixture cost: `open_db` 37.4ms → `open_prepared_clone` 5.1ms (~32ms/call saved).
- Swapped 153 fresh-fixture sites across 25 files to `open_prepared_clone`. EXCLUDED (left on
  `open_db`): `test_db.jl` + `test_pipeline.jl` (legacy-migration testsets — `open_db` is the
  unit-under-test), the equivalence harness, `test_concurrent_writes.jl`, and `test_fast_skip.jl`
  (runs before the helper is defined in include order).
- Full suite: **2170/2170, 192.6s** — i.e. ~same as M1's 187s. The ~5.8s theoretical saving is
  within full-suite run-to-run variance, so M2 shows no standalone wall-clock win at this scale.
  (Its value is real but small; it compounds slightly under M3's smaller parallel buckets.)
