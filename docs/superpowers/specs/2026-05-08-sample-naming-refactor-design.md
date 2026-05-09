# Sample Naming Refactor — Design Spec

**Issue:** #88 (supersedes #83)
**Status:** Approved 2026-05-08
**Worktree:** `.claude/worktrees/sample-naming-refactor` (branch `sample-naming-refactor`)
**Prod-DB safety snapshot:** `/opt/Himalaya.jl/data/himalaya.db.pre-issue88.20260509T012719Z` (139 samples / 678 exposures / 421 auto_peaks; integrity OK)

## 1. Problem

The current `samples` table conflates two roles in the `name` column:

- it is the natural-key upsert target during reingest (`cli.jl:158`),
- it is user-editable via `PATCH /api/samples/:id` (`routes_samples.jl:24`),
- it is the frontend's primary display string.

A UI rename desyncs from the manifest's natural key, breaking the next reingest's upsert (silently creating a duplicate sample). The companion `label` column is similarly muddled — manifest-driven, not user-editable, and redundant in practice.

Issue #83 surfaced one symptom (silent merges on duplicate names) but the deeper issue is wrong column semantics. This spec lands the full fix and is the prerequisite for #89 (slug-based permalinks).

## 2. New data model

| Column | Purpose | Source | Mutability | In URL slug? |
|---|---|---|---|---|
| `samples.name` | Stable scientific identifier | Manifest column 2 (was `label` position) | Set at ingest; never UI-mutable | Yes |
| `samples.display_name` | Friendly user-facing label | Manifest column 3 (was `name` position); user-editable; never clobbered by reingest | User-editable | No |

`samples.label` is dropped.

- `samples.name` is convention-enforced: matches `^[A-Za-z0-9._-]+$`, non-empty, unique within experiment.
- `samples.display_name` is free-form except control characters and tabs/newlines (whitespace and `+`/`&`/etc. are allowed — display strings need to render naturally). Leading/trailing whitespace is stripped at the PATCH-route boundary (`routes_samples.jl`) before the UPDATE so the stored value is always trimmed.

**Why swap which manifest column feeds which DB column.** Conventionally column 2 is a short identifier (`JC001`) and column 3 is descriptive text (`DOPC + cholesterol`). The original schema reversed these. Aligning column 2 = identifier and column 3 = display matches operator expectations and what the data tends to look like.

## 3. Schema migration

Lives in `migrate_schema!` in `db.jl`, alongside the existing additive ALTER list. SQLite 3.35+ supports `DROP COLUMN`; the bundled `SQLite.jl` is current. No FK refs target `samples.label`, so DROP is safe in place — no rename-create-copy-drop pattern needed.

```sql
-- Each statement wrapped individually in try/swallow so a partially-applied
-- legacy DB heals on next open_db.
ALTER TABLE samples ADD COLUMN display_name TEXT;
UPDATE samples SET display_name = name,
                   name         = COALESCE(NULLIF(label, ''), name);
ALTER TABLE samples DROP COLUMN label;
CREATE UNIQUE INDEX IF NOT EXISTS samples_unique_name ON samples(experiment_id, name);
```

**SQLite UPDATE atomicity.** SQLite evaluates all RHS expressions against the OLD row before applying any assignment, so the simultaneous `display_name = name, name = COALESCE(...)` is safe — `display_name` gets the original `name` value, not the rewritten one.

**Column order on fresh vs migrated DBs.** A fresh `create_schema!` will define `samples` with columns in canonical order: `(id, experiment_id, name, display_name, notes)`. SQLite `ALTER TABLE ADD COLUMN` appends, and `DROP COLUMN` does not reorder, so a legacy DB ends up with `(id, experiment_id, name, notes, display_name)`. The frozen-shape test in `cache-shape.test.ts` compares against a `Set`, so this divergence is invisible to the contract. Code that does `SELECT *` is also order-insensitive (rows are accessed via `getproperty`, not indices). The cosmetic mismatch is acceptable; the alternative (full table rebuild for canonical order) costs FK-heal complexity for no functional gain.

**Pre-existing duplicates (warn-and-suffix).** Legacy data may already contain duplicate `(experiment_id, name)` after the COALESCE swap — the missing UNIQUE constraint never blocked it. The migration scans for duplicates **before** creating the unique index and renames the second-and-later occurrences to `<name>-2`, `<name>-3`, …, emitting a `@warn` per renamed sample. This survives every legacy DB without operator intervention and matches the issue's "warn-and-accept" stance for non-conformant convention names.

**Idempotent on re-run.** Each ALTER swallows "duplicate column" / "no such column" errors. The duplicate-suffix pass is a no-op once names are unique. The CREATE UNIQUE INDEX uses `IF NOT EXISTS`.

**Convention compliance is a separate warning pass.** After migration, scan all `samples.name` and emit per-row warnings for names that don't match `^[A-Za-z0-9._-]+$`:

```
WARNING: 12 samples have names that don't match the new naming convention.
They will work but require URL percent-encoding (#89) and will fail
validation if reingested. Offenders: [JC001 alpha, lipid 7, ...]
Fix by editing the manifest's identifier column and running reingest.
```

## 4. Validation module

New file `packages/HimalayaUI/src/validate.jl`. Pure-function, no DB.

```julia
struct ManifestViolation
    kind::Symbol            # :empty_name | :bad_name_chars | :duplicate_name |
                            # :duplicate_filename_in_sample | :overlapping_filenames
    sample_index::Int       # 1-based row in the manifest
    sample_name::String     # may be "" for :empty_name
    detail::String
end

validate_manifest(samples::Vector{ManifestSample}) -> Vector{ManifestViolation}
```

Rules:

1. `name` non-empty
2. `name` matches `^[A-Za-z0-9._-]+$`
3. `name` unique within manifest
4. Filename ranges within a sample have no internal duplicates after expansion
5. Filename ranges across samples don't overlap

Collects all violations (no fail-fast), so the operator sees every fix needed in one pass.

**Surface points:**

- **CLI `init` / `reingest`:** print all violations, indented under a header, exit non-zero, no DB write at all. Validation runs **before** `create_experiment!` / `UPDATE experiments` (currently `parse_manifest` runs after the experiments UPDATE in `_reingest_inner!`).

- **HTTP `POST /api/experiments/:id/reingest`:** new exception type:

  ```julia
  struct ManifestValidationError <: Exception
      violations::Vector{ManifestViolation}
  end
  ```

  Thrown inside the SQLite transaction so a violation rolls back any partial state. Route handler distinguishes:

  ```julia
  try
      res = reingest!(...)
      return 200 ...
  catch e
      if e isa ManifestValidationError
          return 400, JSON3.write(Dict(:error => "manifest_invalid",
                                       :violations => violations_to_json(e.violations)))
      end
      return 500, ...   # genuine internal failures
  end
  ```

- **`cli_init_with_db!` wrapped in `SQLite.transaction`** as part of this work. Currently it isn't; partial inits leak rows on failure (#83 footgun #4). `_reingest_inner!` is already wrapped.

## 5. `experiment.toml` cutover + `migrate-toml` helper

**Hard cutover.** Old keys (`[manifest].label`, with old meaning, and `[manifest].name`, with old meaning of column 3 = identifier) are rejected at parse time:

```
Error: deprecated key '[manifest].label' in /path/to/exp/experiment.toml.
The manifest column meanings were swapped: column 2 is now `name`
(stable identifier), column 3 is now `display_name` (user-friendly label).
Run `himalaya migrate-toml /path/to/exp` to upgrade automatically.
```

Detection: presence of either `[manifest].label` key triggers rejection. Post-migration files have only `name` and `display_name`.

**TOML schema** (in `simple.toml`):

```diff
 [manifest]
 sample_id      = 1
-label          = 2
-name           = 3
+name           = 2
+display_name   = 3
 filenames      = 9
 notes_sample   = 10
 notes_exposure = 11
```

**Parser changes (`config.jl`):**

- `ExperimentConfig` struct: `col_label` → `col_display_name`, both still `Union{Int,String}`.
- `_build_config`: read `[manifest].name` and `[manifest].display_name`. Reject if `[manifest].label` is present.
- `config_to_toml`: emits the new keys.

**Manifest parser (`manifest.jl`):**

- `ManifestSample` struct: `label` → `display_name`.
- `parse_manifest` reads `cfg.col_name` and `cfg.col_display_name` instead of `cfg.col_label` and `cfg.col_name`. **The column-2 cell now flows to `ms.name` (was `ms.label`); the column-3 cell flows to `ms.display_name` (was `ms.name`).** This is a behavior swap, not just a rename.

**`migrate-toml` CLI command** (new `cli_migrate_toml` in `cli.jl`):

- Signature: `himalaya migrate-toml <experiment-dir>`. Refuses without an explicit dir (no globbing across all experiments).
- Reads `<dir>/experiment.toml` raw (string-level, not `TOML.parsefile`) to preserve comments and key order.
- Replaces:
  - the line containing `label` (with original meaning) → `name = N` (where N was the old `label`'s column).
  - the line containing `name` (with original meaning) → `display_name = M` (where M was the old `name`'s column).
- Writes atomically (`mv tmp.toml experiment.toml`).
- Idempotent: if the file already has `display_name` and no `label`, prints `"already migrated"` and exits 0.
- Fails clearly if the file has neither (malformed) or both old + new (manual edit needed).

**Read-only contract.** `migrate-toml` is the second command (alongside `config new --dir`) that writes inside an experiment dir. Both are explicit, operator-invoked one-shots — fits the existing carve-out. The directory-snapshot regression in `test_pipeline.jl` only covers `init` / `analyze` / `reingest` / `serve`, so it stays green.

## 6. Backend route + DB-helper changes

**`PATCH /api/samples/:id`** (`routes_samples.jl:24`):

- Allowlist: `(:name, :notes)` → `(:display_name, :notes)`. `name` becomes immutable via this endpoint.
- `update_sample` event payload now carries `display_name` instead of `name`. `applyRemoteToCache` already spreads the patched fields onto the cached sample, so this is a free swap on the cache side. The SSE contract test (`sseEventPayload.contract.test.ts`) needs the field renamed.

**`PATCH /api/experiments/:id`** (`routes_experiments.jl:42`):

- Removes `:name` from the allowlist. With no other mutable fields today, the route becomes a placeholder that returns 400 for every request. **Decision: keep the route** (defensive surface for future fields). The frontend `updateExperiment` fetcher in `api.ts:96` is typed `Partial<Pick<Experiment, "name" | "data_dir" | ...>>` — drop `name` from the type so callers can't construct a payload that would 400.

**`POST /api/experiments/:id/reingest`** (`routes_experiments.jl:108`):

- Replace catch-all `catch e → 500` with the layered try from §4: `ManifestValidationError → 400 + structured violations`, everything else → 500.

**`create_sample!`** (`db.jl:1119`):

- Drop `label` kwarg, add `display_name` kwarg.
- Call sites:
  - `cli.jl:64` (init insert)
  - `cli.jl:164` (reingest insert path)
  - `cli.jl:170` (reingest update path: `UPDATE samples SET label = ?, notes = ?` → `UPDATE samples SET display_name = ?, notes = ?`)
  - `cli.jl:335,338` (println uses `sample.label` for log messages → `sample.name`)
- `get_samples` and `_experiment_row_to_json` need no changes — they `SELECT *` and the new column shape flows through.

**Reingest upsert key.** Continues to match on `(experiment_id, name)`, but `name` now carries the stable identifier. A sample's identity is preserved across reingests as long as the operator doesn't change the identifier column in the manifest — exactly the contract the issue argues for.

## 7. Frontend changes

**Type definition (`api.ts:26`):**

```diff
 export interface Sample {
   id: number;
   experiment_id: number;
-  label: string | null;
   name: string | null;
+  display_name: string | null;
   notes: string | null;
   tags: SampleTag[];
 }
```

**Fetcher signature (`api.ts:103`):**

```diff
-export const updateSample = (id, patch: { name?: string; notes?: string }, opts?) =>
+export const updateSample = (id, patch: { display_name?: string; notes?: string }, opts?) =>
   request<Sample>("PATCH", `/api/samples/${id}`, patch, opts);
```

**Display priority — single helper.** Today `s.name ?? s.label` is duplicated in 8+ places. Replace with one helper to make the rename auditable and prevent drift:

```ts
// lib/sample/displayName.ts
import type { Sample } from "../../api";
export const sampleDisplayName = (s: Sample): string =>
  s.display_name ?? s.name ?? `Sample #${s.id}`;
```

Call sites updated to import this:

- `queries.ts:208` (comparison label fallback)
- `lib/comparison/labels.ts:38`
- `components/SamplePickerRow.tsx:18`
- `components/NavModal.tsx:127` (primary)
- `components/ExposureListRow.tsx:57`
- `components/ComparisonPickerBody.tsx:119`
- `components/MentionPicker.tsx:46,116`
- `components/MentionChip.tsx:49`
- `components/SampleMetadataCard.tsx:65,70,73,74`

**NavModal secondary line (`NavModal.tsx:128`):** purpose flips. Today shows `label` when `name` differs (because `name` was descriptive). After: show `name` (identifier) when `display_name` exists and differs:

```ts
secondary: s.name && s.display_name && s.name !== s.display_name ? s.name : ""
```

**Search field flip:**
- `NavModal.tsx:109`: `(s.label ?? "")` → `(s.display_name ?? "")`
- `ComparisonPickerBody.tsx:119`: same.

**Frozen-shape contract test (`test/queue/cache-shape.test.ts:53`):**

```diff
-  "id", "experiment_id", "label", "name", "notes", "tags",
+  "id", "experiment_id", "name", "display_name", "notes", "tags",
```

This is the pin that catches drift between SSE payload, cache row, and `applyRemoteToCache` reads. Must change in lockstep with the backend payload rename and the Sample type.

**Test fixtures.** Every fixture builder that constructs a Sample needs `label` removed and `display_name` added. The `"falls back to sample.label when sample.name is missing"` test in `ExposureListRow.test.tsx:56` becomes `"falls back to sample.name when sample.display_name is missing"`.

**No mutation-queue mutator changes needed.** `update_sample` isn't a queue mutator — it goes through the trivial path. `applyRemoteToCache` shape-matches by spreading the payload onto the cached row, so renaming the field within the payload propagates without a registry change.

**No SSE event-kind change.** `update_sample` keeps its kind name; only the payload shape rotates.

## 8. Build sequence

One PR, commit-by-commit, each commit independently green. Order matters because the frontend frozen-shape test ties Sample type, payload, and cache row together.

1. **Backend: schema migration** — `migrate_samples_naming!` in `db.jl` running the four-statement SQL plus the warn-and-suffix duplicate-resolver. Tests in `test_db.jl` for: fresh DB → new shape; legacy DB with mixed label/name → migrated correctly; legacy DB with duplicate names → suffix-renamed with warning; idempotent on second run. The 139-row prod backup is the load-bearing fixture for the duplicate-handling case (sanitized subset can be vendored if needed).
2. **Backend: validation module** — `validate.jl` + `ManifestViolation` + `validate_manifest`. Pure-function tests in new `test_validate.jl`. No call-site changes yet.
3. **Backend: parser cutover** — `config.jl` (struct rename + reject old keys), `manifest.jl` (`ManifestSample.label → display_name`, column-meaning swap in `parse_manifest`), `simple.toml` (key rename). Updates to `test_config.jl` + `test_manifest.jl`.
4. **Backend: migrate-toml CLI** — new `cli_migrate_toml`, dispatch in `main()`, new `test_migrate_toml.jl`.
5. **Backend: route + pipeline cutover** — `create_sample!` kwarg swap, `cli.jl` init/reingest call sites, `routes_samples.jl` allowlist, `routes_experiments.jl` allowlist + 400/500 split, `validate_manifest` integrated into both write paths inside transactions. Updates to `test_routes_samples.jl`, `test_routes_experiments.jl`, `test_pipeline.jl`. **`cli_init_with_db!` wrapped in `SQLite.transaction`** as part of this commit (closes #83 footgun #4).
6. **Frontend: type + helper + atomic rename** — Sample type rename, new `lib/sample/displayName.ts`, all call sites switched in one commit. Updates to `cache-shape.test.ts`, fixture builders, `ExposureListRow` test rename. `tsc --noEmit + vite build + npm test` all green.
7. **Docs** — update `docs/experiment-config.md` with the new `[manifest]` keys and the deprecation message; update `CLAUDE.md` "samples" gotchas if any are now stale; update root `AGENTS.md` if the schema is duplicated there. Verify `docs/superpowers/specs/2026-04-22-himalaya-web-app-design.md` hasn't gone stale on the `samples` shape.

## 9. Testing strategy

Per-commit unit/integration tests are listed inline in §8. Beyond that:

- **Six-layer audit for `update_sample`** per `docs/contract-testing.md`: walk the rename through `route emit (routes_samples.jl) → SSE payload (sseEventPayload.contract.test.ts) → applyRemoteToCache merge → cache row (cache-shape.test.ts) → onMutate / onSuccess`. Each layer must agree the field is `display_name`. The contract tests already exist for it; they update in step 6.
- **Backend Julia full suite** (`Pkg.test("HimalayaUI")`, ~5–10 min). Per CLAUDE.md, capture once and grep the file.
- **Frontend Vitest full suite + Playwright mocked** (`npm test && npm run e2e`).
- **Live smoke against a copy of the prod-backup DB**: copy the backup to a scratch location, run `migrate-toml` against an existing `experiment.toml` copy, do a round-trip `init` / `reingest` against a copy DB, hit the UI, watch for stale-cache or rendering regressions. **Never touch the prod DB or prod experiment dirs.**

## 10. Out of scope (deferred per the issue)

- Schema-level `UNIQUE(name)` on `experiments` (#83 stretch).
- Schema-level `UNIQUE(sample_id, filename)` on `exposures` (#83 stretch).
- Permalink URL scheme (#89 — the consumer of this work).
- Slug column on `comparisons` (`/compare/<id>/<title-slug>` Stack Overflow style).
- Mention-chip click-through navigation.
- Experiment-name editable/stable split (parallel of this work for experiments — smaller because experiments are renamed less often).

## 11. Acceptance criteria

- [ ] `samples` table has `(id, experiment_id, name, display_name, notes)` with `UNIQUE(experiment_id, name)`.
- [ ] Legacy DB with `label`/`name` migrates correctly; integrity check passes; row counts unchanged.
- [ ] Legacy DB with duplicate `(experiment_id, name)` survives migration via suffix-rename with warnings.
- [ ] `experiment.toml` parser rejects old `[manifest].label` key with a clear error pointing at `migrate-toml`.
- [ ] `himalaya migrate-toml <dir>` rewrites the file correctly, is idempotent, fails clearly on malformed input.
- [ ] `validate_manifest` returns all five violation kinds; CLI prints all violations and refuses to write; HTTP returns 400 + structured JSON.
- [ ] `cli_init_with_db!` is wrapped in `SQLite.transaction`; partial inits roll back fully.
- [ ] `PATCH /api/samples/:id` rejects `name`; accepts `display_name` and `notes`.
- [ ] `PATCH /api/experiments/:id` rejects `name`; returns 400 if no other fields are supplied.
- [ ] Frontend renders `sampleDisplayName(s)` everywhere; no remaining `s.label` references.
- [ ] All contract tests green (cache-shape, SSE payload, route response shapes, idempotency replay).
- [ ] Full backend Julia suite green; full frontend Vitest + Playwright mocked green.
- [ ] Live smoke against a prod-backup copy completes without UI regression.
