# Sample Naming Refactor — Design Spec

**Issue:** #88 (supersedes #83)
**Status:** Revised 2026-05-08 (post agent-review pass against `himalaya-reviewer`, `frontend-reviewer`, `queue-reviewer`)
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

**Pre-existing duplicates (warn-and-suffix).** Legacy data may already contain duplicate `(experiment_id, name)` after the COALESCE swap — the missing UNIQUE constraint never blocked it. After the COALESCE pass and **before** creating the unique index, the migration scans for duplicates ordered `BY id ASC` (oldest sample wins the bare name; deterministic across runs against the same DB) and renames the second-and-later occurrences to `<name>-2`, `<name>-3`, …, emitting a `@warn` per renamed sample.

**Atomicity.** All four statements + the duplicate-rename pass + the CREATE UNIQUE INDEX run inside a single `SQLite.transaction(db) do … end`. Without atomic scoping, a Ctrl-C between DROP COLUMN and the duplicate scan leaves a DB where COALESCE has rewritten `name` but no uniqueness exists; the next `open_db` re-runs the COALESCE pass against already-migrated data and double-suffixes correct rows. With atomic scoping, partial state is impossible.

**Migration ordering inside `migrate_schema!`.** Run **after** `migrate_pk_to_autoincrement!` and `_fix_fk_references_after_autoincrement_migration!` so the FK-heal sees a stable `samples` shape before this migration mutates columns. `samples` is in `_AUTOINCREMENT_ENTITIES` so the AUTOINCREMENT rebuild has already run by the time we touch `label`.

**Idempotent on re-run.** Each ALTER swallows "duplicate column" / "no such column" errors. The duplicate-suffix pass is a no-op once names are unique. CREATE UNIQUE INDEX uses `IF NOT EXISTS`. The migration also has a sentinel guard at the top: skip the entire helper if `display_name` already exists AND `label` does not (cheap pragma_table_info check).

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

- **`cli_init_with_db!` wrapped in `SQLite.transaction`** as part of this work, scoped to mirror `_reingest_inner!`: extract the manifest-driven create_sample!/create_exposure! insert phase into `_cli_init_inner!` and wrap that. The auto-analyze step (currently outside the transaction per the function's docstring) **stays outside** — wrapping it would mean a crash mid-analyze rolls back the whole experiment registration. The wrap covers experiment row creation + manifest validation + sample/exposure inserts only.

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
- **Regex-anchored substitutions** so neither inline comments (`label = 2  # …`) nor unrelated occurrences of "label" (in `[manifest]` headers, comments mentioning "axis label units" in `[beamline]`, `display_name` after partial migration, etc.) are misfired:
  - Match `^\s*label\s*=\s*(\S+)` within the `[manifest]` section only (track section context line-by-line) → rewrite to `name = $1` (preserving the indentation and any trailing comment).
  - Match `^\s*name\s*=\s*(\S+)` within the same section → rewrite to `display_name = $1` (also preserving formatting).
- Writes atomically (`mv tmp.toml experiment.toml`).
- Idempotent: if the `[manifest]` section already has `display_name` and no `label`, prints `"already migrated"` and exits 0. Detection uses the same regex-anchored scan over the section, not raw `occursin("label", contents)`.
- Fails clearly if the section has neither (malformed) or both old + new (manual edit needed).

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
  - `cli.jl:311` (`--sample` flag filter: `filter!(sm -> sm.label == sample_filter, samples)` → `sm.name == sample_filter` — operator passes the stable identifier)
  - `cli.jl:335,338` (println uses `sample.label` for log messages → `sample.name` — log shows the stable id)
  - `cli.jl:357,374` (`--sample/-s` flag help text: "sample label" → "sample name (identifier)")
  - `cli.jl:383` (`findfirst(sm -> sm.label == p[:sample], samples)` → `sm.name == p[:sample]`)

- `get_samples` and `_experiment_row_to_json` need no changes — they `SELECT *` and the new column shape flows through.

**`comparisons.jl::picker_samples`** (`packages/HimalayaUI/src/comparisons.jl:507`): the picker query has an explicit column list — `SELECT id, experiment_id, name, label, notes FROM samples …`. Replace `label` with `display_name`. The picker projects to a different shape than the main samples route.

**`routes_export.jl`** (entire wire-contract change):
- Line 13: test fixture builds with `label="D1", name="UX1"` — update test (covered in §9).
- Lines 56, 67, 73, 81: the JSON `Dict(:label => sm.label)` field, the CSV header `["sample_label", "sample_name", …]`, and the per-row writes `sm[:label]`, `sm[:name]`. After the rename: emit `display_name` (or drop the column outright). **Decision: replace `sample_label` with `sample_display_name`** in both the JSON dict key and the CSV header. Downstream consumers (any external pipelines parsing this CSV) will see the column header change — flag in CHANGELOG/release notes.

**`routes_experiments.jl:95`**: error-skipped log string `"$(sm.label)/$(ex.filename): $(sprint(showerror, e))"` → `"$(sm.name)/$(ex.filename): …"`.

**`update_sample` event payload + idempotent_responses migration:** the response body cached by `with_idempotency` for pre-deploy `PATCH /api/samples/:id` calls carries the OLD shape (`{label, name=descriptive}`). Post-deploy, a retried `client_op_id` would replay the old-shape body verbatim and the frontend mutator would write the wrong field into the cache. **Mitigation: `migrate_samples_naming!` runs `DELETE FROM idempotent_responses` as its final step** (the cache is a 1h TTL anyway — purging on the migration deploy aligns the wire format). Document in CHANGELOG so operators know in-flight retries from the deploy window will get fresh executions instead of cached responses.

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

**Display priority — single helper.** Today `s.name ?? s.label` is duplicated across many sites. Replace with one helper to make the rename auditable and prevent drift:

```ts
// lib/sample/displayName.ts
import type { Sample } from "../../api";
// Uses `||` not `??` so an empty-string display_name falls through rather
// than rendering as a blank tile or a leading separator (e.g. " · run.dat"
// in comparison labels). Matches the existing logic in lib/comparison/labels.ts.
export const sampleDisplayName = (s: Pick<Sample, "id" | "name" | "display_name">): string =>
  s.display_name || s.name || `Sample #${s.id}`;
```

Two design choices in the helper:
1. **`||` not `??`** — `lib/comparison/labels.ts:36–37` documents this explicitly: empty-string fields must fall through. Keeping `??` regresses the empty-string-fallback semantics.
2. **`Pick<Sample, …>` not `Sample`** — the mention pickers (`MentionPicker.tsx`, `MentionChip.tsx`) hold tagged-union rows where `sample` is partial; widening the parameter avoids casts at every picker call site.

Comprehensive call-site list (verified by grep against the worktree):

- `queries.ts:208` (comparison label fallback)
- `lib/comparison/labels.ts:38` (inline expression — keep inline, the file already documents the `||`-vs-`??` rationale; just swap `s.name || s.label` to `s.display_name || s.name`)
- `lib/comparison/labels.ts:7,9` (JSDoc fallback chain — update prose)
- `components/SamplePickerRow.tsx:18`
- `components/NavModal.tsx:109` (search needle haystack — see "Search-field UX" below for the bigger change)
- `components/NavModal.tsx:127` (primary line)
- `components/NavModal.tsx:128` (secondary line — see below)
- `components/NavModal.tsx:204` (chip-label generator)
- `components/ExposureListRow.tsx:57`
- `components/ComparisonPickerBody.tsx:119` (search needle — see below)
- `components/ComparisonPickerBody.tsx:16,22,28` (inline mock fixtures — strip `label`)
- `components/MentionPicker.tsx:46,116`
- `components/MentionChip.tsx:49`
- `components/PlotCard.tsx:83`
- `components/TitleButton.tsx:28`
- `components/SampleMetadataCard.tsx:65,70,73,74` — see "SampleMetadataCard rewrite" below; this is the rename UI flow.

**NavModal secondary line (`NavModal.tsx:128`):** purpose flips. Today shows `label` when `name` differs (because `name` was descriptive). After: show `name` (identifier) when `display_name` exists and differs:

```ts
secondary: s.name && s.display_name && s.name !== s.display_name ? s.name : ""
```

**Search-field UX (decision: search both fields).** Today `NavModal.tsx:107–110` searches **both** `s.name` and `s.label`. After the column-meaning swap, the user types either the friendly text (`display_name`) or the stable identifier (`name`) — both must match. Update the haystack to `[s.name ?? "", s.display_name ?? ""].some(h => h.toLowerCase().includes(needle))`. Same change in `ComparisonPickerBody.tsx:119`. The spec previously said "swap `label` to `display_name`" — that was wrong, it would lose identifier search.

**`SampleMetadataCard.tsx` rewrite — the rename UI.**

This card is the only existing rename surface in the app. It binds an `<input>` to `sample.name` and posts `{ name }` on blur (`SampleMetadataCard.tsx:14, 27, 45, 47, 81–88`). After the route allowlist swap, that input silently 400s. Mandatory rewrite:

- Input now binds to `sample.display_name` (controlled `useState` keyed off `sample.display_name`).
- Submit posts `updateSample(id, { display_name })`.
- The breadcrumb at lines 65–76 currently renders `sample.label` as the secondary breadcrumb crumb (which under the old semantics was the friendly text). After the swap, the editable input below shows `display_name` and the breadcrumb above should show the stable `name` (identifier). Otherwise the same string renders twice.
- **Keep the `data-testid="sample-name-input"` selector unchanged** even though it now edits `display_name`. The live spec `e2e/live/sample-rename-preserves-fields.spec.ts:71,116` queries by this testid; renaming it would break the spec for cosmetic reasons. Add an inline comment explaining the testid is historical.

**`mutators/trivial.ts::updateSampleMutator` — load-bearing queue change.**

`packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts` registers `updateSampleMutator` for `kind: "update_sample"` (referenced in `mutatorRegistry.ts:62–63`). Every public surface hardcodes `name` and `label`:

- Line 40: `UpdateSampleInput = { name?: string; notes?: string }` → `{ display_name?: string; notes?: string }`.
- Lines 76–85 (`onSuccess` SSE-wins merge):
  ```ts
  if (response.label !== undefined) patch.label = response.label;   // DELETE this line entirely
  if (response.name !== undefined) patch.name = response.name;       // KEEP — name is now the stable id
  // ADD:
  if (response.display_name !== undefined) patch.display_name = response.display_name;
  ```
- Lines 98–103 (`patchOf` helper): drop `name` from the optimistic patch (immutable now), add `display_name`.

Without these changes the optimistic write puts `{ name: "x" }` into the cache while the server PATCH carries `{ display_name: "x" }` — silent reconciliation drift. **This was missed in the original spec; correcting it is mandatory.**

**`replayCoordinator.ts::synthesizeResponseFromSse` — verified clean.** The default branch `return { ...base, ...payload }` for `update_sample` is payload-shape-agnostic; the rename of payload field flows through. The fix in `trivial.ts` above closes the SSE-wins reconciliation path because `onSuccess` now reads the right field.

**`api.ts::updateExperiment` fetcher type cleanup.** The fetcher at `api.ts:96` types `Partial<Pick<Experiment, "name" | "data_dir" | "analysis_dir" | "manifest_path">>` but the route at `routes_experiments.jl:52` only ever accepted `:name`. Removing `:name` per §6 leaves a fetcher whose every call would 400. Replace with `Partial<Pick<Experiment, "name">>`-then-`Partial<{}>` (i.e. drop the unused fields and `name`). The route still exists as defensive surface.

**Frozen-shape contract test (`test/queue/cache-shape.test.ts:53`):**

```diff
-  "id", "experiment_id", "label", "name", "notes", "tags",
+  "id", "experiment_id", "name", "display_name", "notes", "tags",
```

This is the pin that catches drift between SSE payload, cache row, and `applyRemoteToCache` reads. Must change in lockstep with the backend payload rename and the Sample type.

**Test fixtures (comprehensive sweep).** Every fixture builder that constructs a Sample needs `label` removed and `display_name` added. The `"falls back to sample.label when sample.name is missing"` test in `ExposureListRow.test.tsx:56,60` becomes `"falls back to sample.name when sample.display_name is missing"`. Inventory of files (verified by grep):

Vitest unit tests:
- `test/queue/cache-shape.test.ts:53,329,335,355` — `SAMPLE_KEYS` Set + per-test fixtures + mock response + input field
- `test/queue/sseEventPayload.contract.test.ts:336,342` — fixture + payload field
- `test/queue/rollbackSymmetry.test.ts:56,158–163,213–227` — shared fixture + updateSample case
- `test/queue/authHeaders.test.ts:157–162` — updateSample request mock
- `test/queue/mutatorOnSseWins.test.ts:393–428` — entire SSE-wins block
- `test/queue/treats404AsSuccess.test.tsx:163`
- `test/queries-samples.test.tsx:38,40,55,74`
- `test/SamplePickerRow.test.tsx:7`
- `test/ComparisonPicker.test.tsx:45,50`
- `test/ComparisonPickerPanel.test.tsx:9,16`
- `test/ComparisonPickerBody.test.tsx:12,17,115,123`
- `test/ComparePageEdit.test.tsx:638,644`
- `test/ComparePage.test.tsx:153,296` (legacy comments referencing `sample.label · exposure.filename`)
- `test/NavModal.test.tsx:15,16,20`
- `test/TitleButton.test.tsx:15`
- `test/smoke.test.tsx:56`
- `test/ExposureListRow.test.tsx:36,56,60`

Playwright mocked specs (each hand-rolls Sample objects in `page.route` mocks):
- `e2e/inspect.spec.ts:8`
- `e2e/smoke.spec.ts:8,9`
- `e2e/compare.spec.ts:28,29`
- `e2e/bundle-b-paper-cuts.spec.ts:25,26`
- `e2e/figure-export.spec.ts:31`
- `e2e/multiplayer-replay-rerun.spec.ts:105`

Playwright live specs:
- `e2e/live/sample-rename-preserves-fields.spec.ts` — semantically inverted post-refactor:
  - Line 28 (Sample type with `label`) → drop `label`, add `display_name`.
  - Line 92 (PATCH body `{ name, notes }`) → `{ display_name, notes }`.
  - Lines 116–124 (operates `sample-name-input`, asserts server `name` reflects rename) → asserts on `display_name` instead.
  - The `afterAll` cleanup PATCH at line 92 sends `{ name: fx.originalName, notes }` — must become `{ display_name }` or it silently 400s.

**Mutation queue persisted-state version bump (`lib/queue/persistence.ts:6`).** Bump `SCHEMA_VERSION` from 1 to 2. Pre-deploy clients may have queued `update_sample` ops in `sessionStorage` with `payload: { name: "X" }`; without the bump, they'd rehydrate against the new mutator, the new `patchOf` finds nothing to apply (silent no-op), and the HTTP retry returns 400 from the new route. Bumping to 2 drops the entry as schema-mismatched (counted as `dropped`, not `failed`) and the user can redo the edit.

**No SSE event-kind change.** `update_sample` keeps its kind name; only the payload shape rotates.

## 8. Build sequence

One PR, commit-by-commit, each commit independently green. Order matters because the frontend frozen-shape test ties Sample type, payload, and cache row together.

1. **Backend: schema migration** — `migrate_samples_naming!` in `db.jl`: sentinel-guarded, wrapped in `SQLite.transaction`, runs the four-statement SQL + duplicate-suffix pass (ORDER BY id ASC) + UNIQUE INDEX + final `DELETE FROM idempotent_responses`. Update `create_schema!` so fresh DBs define the canonical column order. Synthetic legacy-shape inline fixtures in `test_db.jl` (matching the existing `migrate_pk_to_autoincrement!` test pattern — no vendored prod DB needed for unit tests). Cases: fresh DB → new shape; legacy DB with mixed label/name → migrated correctly; legacy DB with duplicate names → suffix-renamed with warning, ORDER BY id ASC determinism; idempotent on second run; placement after `migrate_pk_to_autoincrement!` verified.
2. **Backend: validation module** — `validate.jl` + `ManifestViolation` + `validate_manifest`. Pure-function tests in new `test_validate.jl`. No call-site changes yet.
3. **Backend: parser cutover** — `config.jl` (struct rename `col_label → col_display_name` + reject old `[manifest].label` key), `manifest.jl` (`ManifestSample.label → display_name`, column-meaning swap in `parse_manifest`), `simple.toml` (key rename). Updates to `test_config.jl` + `test_manifest.jl`.
4. **Backend: migrate-toml CLI** — new `cli_migrate_toml` with section-aware regex-anchored substitution, dispatch in `main()`, new `test_migrate_toml.jl` (happy path, idempotency, malformed, both old + new keys, inline-comment preservation, the "axis label units" comment in `[beamline]` not misfiring).
5. **Backend: route + pipeline cutover** — `create_sample!` kwarg swap and all 8 call sites in `cli.jl` (incl. `:311, :383` filter + `:357, :374` help text), `routes_samples.jl` allowlist (`display_name` + trim), `routes_experiments.jl` allowlist + 400/500 split, `routes_export.jl` JSON+CSV column rename, `comparisons.jl::picker_samples` SELECT, `routes_experiments.jl:95` log string, `validate_manifest` integrated into both write paths inside transactions. **`_cli_init_inner!` extraction** with `SQLite.transaction` wrap (auto-analyze stays outside). Updates to `test_routes_samples.jl`, `test_routes_experiments.jl`, `test_routes_export.jl`, `test_pipeline.jl`, `test_idempotency_replay_invariant.jl` (every testset's `create_sample!(label="D1")` fixture call + the `:name` PATCH body at line 207).
6. **Frontend: type + helper + atomic rename** — `Sample` type rename, new `lib/sample/displayName.ts` (using `||`, parametrised on `Pick<Sample, …>`), `lib/queue/mutators/trivial.ts::updateSampleMutator` rewrite (UpdateSampleInput, onSuccess, patchOf), `SampleMetadataCard.tsx` full rewrite (input now edits display_name, breadcrumb shows stable name), `lib/queue/persistence.ts::SCHEMA_VERSION` bump 1→2, search both fields in NavModal + ComparisonPickerBody, all 17 component sites + 18 Vitest fixtures + 6 mocked Playwright specs + 1 live Playwright spec switched in one commit. `tsc --noEmit + vite build + npm test + npm run e2e` all green.
7. **Docs** — update `docs/experiment-config.md` with the new `[manifest]` keys and the deprecation message; update `CLAUDE.md` "samples" gotchas if any are now stale; update root `AGENTS.md` if the schema is duplicated there. Verify `docs/superpowers/specs/2026-04-22-himalaya-web-app-design.md` hasn't gone stale on the `samples` shape. CHANGELOG entry calling out: (a) wire-format change to `routes_export.jl` CSV column name, (b) one-shot purge of `idempotent_responses` on first migrated boot.

## 9. Testing strategy

Per-commit unit/integration tests are listed inline in §8. Beyond that:

- **Six-layer audit for `update_sample`** per `docs/contract-testing.md`: every layer touched and pinned.
  - L1 route emit (`routes_samples.jl:42`): payload field is `display_name`.
  - L2 SSE payload (`sseEventPayload.contract.test.ts:336,342`): fixture and assertion updated.
  - L3 `applyRemoteToCache` merge: generic spread; verified payload-shape-agnostic — no code change.
  - L4 cache row (`cache-shape.test.ts:53,329,335,355`): `SAMPLE_KEYS` + per-test fixtures updated.
  - L5 `onMutate` (`mutators/trivial.ts:98–103` `patchOf`): rewritten — was missed in the original spec.
  - L6 `onSuccess` (`mutators/trivial.ts:76–85` SSE-wins merge): rewritten — was missed in the original spec.
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
- [ ] Frontend renders `sampleDisplayName(s)` everywhere; no remaining `s.label` or `sample.label` references in `src/`, `test/`, or `e2e/`.
- [ ] `mutators/trivial.ts::updateSampleMutator` reads/writes `display_name`; no `response.label` reference remains.
- [ ] `lib/queue/persistence.ts::SCHEMA_VERSION === 2`; pre-deploy queued ops are dropped on rehydrate.
- [ ] `routes_export.jl` CSV header reads `sample_display_name`; downstream-consumer note in CHANGELOG.
- [ ] `comparisons.jl::picker_samples` SELECT projects `display_name` (no `label`).
- [ ] All contract tests green (cache-shape, SSE payload, route response shapes, idempotency replay incl. fixture bodies).
- [ ] Full backend Julia suite green (`Pkg.test("HimalayaUI")`); full frontend Vitest + Playwright mocked green.
- [ ] Live smoke against a prod-backup copy completes without UI regression; `e2e/live/sample-rename-preserves-fields.spec.ts` updated and green against a migrated dev backend.
