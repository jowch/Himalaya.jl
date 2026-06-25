# Ingestion Redesign — Phase E1: Frontend Shell, IA, Routing, API Client, Queue Wiring, Reusable Dropdown

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. TDD is mandatory: write the failing test, run it red, write the minimal implementation, run it green, commit.

**Goal:** Build the frontend skeleton for the ingestion redesign — the new information architecture (top nav `Experiments | Series`; per-experiment tabs `Corpus | Configuration`; grouping surfaced via a banner, not a tab), the `/experiments` routing tree under a dedicated `ExperimentShell` that renders its own chrome, the page skeletons (`ExperimentsHomePage`, `NewExperimentPage` with a directory-path picker, `ExperimentCorpusPage`, `ConfigurationPage` shell), the `api.ts` additions for scan/ingest/experiment endpoints, the `ingest_*` SSE receive path (broadcast-only progress) plus the SSE receive path for the structural event kinds, and **two new reusable `src/print/ui/` primitives** — `Dropdown` (a labelled select-trigger + Menu popover, consumed by the directory picker here and exported for E2) and `StatBar` (the hairline-divided stat ledger).

**Architecture:** Everything lands under `packages/HimalayaUI/frontend/src/`. The new `/experiments/*` routes sit **outside** `CorpusShell` (so the experiment header chrome never stacks on `CorpusTopbar`); `ExperimentShell` renders its own top chrome (wordmark + `Experiments | Series` nav + experiment header + `Corpus | Configuration` tab bar + `<Outlet>`). Server state stays in TanStack Query (`queries.ts` hooks + `api.ts` fetchers); client state stays in Zustand (`state.ts` named actions); the mutation queue's `applyRemoteToCache` stays **pure** (cache-only) while ephemeral `ingestInFlight` store writes live in a **separate `App.tsx` SSE listener** per spec §9.3/§9.6. Appearance lives only in `src/print/ui/` primitives; consumers pass placement-only `className` (mechanically enforced by `scripts/check-design.mjs`).

**Tech Stack:** React 18 + Vite + TypeScript strict (`exactOptionalPropertyTypes`) + TailwindCSS 4, TanStack Query, Zustand, react-router-dom v6, Vitest + React Testing Library + jsdom. Frontend package at `packages/HimalayaUI/frontend/`.

**Spec:** `docs/superpowers/specs/2026-06-15-ingestion-redesign-design.md` §7 (IA), §8.6/§8.7 (component reuse + new components), §9.2 (REST endpoints), §9.3 (event kinds + `broadcast_progress!` frame shape), §9.6 (frontend wiring), §13 (open questions). Visual + interaction reference: `docs/redesign-mockups/ingestion-prototype.html` (canonical). Read both before starting.

**Backend contract this plan calls:** `docs/superpowers/plans/2026-06-18-ingestion-phase-a-schema.md` defines the data shapes the API returns — `loads` table, `exposures.experiment_id`, typed geometry columns + per-field `*_source`, `last_scanned_at`/`scan_signature`/`ingest_status`, and the `samples` `display_name → name` collapse. This plan's `api.ts` types mirror those.

---

## Scope boundary (E1 vs E2)

**E1 (this plan) owns:** IA/routing/chrome, the four page skeletons, the directory picker field, `Dropdown` + `StatBar` primitives, `api.ts` fetchers + the **canonical nested `Load`/`LoadSample`/`LoadExposure`/`GroupingFlag` types + `queryKeys.loads` + `useLoads`** (E2 imports these, must NOT redefine), the `display_name → name` collapse + consumer sweep (Task 1b), the `ingestInFlight` store slice, the `ingest_*` SSE handling (App-listener short-circuit: store write + shared `invalidateIngestFrameCache` helper call; the canonical arms in `applyRemoteToCache` are the tested path for non-App callers), and the SSE RECEIVE path for the structural event kinds (`sample_renamed`, `exposure_moved`, `sample_created`, `sample_split`, `grouping_flag_dismissed` — **there is no `sample_merged`**).

**E2 (a later plan) owns:** the grouping-review surface internals — `GroupingReviewPage`, `LoadFold`/`SampleFold`/`ExposureLeaf`, the bulk-select bar wiring, `GeometryLedger`, `AcquisitionTimeline`, `SourcesCard`, the `useUndoStack` extraction, and the structural MUTATIONS (merge/split/rename/move mutators + `useQueueMutation` wiring). Where this plan slots an E2 component, it renders a clearly-labelled placeholder and cross-references the E2 component by name. **Do not build E2's components here.**

**E1 also owns (backend, Task 1c):** the additive migration (`description`/`image_pattern`/`metadata_pattern`/`integration_pattern` columns on experiments) + widening `PATCH /api/experiments/:id` to accept `name`/`description`/patterns (dropping them from `_READONLY_FIELDS`; rescan-on-pattern-change via `scan_signature = NULL`). This is the minimum backend that makes the ExperimentShell name/description edit-in-place and the E2 SourcesCard pattern rows work against a real backend. Geometry ×6 were already built in Phase C.

**Out of scope entirely (other phases):** `prp.jl`/`geometry.jl`/`grouping.jl`/`scan_and_group!` = Phase B; scan/rescan routes + SSE + scheduler = Phase C; structural-edit event kinds backend = Phase D; the Phase A schema migrations (loads table, exposure grouping cols, dedup key). This plan consumes those endpoints/frames; it does not implement them. Against a backend that has not shipped them, the new hooks simply return empty/loading — the tests here mock `fetch`/the API layer, so they never depend on a live backend.

---

## File map

| File | Responsibility | This plan |
|---|---|---|
| `packages/HimalayaUI/src/db.jl` (or migration runner) | additive `description`/`*_pattern` columns | MODIFY (Task 1c) |
| `packages/HimalayaUI/src/routes_experiments.jl` | widen PATCH: drop name/description from `_READONLY_FIELDS`; pattern rescan | MODIFY (Task 1c) |
| `packages/HimalayaUI/src/ingest.jl` (or scan entry) | scanner reads pattern cols with TOML fallback | MODIFY (Task 1c) |
| `packages/HimalayaUI/test/test_experiments_patch.jl` | backend PATCH widening tests | CREATE (Task 1c) |
| `src/api.ts` | fetchers + types (incl. nested `Load`/`LoadSample`/`LoadExposure`/`GroupingFlag`) | MODIFY (Tasks 1, 1b, 2) |
| `src/lib/sample/displayName.ts`, `lib/comparison/labels.ts`, `lib/scoping/proposeOrdering.ts`, `lib/queue/mutators/trivial.ts`, `print/components/AddSamplePicker.tsx`, `print/pages/{FocusPage,LoupePage,SeriesBuilderPage,SeriesScopingPage,builderAdapters}.tsx`, `print/shell/NavModal.tsx` | `display_name → name` consumer sweep | MODIFY (Task 1b) |
| `src/state.ts` | Zustand store | MODIFY (Task 3) |
| `src/queries.ts` | TanStack hooks + `queryKeys` | MODIFY (Task 4) |
| `src/print/ui/Dropdown.tsx` | new primitive | CREATE (Task 5) |
| `src/print/ui/StatBar.tsx` | new primitive | CREATE (Task 6) |
| `src/print/ui/index.ts` | primitive barrel | MODIFY (Tasks 5, 6) |
| `src/print/components/PageFrame.tsx` | width keys | MODIFY (Task 7) |
| `src/print/components/DirectoryPickerField.tsx` | path picker | CREATE (Task 8) |
| `src/print/shell/ExperimentTopNav.tsx` | `Experiments \| Series` nav | CREATE (Task 9) |
| `src/print/shell/ExperimentShell.tsx` | `/experiments/:id` layout chrome | CREATE (Task 10) |
| `src/print/pages/ExperimentsHomePage.tsx` | `/experiments` gallery | CREATE (Task 11) |
| `src/print/pages/NewExperimentPage.tsx` | `/experiments/new` | CREATE (Task 12) |
| `src/print/pages/ExperimentCorpusPage.tsx` | Corpus tab body | CREATE (Task 13) |
| `src/print/pages/ConfigurationPage.tsx` | Configuration tab body (shell) | CREATE (Task 14) |
| `src/print/shell/AppRoutes.tsx` | route table | MODIFY (Task 15) |
| `src/print/shell/CorpusTopbar.tsx` | add Experiments stage | MODIFY (Task 16) |
| ~~`src/print/shell/NavModal.tsx` (home phase bar)~~ | — | Task 17 DROPPED (no such element exists; NavModal's `display_name` reads are in Task 1b) |
| `src/lib/queue/applyRemoteToCache.ts` | `ingest_*` + structural receive | MODIFY (Tasks 18, 19) — **no `types.ts` change** (`SseEvent.kind` is already `string`) |
| `src/print/App.tsx` | `ingestInFlight` store write + `ingest_*` short-circuit | MODIFY (Task 20) |
| `src/lib/queue/persistence.ts` | `SCHEMA_VERSION` 4→5 | MODIFY (Task 1b; Task 21 folded) |
| `test/*.test.ts(x)` | Vitest unit tests | CREATE (every task) |

**Test convention:** Vitest tests live in `packages/HimalayaUI/frontend/test/` (flat, e.g. `test/AppRoutes.test.tsx`). Run a single file during TDD from `packages/HimalayaUI/frontend/`:

```bash
npm test -- test/<file>.test.ts
```

`npm test` is one-shot (the repo's settings.json adds Vitest's `--run` flag). Commit after each task once its test passes. Do **not** assert on Tailwind class strings — use `data-testid`/`data-*`/roles (AGENTS.md anti-pattern). The final build gate (`npm run build` = `tsc --noEmit` + `lint:design` + `vite build`) runs in Task 22.

---

## Task 1: `api.ts` — extend `Experiment`, add `Sample.name` collapse, add `Load`

Mirror the Phase-A schema: typed geometry + per-field `*_source` + scan columns on `Experiment`; drop `display_name` from `Sample`; add a `Load` interface. Type-only changes plus a label fallback; no fetchers yet (Task 2).

**Files:**
- Modify: `src/api.ts` (`Experiment` interface ~8-22, `Sample` interface ~31-38)
- Test: `test/apiTypes.ingestion.test.ts` (CREATE)

- [ ] **Step 1: Write the failing test**

```ts
// test/apiTypes.ingestion.test.ts
import { describe, it, expect } from "vitest";
import type { Experiment, Sample, Load, GroupingFlag } from "../src/api";

describe("ingestion api types (Phase E1)", () => {
  it("Experiment carries typed geometry source tags + scan columns + E1 additive fields", () => {
    const e: Experiment = {
      id: 1, name: "SSRL", description: "AgBe SAXS run", path: "/d",
      data_dir: "/d/data", analysis_dir: "/d/an",
      manifest_path: null, created_at: "2026-04-12", q_units: "A^-1",
      beam_center_x: 421.4, beam_center_y: 836.9, pixel_size_um: 172, energy_kev: 9,
      flight_path_m: 1.8095,
      energy_kev_source: "prp", flight_path_m_source: "setup",
      beam_center_x_source: "setup", beam_center_y_source: "setup",
      pixel_size_um_source: "prp", q_units_source: "prp",
      last_scanned_at: "2026-04-12T10:00:00", scan_signature: "sig", ingest_status: "idle",
      image_pattern: "*.tif", metadata_pattern: null, integration_pattern: null,
    };
    expect(e.flight_path_m_source).toBe("setup");
    expect(e.ingest_status).toBe("idle");
    expect(e.description).toBe("AgBe SAXS run");
    expect(e.image_pattern).toBe("*.tif");
  });

  it("Sample has a single non-null `name` label and no `display_name`", () => {
    const s: Sample = {
      id: 1, experiment_id: 1, name: "HA85 (S01P15)", notes: null, tags: [],
    };
    expect(s.name).toBe("HA85 (S01P15)");
    // @ts-expect-error display_name is removed from Sample
    expect(s.display_name).toBeUndefined();
    // @ts-expect-error name is non-null after the collapse — null is not assignable
    const _nn: Sample = { id: 2, experiment_id: 1, name: null, notes: null, tags: [] };
    void _nn;
  });

  it("Load is the NESTED roll-up (Load ▸ Sample ▸ Exposures)", () => {
    const l: Load = {
      load_id: 3, load_index: 1, session_id: 1,
      start_time: "2026-04-12T10:00:00", end_time: "2026-04-12T10:10:00",
      frame_count: 64, note: null,
      samples: [
        {
          sample_id: 9, name: "HA85 (S01P15)", slot_index: 15,
          grouping_source: "computed", name_source: "computed",
          merged_into_id: null, flag: null,
          exposures: [
            { id: 51, filename: "ha85_001.tif", horizontal_position: 12.4, timestamp: "2026-04-12T10:00:01" },
          ],
        },
      ],
    };
    expect(l.load_index).toBe(1);
    expect(l.samples[0]!.exposures[0]!.id).toBe(51);
    expect(l.samples[0]!.flag).toBeNull();
  });

  it("a Load sample can carry a merge or split GroupingFlag", () => {
    const merge: GroupingFlag = { kind: "merge", merge_with_sample_id: 4, merge_with_label: "HA85 (S01P14)" };
    const split: GroupingFlag = { kind: "split", split_at_index: 32, jump_from: 12.4, jump_to: 48.1 };
    expect(merge.kind).toBe("merge");
    expect(split.kind).toBe("split");
  });
});
```

> The test imports `GroupingFlag` alongside `Experiment, Sample, Load` (add it to the import line). `Load` is the **nested** roll-up E2 consumes — E1 OWNS this shape; E2 must NOT redefine it. The flat-row shape an earlier draft used is wrong (it can't drive `LoadFold`/`SampleFold`/`ExposureLeaf`); the nested shape below is canonical (pinned in spec §8.8).

- [ ] **Step 2: Run it red** — `npm test -- test/apiTypes.ingestion.test.ts` → FAIL (`Load` not exported; `*_source`/`ingest_status`/`description`/`image_pattern` missing on `Experiment`; `name` not required / `display_name` still present on `Sample`; the `@ts-expect-error` may even fail to compile the suite).

- [ ] **Step 3: Implement**

Replace the `Experiment` interface (currently `src/api.ts:8-22`):

```ts
export type GeometrySource = "prp" | "setup" | "user" | "default";
export type IngestStatus = "idle" | "scanning" | "analyzing" | "complete" | "failed";

export interface Experiment {
  id: number;
  name: string | null;
  description: string | null;     // Phase E1 additive column (Task 1c migration)
  path: string;
  data_dir: string;
  analysis_dir: string;
  manifest_path: string | null;
  created_at: string;
  q_units: string | null;
  beam_center_x: number | null;
  beam_center_y: number | null;
  pixel_size_um: number | null;
  energy_kev: number | null;
  flight_path_m: number | null;
  // Phase A typed-geometry per-field provenance + scan bookkeeping.
  energy_kev_source: GeometrySource;
  flight_path_m_source: GeometrySource;
  beam_center_x_source: GeometrySource;
  beam_center_y_source: GeometrySource;
  pixel_size_um_source: GeometrySource;
  q_units_source: GeometrySource;
  last_scanned_at: string | null;
  scan_signature: string | null;
  ingest_status: IngestStatus;
  // Phase E1 additive columns: editable file-pattern globs (Task 1c migration).
  // NULL = use legacy experiment.toml fallback.
  image_pattern: string | null;
  metadata_pattern: string | null;
  integration_pattern: string | null;
}
```

Replace the `Sample` interface (currently `src/api.ts:31-38`) — drop `display_name`, make `name` non-null (the display_name→name collapse always populates `name` server-side; spec §8.8):

```ts
export interface Sample {
  id: number;
  experiment_id: number;
  name: string;            // non-null after the collapse (was `string | null` + `display_name`)
  notes: string | null;
  tags: SampleTag[];
}
```

Add the **nested** `Load` roll-up cluster immediately after `Sample` (before `CorpusSample`). **E1 OWNS these types; E2 imports them and must NOT redefine them** (pinned in spec §8.8):

```ts
/** A merge/split discrepancy flag the auto-grouper raised on a slot (spec §8.8).
 *  `null` when the slot is clean. The structural-edit dismissal arm
 *  (`grouping_flag_dismissed`, Phase D) clears it. */
export type GroupingFlag =
  | { kind: "merge"; merge_with_sample_id: number; merge_with_label: string }
  | { kind: "split"; split_at_index: number; jump_from: number; jump_to: number }
  | null;

/** One exposure leaf under a load's sample slot. */
export interface LoadExposure {
  id: number;
  filename: string;
  horizontal_position: number | null;
  timestamp: string | null;
}

/** One sample slot inside a load (a (load, slot) coordinate). `name_source`/
 *  `grouping_source` are provenance tags ("user" | "computed" | …);
 *  `merged_into_id` is non-null when this slot was merged into a sibling. */
export interface LoadSample {
  sample_id: number;
  name: string;
  slot_index: number;
  grouping_source: string;
  name_source: string;
  merged_into_id: number | null;
  flag: GroupingFlag;
  exposures: LoadExposure[];
}

/** One rack-load roll-up (Phase A `loads` table), returned NESTED by
 *  `GET /api/experiments/:id/loads` (see Task 2): Load ▸ Sample ▸ Exposures.
 *  Drives E2's LoadFold/SampleFold/ExposureLeaf + the grouping-review count
 *  (samples whose `flag` is non-null). */
export interface Load {
  load_id: number;
  load_index: number;
  session_id: number | null;
  start_time: string | null;
  end_time: string | null;
  frame_count: number;
  note: string | null;
  samples: LoadSample[];
}
```

> **Label sweep note (spec §8.8 / §9.6):** making `name` non-null and dropping `Sample.display_name` is a type change that will NOT compile until every consumer is converted. The full atomic sweep is **Task 1b** (the next task, sequenced immediately after this one) — it converts all consumers AND bumps the queue-op `SCHEMA_VERSION` 4→5 so the repo's `tsc --noEmit` is green again at Task 1b's commit. This task (Task 1) intentionally leaves the tree red between its own commit and Task 1b's; do them back-to-back. (The verified consumer list lives in Task 1b.)

- [ ] **Step 4: Flip `updateSample`'s patch key**

In `src/api.ts:134`, change the patch type from `{ display_name?: string; notes?: string }` to `{ name?: string; notes?: string }`:

```ts
export const updateSample = (id: number, patch: { name?: string; notes?: string }, opts?: AuthOpts) =>
  request<Sample>("PATCH", `/api/samples/${id}`, patch, opts);
```

- [ ] **Step 5: Run it green** — `npm test -- test/apiTypes.ingestion.test.ts` → PASS. (The broader suite + `tsc` are RED until Task 1b lands — that is expected; this commit is intentionally a non-green intermediate paired with Task 1b.)

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/api.ts packages/HimalayaUI/frontend/test/apiTypes.ingestion.test.ts
git commit -m "feat(api): typed geometry source tags + scan cols on Experiment; Sample.name collapse; nested Load roll-up"
```

---

## Task 1b: `display_name → name` consumer sweep + `SCHEMA_VERSION` 4→5

Task 1 made `Sample.name` non-null and removed `display_name`, which breaks `tsc` across every consumer that read `display_name` or treated `name` as nullable. This task converts them all in ONE commit so `tsc --noEmit` is green again, AND bumps the queue-op `SCHEMA_VERSION` 4→5 so a pre-deploy `update_sample` op carrying the old `display_name` key DROPS on rehydrate (the existing "N edits couldn't be restored" toast).

**Files (verified live via `git grep -n display_name src/`, 2026-06-18 — all non-test occurrences):**
- `src/api.ts` — `Sample.display_name` (line 35, removed in Task 1); `updateSample` patch key (134, flipped in Task 1).
- `src/lib/sample/displayName.ts` — `sampleDisplayName`'s `Pick<Sample, "id" | "name" | "display_name">` → `Pick<Sample, "id" | "name">`; body `s.display_name || s.name || \`Sample #${s.id}\`` → `s.name || \`Sample #${s.id}\``.
- `src/lib/comparison/labels.ts` — `sample.display_name || sample.name` (line 38) → `sample.name`; drop the `display_name` mentions in the doc comment (7/9).
- `src/lib/scoping/proposeOrdering.ts` — `s.sample.display_name ?? s.sample.name ?? ""` (55) → `s.sample.name ?? ""` (name is now non-null, so the `?? ""` is dead but harmless; can simplify to `s.sample.name`).
- `src/lib/queue/mutators/trivial.ts` — `UpdateSampleInput` (41) `{ display_name?; notes? }` → `{ name?; notes? }`; `patchOf` (99-104) `display_name` → `name`; `onSuccess` (84-86) drop the `if (response.display_name !== undefined) patch.display_name = …` line (the `response.name` line at 84 already handles it).
- `src/print/components/AddSamplePicker.tsx` (17) — `s.display_name ?? s.name ?? \`Sample ${s.id}\`` → `s.name`.
- `src/print/pages/FocusPage.tsx` (344, 348) — `corpusSample?.display_name ?? corpusSample?.name ?? …` → `corpusSample?.name ?? …` (corpusSample is optional, so keep the `?.` + fallback).
- `src/print/pages/LoupePage.tsx` (397, 456) — `sample?.display_name ?? sample?.name ?? …` → `sample?.name ?? …`.
- `src/print/pages/SeriesBuilderPage.tsx` (221, 1109) — `… .display_name ?? … .name ?? \`Sample …\`` → `… .name ?? \`Sample …\``.
- `src/print/pages/SeriesScopingPage.tsx` (332) — `s.sample.display_name ?? s.sample.name ?? ""` → `s.sample.name`.
- `src/print/pages/builderAdapters.ts` (84) — doc comment only ("id → display_name ?? name ?? …" → "id → name").
- `src/print/shell/NavModal.tsx` (117, 128, 147, 148, 257, 371) — every `display_name` read collapses to `name`. Note 147/148 build a primary/secondary pair from `display_name` vs `name`; with one field, the secondary (`s.name && s.display_name && s.name !== s.display_name ? s.name : ""`) becomes `""` (drop the secondary, or keep an empty string). Simplify: `primary: s.name || \`Sample ${s.id}\``, `secondary: ""`.
- `src/queries.ts` (347) — doc comment only ("`${sample.display_name || sample.name} · …`" → "`${sample.name} · …`").
- `src/lib/queue/persistence.ts` (14) — `SCHEMA_VERSION = 4` → `5`. Update the comment with a "Bumped 4 → 5: samples.display_name collapsed to name; pre-deploy update_sample ops carrying display_name DROP on rehydrate" line.

> **Re-verify before editing:** run `git grep -n display_name packages/HimalayaUI/frontend/src` once more at implementation time and reconcile against this list (a sibling branch may have touched a file). The list above was captured from live source on 2026-06-18.

- [ ] **Step 1: Write the failing test** — pin the two load-bearing behaviours (helper collapse + schema bump):

```ts
// test/displayNameCollapse.test.ts
import { describe, it, expect } from "vitest";
import { sampleDisplayName } from "../src/lib/sample/displayName";
import { SCHEMA_VERSION } from "../src/lib/queue/persistence";

describe("display_name → name collapse (Phase E1, Task 1b)", () => {
  it("sampleDisplayName reads name only", () => {
    expect(sampleDisplayName({ id: 7, name: "HA85 (S01P15)" })).toBe("HA85 (S01P15)");
    expect(sampleDisplayName({ id: 7, name: "" })).toBe("Sample #7");
  });
  it("queue-op SCHEMA_VERSION is bumped to 5 for the label collapse", () => {
    expect(SCHEMA_VERSION).toBe(5);
  });
});
```

> The `sampleDisplayName` call drops the `display_name` field from its argument — that compiles only after the `Pick<>` is narrowed, so this test red-then-greens the helper signature too. (This subsumes the standalone `persistenceVersion.ingestion.test.ts` Task 21 used to add — Task 21 is folded here and removed.)

- [ ] **Step 2: Run it red** — `npm test -- test/displayNameCollapse.test.ts` → FAIL (the arg without `display_name` won't compile / `SCHEMA_VERSION === 4`).

- [ ] **Step 3: Convert every consumer** in the file list above, then bump `SCHEMA_VERSION`.

- [ ] **Step 4: Run it green** + `tsc` — `npm test -- test/displayNameCollapse.test.ts` → PASS, then `npx tsc --noEmit -p tsconfig.json` (or `tsconfig.build.json` with `--noEmit`) → no `display_name` / nullable-`name` errors. Fix any straggler the grep missed.

- [ ] **Step 5: Commit**

```bash
git add -u packages/HimalayaUI/frontend/src packages/HimalayaUI/frontend/test/displayNameCollapse.test.ts
git commit -m "refactor: collapse samples.display_name into non-null name (all consumers) + bump queue SCHEMA_VERSION 4->5"
```

---

## Task 1c: Backend — additive migration + widen `PATCH /api/experiments/:id`

**Before writing any frontend code for name/description/pattern editing (Tasks 2 and 10), the backend must be widened here.** This task is pure Julia; it touches no frontend files.

**Background:** Phase C built `PATCH /api/experiments/:id` for geometry-only overrides and put `name` and `description` in `_READONLY_FIELDS` because no `description` column existed and the Phase-E1 migration hadn't landed. This task:

1. Adds the missing columns via an additive migration (geometry-cols pattern — safe, backward-compat).
2. Drops `name` and `description` from `_READONLY_FIELDS` so they are plain-writable (no `*_source` stamp — these are identity/documentation fields, not geometry).
3. Makes `image_pattern`, `metadata_pattern`, `integration_pattern` patchable, where each write invalidates `scan_signature` (setting it to `NULL`) so the next scan re-discovers with the new glob.
4. Keeps `data_dir`, `analysis_dir`, `manifest_path`, `path`, `id`, `created_at` read-only.

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl` (or the migration runner, whichever applies the additive DDL) — add `description TEXT`, `image_pattern TEXT`, `metadata_pattern TEXT`, `integration_pattern TEXT` columns to `experiments`
- Modify: `packages/HimalayaUI/src/routes_experiments.jl` — drop `"name"`, `"description"` from `_READONLY_FIELDS`; add `"image_pattern"`, `"metadata_pattern"`, `"integration_pattern"` as patchable with rescan-on-change (set `scan_signature = NULL`); keep `"data_dir"`, `"analysis_dir"`, `"manifest_path"`, `"path"`, `"id"`, `"created_at"` read-only
- Modify: `packages/HimalayaUI/src/ingest.jl` (or the scan entry point) — scanner reads `image_pattern`/`metadata_pattern`/`integration_pattern` from the DB row first, falling back to the legacy `config` TOML blob when the column is `NULL` (experiment.toml is being retired)
- Test: `packages/HimalayaUI/test/test_experiments_patch.jl` (CREATE or extend existing) — register in `runtests.jl` `ALL_ORDER` + one `GROUPS` bucket

> **Test harness (as-built, must match):** route tests dispatch in-process: `with_inproc_routes(db) do call ... end`; fresh DB from `open_prepared_clone(tmp)`. Body bytes: `Vector{UInt8}(JSON3.write(Dict(...)))`. Headers: `["Content-Type"=>"application/json", "X-Username"=>"alice"]`. `call` never throws on 4xx/5xx; no `status_exception=false` needed.

> **Geometry ×6 already work** — Phase C built those (field + `*_source = 'user'` write). This task only widens the route further; do NOT re-spec or re-implement the geometry arm.

- [ ] **Step 1: Write the failing test**

```julia
# packages/HimalayaUI/test/test_experiments_patch.jl
using Test, JSON3, HimalayaUI
using HimalayaUI: with_inproc_routes, open_prepared_clone

@testset "PATCH /api/experiments/:id — widened fields" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        with_inproc_routes(db) do call
            # seed a minimal experiment row
            body_create = Vector{UInt8}(JSON3.write(Dict("path" => tmp, "name" => "test")))
            r = call("POST", "/api/experiments"; headers=["Content-Type"=>"application/json","X-Username"=>"alice"], body=body_create)
            exp_id = JSON3.read(r.body)["id"]

            # name is patchable (not in _READONLY_FIELDS)
            r2 = call("PATCH", "/api/experiments/$exp_id";
                headers=["Content-Type"=>"application/json","X-Username"=>"alice"],
                body=Vector{UInt8}(JSON3.write(Dict("name" => "renamed"))))
            @test r2.status == 200
            r2b = call("GET", "/api/experiments/$exp_id")
            @test JSON3.read(r2b.body)["name"] == "renamed"

            # description is patchable
            r3 = call("PATCH", "/api/experiments/$exp_id";
                headers=["Content-Type"=>"application/json","X-Username"=>"alice"],
                body=Vector{UInt8}(JSON3.write(Dict("description" => "AgBe SAXS run"))))
            @test r3.status == 200
            r3b = call("GET", "/api/experiments/$exp_id")
            @test JSON3.read(r3b.body)["description"] == "AgBe SAXS run"

            # image_pattern is patchable AND resets scan_signature
            r4 = call("PATCH", "/api/experiments/$exp_id";
                headers=["Content-Type"=>"application/json","X-Username"=>"alice"],
                body=Vector{UInt8}(JSON3.write(Dict("image_pattern" => "*.tif"))))
            @test r4.status == 200
            r4b = call("GET", "/api/experiments/$exp_id")
            exp4 = JSON3.read(r4b.body)
            @test exp4["image_pattern"] == "*.tif"
            @test ismissing(exp4["scan_signature"]) || exp4["scan_signature"] === nothing

            # data_dir stays read-only (400)
            r5 = call("PATCH", "/api/experiments/$exp_id";
                headers=["Content-Type"=>"application/json","X-Username"=>"alice"],
                body=Vector{UInt8}(JSON3.write(Dict("data_dir" => "/other"))))
            @test r5.status == 400
        end
    end
end
```

- [ ] **Step 2: Run it red** — `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' 2>&1 | grep -A5 "widened"` → FAIL (name/description 400, pattern fields unknown column or 400, scan_signature not reset).

- [ ] **Step 3: Implement**

**Migration (additive, geometry-cols precedent):** In the migration runner / `db.jl`, add:

```julia
# Phase E1: editable config columns on experiments
# Additive — safe to run against old DBs (ADD COLUMN does not rewrite existing rows).
for col in [
    "description TEXT",
    "image_pattern TEXT",
    "metadata_pattern TEXT",
    "integration_pattern TEXT",
]
    try
        DBInterface.execute(db, "ALTER TABLE experiments ADD COLUMN $col")
    catch
        # column already exists (idempotent)
    end
end
```

**Route widening (`routes_experiments.jl`):** Change `_READONLY_FIELDS` from:

```julia
const _READONLY_FIELDS = ["data_dir", "analysis_dir", "manifest_path", "path",
                           "id", "created_at", "name", "description"]
```

to (drop `"name"` and `"description"`):

```julia
const _READONLY_FIELDS = ["data_dir", "analysis_dir", "manifest_path", "path",
                           "id", "created_at"]
```

**Extend the SET-clause builder.** The live handler (`routes_experiments.jl:222-248`) builds `set_clauses::Vector{String}` + `params::Vector{Any}` by looping **only over `_GEOMETRY_PATCH_FIELDS`** (geometry value `+` `<field>_source = 'user'`), then 400s if `isempty(set_clauses)`. Dropping `name`/`description` from `_READONLY_FIELDS` is **necessary but not sufficient** — without writing them into `set_clauses` they fall through to that 400. Add two more module consts and two append loops **after** the geometry loop, writing into the **same** `set_clauses`/`params` arrays:

```julia
# (module top, beside _GEOMETRY_PATCH_FIELDS / _READONLY_FIELDS)
const _PLAIN_PATCH_FIELDS = ["name", "description"]                                # plain write, NO *_source, NO rescan
const _PATTERN_FIELDS     = ["image_pattern", "metadata_pattern", "integration_pattern"]  # plain write + rescan trigger
```

```julia
# ... after the existing `for field in _GEOMETRY_PATCH_FIELDS ... end` loop,
#     reusing the SAME set_clauses / params that the geometry loop appended to:

# Plain identity fields: no *_source stamp, no rescan.
for field in _PLAIN_PATCH_FIELDS
    val = get(body, Symbol(field), get(body, field, nothing))
    val === nothing && continue
    push!(set_clauses, "$field = ?"); push!(params, val)
end

# Pattern fields: plain write + invalidate scan_signature so the next scan re-discovers.
pattern_touched = false
for field in _PATTERN_FIELDS
    val = get(body, Symbol(field), get(body, field, nothing))
    val === nothing && continue
    push!(set_clauses, "$field = ?"); push!(params, val)
    pattern_touched = true
end
pattern_touched && push!(set_clauses, "scan_signature = NULL")
```

Now a `name`/`description`/pattern-only body produces non-empty `set_clauses`, so the existing `isempty(set_clauses) && return 400` check passes and the existing `UPDATE experiments SET $(join(set_clauses, ", ")) WHERE id = ?` writes the new fields under `_DB_WRITE_LOCK` + transaction — no other handler changes needed. Also widen the no-patchable-fields 400 message to mention the new accepted fields (currently it names only the geometry fields).

> `name` and `description` are plain writes with no `*_source` stamp and NO rescan trigger; the 3 pattern fields are plain writes that additionally NULL `scan_signature` (one trigger regardless of how many patterns changed). Geometry fields keep their `<field>_source = 'user'` stamp (unchanged).

**Scanner fallback (`ingest.jl` / `scan.jl`):** When reading glob patterns, prefer the DB columns over the legacy TOML blob:

```julia
image_pat = coalesce(exp.image_pattern, get(config, "image_pattern", "*.tif"))
metadata_pat = coalesce(exp.metadata_pattern, get(config, "metadata_pattern", "*.prp"))
integration_pat = coalesce(exp.integration_pattern, get(config, "integration_pattern", "*.dat"))
```

- [ ] **Step 4: Run it green** — suite green on the new test.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/ packages/HimalayaUI/test/test_experiments_patch.jl packages/HimalayaUI/test/runtests.jl
git commit -m "feat(backend): additive migration + widen PATCH /api/experiments/:id (name/description/patterns; rescan-on-pattern-change)"
```

---

## Task 2: `api.ts` — scan/ingest/experiment/loads fetchers + `ExperimentPatch`

Add the create-from-directory, rescan, loads roll-up, experiment PATCH, and directory-picker fetchers (spec §9.2). Widen `updateExperiment` from `Record<string, never>` to the canonical typed `ExperimentPatch`. **E2 imports `ExperimentPatch` from here and must NOT redefine it.**

**Files:**
- Modify: `src/api.ts` (after the existing Experiments block ~118-127)
- Test: `test/apiFetchers.ingestion.test.ts` (CREATE)

- [ ] **Step 1: Write the failing test**

```ts
// test/apiFetchers.ingestion.test.ts
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import * as api from "../src/api";

function mockFetchOnce(body: unknown, ok = true, status = 200): void {
  vi.spyOn(globalThis, "fetch").mockResolvedValueOnce({
    ok, status,
    text: async () => JSON.stringify(body),
    json: async () => body,
  } as Response);
}

describe("ingestion fetchers (Phase E1)", () => {
  beforeEach(() => vi.restoreAllMocks());
  afterEach(() => vi.restoreAllMocks());

  it("createExperiment POSTs {path,name?,patterns?} to /api/experiments", async () => {
    const spy = vi.spyOn(globalThis, "fetch").mockResolvedValueOnce({
      ok: true, status: 200,
      text: async () => JSON.stringify({ id: 7 }),
      json: async () => ({ id: 7 }),
    } as Response);
    const out = await api.createExperiment({ path: "/d", name: "X" }, { username: "u" });
    expect(out.id).toBe(7);
    const [url, init] = spy.mock.calls[0]!;
    expect(url).toBe("/api/experiments");
    expect(init!.method).toBe("POST");
    expect(JSON.parse(init!.body as string)).toEqual({ path: "/d", name: "X" });
    expect((init!.headers as Record<string, string>)["X-Username"]).toBe("u");
  });

  it("triggerScan POSTs to /api/experiments/:id/scan", async () => {
    const spy = mockFetchSpy({ id: 7, ingest_status: "scanning" });
    await api.triggerScan(7);
    expect(spy.mock.calls[0]![0]).toBe("/api/experiments/7/scan");
    expect((spy.mock.calls[0]![1] as RequestInit).method).toBe("POST");
  });

  it("listLoads GETs /api/experiments/:id/loads (nested roll-up)", async () => {
    const spy = mockFetchSpy([{
      load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null,
      frame_count: 0, note: null, samples: [],
    }]);
    const loads = await api.listLoads(7);
    expect(loads).toHaveLength(1);
    expect(loads[0]!.samples).toEqual([]);
    expect(spy.mock.calls[0]![0]).toBe("/api/experiments/7/loads");
  });

  it("updateExperiment PATCHes a geometry partial", async () => {
    const spy = mockFetchSpy({ id: 7, flight_path_m: 1.81, flight_path_m_source: "user" });
    const patch: api.ExperimentPatch = { flight_path_m: 1.81 };
    await api.updateExperiment(7, patch, { username: "u" });
    expect(spy.mock.calls[0]![0]).toBe("/api/experiments/7");
    expect((spy.mock.calls[0]![1] as RequestInit).method).toBe("PATCH");
    expect(JSON.parse((spy.mock.calls[0]![1] as RequestInit).body as string)).toEqual({ flight_path_m: 1.81 });
  });

  it("updateExperiment accepts name + description + file patterns", async () => {
    const spy = mockFetchSpy({ id: 7, name: "renamed" });
    const patch: api.ExperimentPatch = {
      name: "renamed", description: "AgBe run",
      image_pattern: "*.tif", metadata_pattern: "*.prp", integration_pattern: "*.dat",
    };
    await api.updateExperiment(7, patch, { username: "u" });
    const sent = JSON.parse((spy.mock.calls[0]![1] as RequestInit).body as string);
    expect(sent.name).toBe("renamed");
    expect(sent.image_pattern).toBe("*.tif");
  });

  it("suggestPaths GETs /api/fs/suggest?prefix=…", async () => {
    const spy = mockFetchSpy({ suggestions: ["/Volumes/data/ssrl/2026_04/1p7m"] });
    const out = await api.suggestPaths("/Volumes/data/ssrl/2026_04/1");
    expect(out.suggestions[0]).toContain("1p7m");
    expect(spy.mock.calls[0]![0]).toContain("/api/fs/suggest?");
    expect(spy.mock.calls[0]![0]).toContain("prefix=");
  });

  it("validatePath GETs /api/fs/validate?path=… and returns counts", async () => {
    const spy = mockFetchSpy({ ok: true, matched: 682, scanned: 700, message: null });
    const out = await api.validatePath("/Volumes/data/ssrl/2026_04/1p7m");
    expect(out.matched).toBe(682);
    expect(spy.mock.calls[0]![0]).toContain("/api/fs/validate?");
  });
});

function mockFetchSpy(body: unknown) {
  return vi.spyOn(globalThis, "fetch").mockResolvedValueOnce({
    ok: true, status: 200,
    text: async () => JSON.stringify(body),
    json: async () => body,
  } as Response);
}
```

- [ ] **Step 2: Run it red** — `npm test -- test/apiFetchers.ingestion.test.ts` → FAIL (none of these fetchers/types exist).

- [ ] **Step 3: Implement** — append to `src/api.ts` (after the Experiments block, ~line 127). Replace the existing `updateExperiment` (lines 123-127) and add the rest:

```ts
/** Canonical PATCH body for `PATCH /api/experiments/:id`. **E1 DEFINES,
 *  E2 IMPORTS — never redefine.** All fields are optional; the backend
 *  writes what is present:
 *  - name/description: plain writes, NO *_source stamp, NO rescan.
 *  - Geometry ×6: each field written + *_source stamped 'user' server-side
 *    (already built in Phase C — this widens the same route).
 *  - File patterns ×3: plain write + scan_signature invalidated server-side
 *    so the next scan re-discovers with the new glob.
 *  - data_dir/analysis_dir/path are READ-ONLY (400 if sent). */
export interface ExperimentPatch {
  name?: string;
  description?: string;
  energy_kev?: number;
  flight_path_m?: number;
  beam_center_x?: number;
  beam_center_y?: number;
  pixel_size_um?: number;
  q_units?: string;
  image_pattern?: string;
  metadata_pattern?: string;
  integration_pattern?: string;
}

export const updateExperiment = (
  id: number,
  patch: ExperimentPatch,
  opts?: AuthOpts,
) => request<Experiment>("PATCH", `/api/experiments/${id}`, patch, opts);

/** Create-from-directory (spec §9.2). Returns the new experiment id
 *  immediately; the first scan runs async with progress over SSE. */
export interface CreateExperimentBody {
  path: string;
  name?: string;
  patterns?: { image?: string; metadata?: string; integration?: string };
}
export const createExperiment = (body: CreateExperimentBody, opts?: AuthOpts) =>
  request<Experiment>("POST", "/api/experiments", body, opts);

/** Rescan: cheap change-check then additive ingest of new files. Idempotent. */
export const triggerScan = (id: number, opts?: AuthOpts) =>
  request<Experiment>("POST", `/api/experiments/${id}/scan`, {}, opts);

/** The Load ▸ Sample ▸ Exposures roll-up for the grouping-review surface
 *  (spec §9.2 — a dedicated endpoint, distinct from the flat corpus samples). */
export const listLoads = (id: number) =>
  request<Load[]>("GET", `/api/experiments/${id}/loads`);

/** Directory-picker path autocomplete (spec §9.2, read-only). */
export interface PathSuggestResponse { suggestions: string[] }
export const suggestPaths = (prefix: string) =>
  request<PathSuggestResponse>(
    "GET", `/api/fs/suggest?prefix=${encodeURIComponent(prefix)}`);

/** Directory-picker validate-path probe (spec §9.2). `matched`/`scanned` drive
 *  the validation line; `ok=false` + `message` powers the failed-scan preview. */
export interface ValidatePathResponse {
  ok: boolean;
  matched: number;
  scanned: number;
  message: string | null;
}
export const validatePath = (path: string) =>
  request<ValidatePathResponse>(
    "GET", `/api/fs/validate?path=${encodeURIComponent(path)}`);
```

- [ ] **Step 4: Run it green** — `npm test -- test/apiFetchers.ingestion.test.ts` → PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/api.ts packages/HimalayaUI/frontend/test/apiFetchers.ingestion.test.ts
git commit -m "feat(api): scan/ingest/loads/geometry-patch + directory-picker fetchers"
```

---

## Task 3: `state.ts` — ephemeral `ingestInFlight` slice

Per spec §9.6: an ephemeral `ingestInFlight: Record<number, IngestProgress> | null`, named setters, **NOT** in `partialize`, **NO** persist `version` bump (UI-only state). Also persist `activeExperimentId` is already there — no change needed.

**Files:**
- Modify: `src/state.ts` (`AppState` interface; the store body; named actions)
- Test: `test/stateIngest.test.ts` (CREATE)

- [ ] **Step 1: Write the failing test**

```ts
// test/stateIngest.test.ts
import { describe, it, expect, beforeEach } from "vitest";
import { useAppState } from "../src/state";

describe("ingestInFlight store slice (Phase E1)", () => {
  beforeEach(() => {
    useAppState.setState({ ingestInFlight: null });
  });

  it("starts null", () => {
    expect(useAppState.getState().ingestInFlight).toBeNull();
  });

  it("setIngestProgress upserts one experiment's progress", () => {
    useAppState.getState().setIngestProgress(7, { processed: 100, total: 680, status: "scanning" });
    expect(useAppState.getState().ingestInFlight).toEqual({
      7: { processed: 100, total: 680, status: "scanning" },
    });
    useAppState.getState().setIngestProgress(7, { processed: 300, total: 680, status: "scanning" });
    expect(useAppState.getState().ingestInFlight![7]!.processed).toBe(300);
  });

  it("clearIngestProgress removes one experiment, nulling the map when empty", () => {
    useAppState.getState().setIngestProgress(7, { processed: 1, total: 2, status: "scanning" });
    useAppState.getState().clearIngestProgress(7);
    expect(useAppState.getState().ingestInFlight).toBeNull();
  });

  it("is NOT included in the persisted partition", () => {
    // partialize (state.ts:500-508) must not surface ingestInFlight (ephemeral).
    // "himalaya-ui:state" is the real persist key — it is `LS_KEY` (state.ts:24),
    // the `name` passed to zustand persist (state.ts:489). Assert the serialized
    // blob never carries the key.
    useAppState.getState().setIngestProgress(7, { processed: 1, total: 2, status: "scanning" });
    const raw = localStorage.getItem("himalaya-ui:state");
    expect(raw === null || !raw.includes("ingestInFlight")).toBe(true);
  });
});
```

- [ ] **Step 2: Run it red** — `npm test -- test/stateIngest.test.ts` → FAIL (`setIngestProgress`/`ingestInFlight` undefined).

- [ ] **Step 3: Implement**

Add the type + interface fields near the other ephemeral state in `AppState` (after `previewIndexId`, ~line 83):

```ts
  /** Per-experiment live-ingest progress (spec §9.3/§9.6). Ephemeral — written
   *  by the App.tsx SSE listener on `ingest_*` frames, read by the experiment
   *  header + LiveIngestUnfold. NOT persisted (omitted from partialize); a
   *  reload re-derives terminal state from the experiment's `ingest_status`. */
  ingestInFlight: Record<number, IngestProgress> | null;
```

Add the `IngestProgress` type near the top of `state.ts` (after the `NavModalStep` type, ~line 37):

```ts
export type IngestProgressStatus = "scanning" | "analyzing" | "complete" | "failed";
export interface IngestProgress {
  processed: number;
  total: number;
  status: IngestProgressStatus;
}
```

Add the setter signatures to the interface (in the setters block, after `setPreviewIndex`, ~line 172):

```ts
  setIngestProgress: (experimentId: number, progress: IngestProgress) => void;
  clearIngestProgress: (experimentId: number) => void;
```

Add the initial value in the store body (after `previewIndexId: undefined,`, ~line 287):

```ts
        ingestInFlight: null,
```

Add the actions in the store body (after `setPreviewIndex`, ~line 335):

```ts
        setIngestProgress: (experimentId, progress) =>
          set((s) => ({
            ingestInFlight: { ...(s.ingestInFlight ?? {}), [experimentId]: progress },
          })),
        clearIngestProgress: (experimentId) =>
          set((s) => {
            if (s.ingestInFlight === null) return {};
            const next = { ...s.ingestInFlight };
            delete next[experimentId];
            return { ingestInFlight: Object.keys(next).length === 0 ? null : next };
          }),
```

> **Do NOT** add `ingestInFlight` to `partialize` (lines 500-508) and **do NOT** bump `version` (stays 5). The slice is ephemeral by design (spec §9.6: bumping `version` without a key-preserving migrate risks wiping the persisted blob).

- [ ] **Step 4: Run it green** — `npm test -- test/stateIngest.test.ts` → PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/state.ts packages/HimalayaUI/frontend/test/stateIngest.test.ts
git commit -m "feat(state): ephemeral ingestInFlight slice + setters (not persisted)"
```

---

## Task 4: `queries.ts` — `queryKeys` + read hooks for experiment detail / loads

Add `queryKeys.loads(id)` (distinct from `queryKeys.samples(experimentId)`, and shaped like it — `id ?? "none"`-tolerant), plus the read hook `useLoads`. (Mutation hooks — `useStartIngest`, `useSaveGeometryOverride`, the structural-edit hooks — are E2 / deferred; this plan adds only the read hook the skeleton pages need. The Configuration page calls `useExperiment` directly.)

> **Spec deviation — `queryKeys.experimentConfig` / `useExperimentConfig` / `useExperimentDetail` are NOT implemented (spec §9.6 names them).** The Configuration page reads from `useExperiment` (= `queryKeys.experiment(id)`, the same experiment detail row), not a separate `/config` resource and not an alias wrapper. There is no `/api/experiments/:id/config` endpoint in the as-built backend, and nothing invalidates a `config` key separately. Minting `queryKeys.experimentConfig` would produce a dead cache key that is never populated or invalidated. Decision: the Configuration page and ExperimentShell call `useExperiment` directly (same key, same `enabled: id > 0` gate). If E2 later requires a distinct config resource (e.g. pattern-set read from a new endpoint), it can mint the key then. The spec §9.6 references to `queryKeys.experimentConfig` and `useExperimentDetail` are superseded by this plan's scope.

**Files:**
- Modify: `src/queries.ts` (`queryKeys` object ~48-103; add hooks after `useExperiment` ~120)
- Test: `test/queriesIngestion.test.tsx` (CREATE)

- [ ] **Step 1: Write the failing test**

```tsx
// test/queriesIngestion.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { renderHook, waitFor } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { queryKeys, useLoads } from "../src/queries";
import * as api from "../src/api";

function wrapper() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return ({ children }: { children: React.ReactNode }) => (
    <QueryClientProvider client={qc}>{children}</QueryClientProvider>
  );
}

describe("ingestion queryKeys + hooks (Phase E1)", () => {
  beforeEach(() => vi.restoreAllMocks());
  afterEach(() => vi.restoreAllMocks());

  it("loads key is distinct from samples key, and `undefined`-tolerant like samples", () => {
    expect(queryKeys.loads(7)).toEqual(["experiment", 7, "loads"]);
    expect(queryKeys.samples(7)).toEqual(["experiment", 7, "samples"]);
    // Mirrors queryKeys.samples' `id ?? "none"` shape (queries.ts:51) so a
    // disabled (undefined-id) loads query never prefix-collides with an
    // enabled one.
    expect(queryKeys.loads(undefined)).toEqual(["experiment", "none", "loads"]);
  });

  it("useLoads fetches listLoads", async () => {
    const spy = vi.spyOn(api, "listLoads").mockResolvedValue([]);
    const { result } = renderHook(() => useLoads(7), { wrapper: wrapper() });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(spy).toHaveBeenCalledWith(7);
  });
});
```

- [ ] **Step 2: Run it red** — `npm test -- test/queriesIngestion.test.tsx` → FAIL (`queryKeys.loads`/`useLoads` undefined).

- [ ] **Step 3: Implement**

Add to the `queryKeys` object (after `samples:` ~line 52, mirroring its `id ?? "none"` signature so the disabled-query collision guard holds):

```ts
  loads:      (id: number | undefined) =>
    ["experiment", id ?? "none", "loads"] as const,
```

Add `useLoads` after `useExperiment` (~line 121). Callers that previously used `useExperimentDetail` call `useExperiment` directly — no alias is needed:

```ts
/** The Load ▸ Sample ▸ Exposures roll-up (grouping-review + Configuration
 *  acquisition timeline read off this; spec §9.2/§9.6). Gated id>0 like
 *  useExperiment so a zero/undefined experiment never hits /loads. */
export function useLoads(id: number) {
  return useQuery({
    queryKey: queryKeys.loads(id),
    queryFn: () => api.listLoads(id),
    enabled: id > 0,
  });
}
```

- [ ] **Step 4: Run it green** — `npm test -- test/queriesIngestion.test.tsx` → PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/queries.ts packages/HimalayaUI/frontend/test/queriesIngestion.test.tsx
git commit -m "feat(queries): loads key + useLoads hook"
```

---

## Task 5: `Dropdown` primitive (`src/print/ui/Dropdown.tsx`)

A reusable labelled select-trigger + Menu popover. Built on the existing `Menu<T>` primitive (it owns the popover list look + roving focus). `Dropdown` adds the trigger button (current-value display + chevron), open/close state, outside-click close, and focus-return to the trigger on close (which `Menu` leaves to the owner). Consumed by `DirectoryPickerField` (Task 8) and exported for E2's grouping/sources menus. Appearance lives here (`src/print/ui/`); consumer `className` is placement-only.

**Files:**
- Create: `src/print/ui/Dropdown.tsx`
- Modify: `src/print/ui/index.ts` (export)
- Test: `test/Dropdown.test.tsx` (CREATE)

- [ ] **Step 1: Write the failing test**

```tsx
// test/Dropdown.test.tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { Dropdown } from "../src/print/ui/Dropdown";

const OPTIONS = [
  { value: "a", label: "Alpha" },
  { value: "b", label: "Beta" },
] as const;

describe("Dropdown (Phase E1 primitive)", () => {
  it("renders the active option's label on the trigger", () => {
    render(<Dropdown aria-label="Pick" options={OPTIONS} value="b" onChange={() => {}} />);
    expect(screen.getByTestId("dropdown-trigger")).toHaveTextContent("Beta");
  });

  it("opens the menu on trigger click and selects an option", () => {
    const onChange = vi.fn();
    render(<Dropdown aria-label="Pick" options={OPTIONS} value="a" onChange={onChange} />);
    expect(screen.queryByRole("menu")).toBeNull();
    fireEvent.click(screen.getByTestId("dropdown-trigger"));
    expect(screen.getByRole("menu")).toBeInTheDocument();
    fireEvent.click(screen.getByText("Beta"));
    expect(onChange).toHaveBeenCalledWith("b");
    // closes after select
    expect(screen.queryByRole("menu")).toBeNull();
  });

  it("renders a placeholder when value matches no option", () => {
    render(
      <Dropdown aria-label="Pick" options={OPTIONS} value={undefined}
        placeholder="Choose…" onChange={() => {}} />,
    );
    expect(screen.getByTestId("dropdown-trigger")).toHaveTextContent("Choose…");
  });

  it("trigger carries aria-haspopup=menu and aria-expanded", () => {
    render(<Dropdown aria-label="Pick" options={OPTIONS} value="a" onChange={() => {}} />);
    const trigger = screen.getByTestId("dropdown-trigger");
    expect(trigger).toHaveAttribute("aria-haspopup", "menu");
    expect(trigger).toHaveAttribute("aria-expanded", "false");
    fireEvent.click(trigger);
    expect(trigger).toHaveAttribute("aria-expanded", "true");
  });
});
```

- [ ] **Step 2: Run it red** — `npm test -- test/Dropdown.test.tsx` → FAIL (no `Dropdown`).

- [ ] **Step 3: Implement** — create `src/print/ui/Dropdown.tsx`:

```tsx
import { useRef, useState } from "react";
import type { ReactNode } from "react";
import { Menu } from "./Menu";
import type { MenuOption } from "./Menu";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface DropdownProps<T extends string> {
  /** Required: names the control for assistive tech. */
  "aria-label": string;
  options: ReadonlyArray<MenuOption<T>>;
  /** The current value. `undefined` shows the placeholder. */
  value: T | undefined;
  onChange: (value: T) => void;
  /** Shown on the trigger when no option matches `value`. */
  placeholder?: ReactNode;
  /** Disables the trigger. */
  disabled?: boolean;
  /** PLACEMENT-ONLY: margin / width / position. No appearance utils. */
  className?: string;
}

/**
 * Dropdown<T> — a labelled select-trigger that opens a {@link Menu} value-
 * selector popover (closed look / open placement). The trigger shows the active
 * option's label (or `placeholder`) + a chevron; clicking toggles the Menu,
 * which owns the list look + roving focus + Escape. Selecting closes the menu
 * and returns focus to the trigger (Menu leaves focus-return to the owner).
 *
 * The reusable replacement for ad-hoc `<select>` + `Menu` wirings (spec: "turn
 * the dropdown into a component and reuse it"). Consumed by DirectoryPickerField
 * here and the grouping/sources menus in E2.
 */
export function Dropdown<T extends string>({
  "aria-label": ariaLabel,
  options,
  value,
  onChange,
  placeholder,
  disabled = false,
  className = "",
}: DropdownProps<T>): JSX.Element {
  const [open, setOpen] = useState(false);
  const triggerRef = useRef<HTMLButtonElement>(null);

  const active = options.find((o) => o.value === value);
  const triggerLabel = active ? active.label : (placeholder ?? "");

  const close = (): void => {
    setOpen(false);
    triggerRef.current?.focus();
  };

  return (
    <div className={cx("relative inline-flex", className)}>
      <button
        ref={triggerRef}
        type="button"
        data-testid="dropdown-trigger"
        aria-haspopup="menu"
        aria-expanded={open}
        aria-label={ariaLabel}
        disabled={disabled}
        onClick={() => setOpen((v) => !v)}
        className={cx(
          "inline-flex items-center justify-between gap-2 min-w-0",
          "rounded-sm border border-hair-strong bg-plate px-2.5 py-1.5",
          "text-sm text-ink transition-colors hover:border-accent",
          "focus-visible:outline focus-visible:outline-2 focus-visible:outline-accent focus-visible:outline-offset-2",
          "disabled:opacity-40 disabled:cursor-not-allowed",
        )}
      >
        <span className={cx("truncate", active ? "text-ink" : "text-ink-soft")}>
          {triggerLabel}
        </span>
        <span aria-hidden="true" className="text-ink-faint">▾</span>
      </button>
      <Menu<T>
        open={open}
        options={options}
        activeValue={value}
        onSelect={onChange}
        onClose={close}
        aria-label={ariaLabel}
        className="left-0 top-full w-full"
      />
    </div>
  );
}
```

> Appearance (border/bg/text/radius/focus ring) lives here — this file is `src/print/ui/**`, exempt from the design guard. `rounded-sm` is the one 5px step. The Menu is anchored by the `relative` wrapper; `left-0 top-full w-full` is placement-only on the Menu.

- [ ] **Step 4: Export it** — add to `src/print/ui/index.ts` (alphabetically near `Dot`/`EmptyState`):

```ts
export { Dropdown } from "./Dropdown";
export type { DropdownProps } from "./Dropdown";
```

- [ ] **Step 5: Run it green** — `npm test -- test/Dropdown.test.tsx` → PASS.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/ui/Dropdown.tsx packages/HimalayaUI/frontend/src/print/ui/index.ts packages/HimalayaUI/frontend/test/Dropdown.test.tsx
git commit -m "feat(ui): reusable Dropdown primitive (trigger + Menu popover)"
```

---

## Task 6: `StatBar` primitive (`src/print/ui/StatBar.tsx`)

The hairline-divided mono stat band (`loads · samples · exposures · sessions`/date-range). Per Jonathan's note, its own component. Each cell shows a mono value + an uppercase caption; a thin left hairline divides cells (first cell has none). A cell may be `pending` (italic faint placeholder) for the live-ingest "discovered-so-far" state (spec §8.3).

**Files:**
- Create: `src/print/ui/StatBar.tsx`
- Modify: `src/print/ui/index.ts`
- Test: `test/StatBar.test.tsx` (CREATE)

- [ ] **Step 1: Write the failing test**

```tsx
// test/StatBar.test.tsx
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { StatBar } from "../src/print/ui/StatBar";

describe("StatBar (Phase E1 primitive)", () => {
  it("renders one cell per stat with value + caption", () => {
    render(
      <StatBar
        aria-label="Experiment stats"
        stats={[
          { key: "loads", caption: "Loads", value: "13" },
          { key: "samples", caption: "Samples", value: "170" },
          { key: "exposures", caption: "Exposures", value: "682" },
          { key: "sessions", caption: "Sessions", value: "4" },
        ]}
      />,
    );
    const cells = screen.getAllByTestId("statbar-cell");
    expect(cells).toHaveLength(4);
    expect(cells[0]).toHaveTextContent("Loads");
    expect(cells[0]).toHaveTextContent("13");
  });

  it("marks a pending cell for assistive tech + a data flag", () => {
    render(
      <StatBar
        aria-label="Experiment stats"
        stats={[{ key: "span", caption: "Span", value: "pending", pending: true }]}
      />,
    );
    expect(screen.getByTestId("statbar-cell")).toHaveAttribute("data-pending", "true");
  });

  it("exposes the aria-label on the band", () => {
    render(<StatBar aria-label="X" stats={[{ key: "a", caption: "A", value: "1" }]} />);
    expect(screen.getByTestId("statbar")).toHaveAttribute("aria-label", "X");
  });
});
```

- [ ] **Step 2: Run it red** — `npm test -- test/StatBar.test.tsx` → FAIL (no `StatBar`).

- [ ] **Step 3: Implement** — create `src/print/ui/StatBar.tsx`:

```tsx
import type { ReactNode } from "react";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface StatBarStat {
  /** Stable react key + test handle. */
  key: string;
  /** Uppercase micro-caption under the value. */
  caption: string;
  /** The value (already formatted by the caller; mono unless pending). */
  value: ReactNode;
  /** When true, the value renders as a faint italic placeholder
   *  ("discovered-so-far" / span-pending; spec §8.3). */
  pending?: boolean;
}

export interface StatBarProps {
  /** Required: names the band for assistive tech. */
  "aria-label": string;
  stats: ReadonlyArray<StatBarStat>;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/**
 * StatBar — the quiet stat ledger (DESIGN: hairline-divided cells, no heavy
 * rule). Each cell is a mono value over an uppercase micro-caption; a thin left
 * hairline divides cells (first has none). A `pending` cell italicizes a faint
 * placeholder for the live-ingest discovered-so-far state. Closed look; the
 * consumer's `className` is placement-only.
 */
export function StatBar({
  "aria-label": ariaLabel,
  stats,
  className = "",
}: StatBarProps): JSX.Element {
  return (
    <div
      data-testid="statbar"
      role="group"
      aria-label={ariaLabel}
      className={cx("flex", className)}
    >
      {stats.map((s, i) => (
        <div
          key={s.key}
          data-testid="statbar-cell"
          data-pending={s.pending ? "true" : undefined}
          className={cx(
            "flex flex-col gap-1 px-7 py-0.5",
            i === 0 ? "pl-0" : "border-l border-hair",
          )}
        >
          <span
            className={cx(
              "leading-none",
              s.pending
                ? "text-sm italic text-ink-faint"
                : "font-mono text-xl text-ink",
            )}
          >
            {s.value}
          </span>
          <span className="text-xs font-bold uppercase tracking-wide text-ink-faint">
            {s.caption}
          </span>
        </div>
      ))}
    </div>
  );
}
```

- [ ] **Step 4: Export it** — add to `src/print/ui/index.ts`:

```ts
export { StatBar } from "./StatBar";
export type { StatBarProps, StatBarStat } from "./StatBar";
```

- [ ] **Step 5: Run it green** — `npm test -- test/StatBar.test.tsx` → PASS.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/ui/StatBar.tsx packages/HimalayaUI/frontend/src/print/ui/index.ts packages/HimalayaUI/frontend/test/StatBar.test.tsx
git commit -m "feat(ui): StatBar primitive (hairline-divided stat ledger)"
```

---

## Task 7: `PageFrame` — add `home` and `experiment` width keys

Per spec §9.6: add `home` (~1080px) and `experiment` (~1280px) keys to `PAGE_WIDTHS`.

**Files:**
- Modify: `src/print/components/PageFrame.tsx` (`PAGE_WIDTHS` ~7-13)
- Test: `test/PageFrame.ingestion.test.tsx` (CREATE)

- [ ] **Step 1: Write the failing test**

```tsx
// test/PageFrame.ingestion.test.tsx
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { PageFrame, PAGE_WIDTHS } from "../src/print/components/PageFrame";

describe("PageFrame ingestion widths (Phase E1)", () => {
  it("exposes home + experiment width keys", () => {
    expect(PAGE_WIDTHS.home).toBeDefined();
    expect(PAGE_WIDTHS.experiment).toBeDefined();
  });
  it("renders with the new keys", () => {
    render(<PageFrame width="experiment">x</PageFrame>);
    expect(screen.getByTestId("page-frame")).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run it red** — `npm test -- test/PageFrame.ingestion.test.tsx` → FAIL (`home`/`experiment` keys missing; the `PAGE_WIDTHS.home` access is `undefined`).

- [ ] **Step 3: Implement** — extend the `PAGE_WIDTHS` map in `src/print/components/PageFrame.tsx`:

```ts
export const PAGE_WIDTHS = {
  loupe: "max-w-[1100px]",
  sheet: "max-w-[1240px]",
  folio: "max-w-[1380px]",
  scoping: "max-w-[760px]",
  builder: "max-w-[1180px]",
  home: "max-w-[1080px]",
  experiment: "max-w-[1280px]",
} as const;
```

- [ ] **Step 4: Run it green** — `npm test -- test/PageFrame.ingestion.test.tsx` → PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/components/PageFrame.tsx packages/HimalayaUI/frontend/test/PageFrame.ingestion.test.tsx
git commit -m "feat(ui): PageFrame home + experiment width keys"
```

---

## Task 8: `DirectoryPickerField` (`src/print/components/DirectoryPickerField.tsx`)

The directory-path picker: a path `Input` + a live suggestion dropdown (built on `Popover` + the suggestion list), Tab-to-complete the top suggestion, ↑/↓ to choose, and a validation line driven by `validatePath`. **No full file browser** (spec §8.7). It is a composite (composes `ui/` primitives), so it lives under `src/print/components/`; it carries no inline appearance utilities (appearance comes from `Input`/`Popover`/tokens). Slots into `NewExperimentPage` (Task 12).

**Files:**
- Create: `src/print/components/DirectoryPickerField.tsx`
- Test: `test/DirectoryPickerField.test.tsx` (CREATE)

- [ ] **Step 1: Write the failing test**

```tsx
// test/DirectoryPickerField.test.tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { DirectoryPickerField } from "../src/print/components/DirectoryPickerField";

describe("DirectoryPickerField (Phase E1)", () => {
  it("shows suggestions and Tab completes the top one", () => {
    const onChange = vi.fn();
    render(
      <DirectoryPickerField
        value="/Volumes/data/ssrl/2026_04/1"
        onChange={onChange}
        suggestions={["/Volumes/data/ssrl/2026_04/1p7m", "/Volumes/data/ssrl/2026_04/10x"]}
        validation={null}
      />,
    );
    // Input puts data-testid on its WRAPPER; onKeyDown (via ...rest) lands on
    // the inner <input>. Fire on the inner input so the handler runs.
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.keyDown(input, { key: "Tab" });
    expect(onChange).toHaveBeenCalledWith("/Volumes/data/ssrl/2026_04/1p7m");
  });

  it("arrow-down moves the active suggestion, Enter completes it", () => {
    const onChange = vi.fn();
    render(
      <DirectoryPickerField
        value="/Volumes/data/ssrl/2026_04/1"
        onChange={onChange}
        suggestions={["/a/one", "/a/two"]}
        validation={null}
      />,
    );
    const input = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.keyDown(input, { key: "ArrowDown" }); // moves to index 1
    fireEvent.keyDown(input, { key: "Enter" });
    expect(onChange).toHaveBeenCalledWith("/a/two");
  });

  it("renders a positive validation line", () => {
    render(
      <DirectoryPickerField
        value="/d" onChange={() => {}} suggestions={[]}
        validation={{ ok: true, matched: 682, scanned: 700, message: null }}
      />,
    );
    expect(screen.getByTestId("dirpicker-validation")).toHaveTextContent("682");
  });

  it("renders a failure validation line with the message", () => {
    render(
      <DirectoryPickerField
        value="/d" onChange={() => {}} suggestions={[]}
        validation={{ ok: false, matched: 0, scanned: 0, message: "No exposures found" }}
      />,
    );
    const line = screen.getByTestId("dirpicker-validation");
    expect(line).toHaveAttribute("data-ok", "false");
    expect(line).toHaveTextContent("No exposures found");
  });
});
```

- [ ] **Step 2: Run it red** — `npm test -- test/DirectoryPickerField.test.tsx` → FAIL (no component).

- [ ] **Step 3: Implement** — create `src/print/components/DirectoryPickerField.tsx`:

```tsx
import { useRef, useState } from "react";
import type { KeyboardEvent } from "react";
import { Input } from "../ui/Input";
import { Popover } from "../ui/Popover";
import type { ValidatePathResponse } from "../../api";

export interface DirectoryPickerFieldProps {
  value: string;
  onChange: (path: string) => void;
  /** Live autocomplete suggestions for the current `value` (caller fetches via
   *  `suggestPaths`). */
  suggestions: ReadonlyArray<string>;
  /** Live validation probe result, or null when not yet probed. */
  validation: ValidatePathResponse | null;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/**
 * DirectoryPickerField — path input + suggestion list (via Popover) +
 * Tab/↑↓/↵ completion + a validation line. No full file browser (spec §8.7).
 * Composes `Input`, `Popover`, and design-system token utilities only
 * (no inline appearance); the suggestion list is anchored inside Popover,
 * which provides outside-pointerdown dismiss, Escape close, and focus-return
 * to the trigger — the functional gaps of a hand-rolled `<ul role=listbox>`.
 *
 * Popover trigger: the Input itself (its inner `<input>` receives the
 * `onValueChange` + keyboard handlers via `...rest`; the Popover clone adds
 * `aria-haspopup="dialog"` + `aria-expanded` without conflicting with those).
 * The popover stays open whenever `suggestions.length > 0` — Popover's
 * outside-pointerdown close matches the user expectation of clicking away to
 * dismiss. Escape resets focus to the input (Popover handles it).
 */
export function DirectoryPickerField({
  value,
  onChange,
  suggestions,
  validation,
  className = "",
}: DirectoryPickerFieldProps): JSX.Element {
  const [active, setActive] = useState(0);
  const inputFocusRef = useRef<HTMLInputElement>(null);

  const complete = (path: string): void => {
    onChange(path);
    setActive(0);
  };

  const onKeyDown = (e: KeyboardEvent<HTMLInputElement>): void => {
    if (suggestions.length === 0) return;
    if (e.key === "Tab") {
      // Tab completes the TOP suggestion (spec §8.7). preventDefault so focus
      // doesn't leave the field on the completing keystroke.
      e.preventDefault();
      complete(suggestions[0]!);
    } else if (e.key === "ArrowDown") {
      e.preventDefault();
      setActive((i) => Math.min(i + 1, suggestions.length - 1));
    } else if (e.key === "ArrowUp") {
      e.preventDefault();
      setActive((i) => Math.max(i - 1, 0));
    } else if (e.key === "Enter") {
      e.preventDefault();
      complete(suggestions[active] ?? suggestions[0]!);
    }
  };

  // The Input is the Popover trigger. Popover's cloneElement adds
  // aria-haspopup/aria-expanded/onClick to it; onValueChange and onKeyDown
  // pass through Input's ...rest spread onto the inner <input>. The Popover
  // open state is driven externally by whether suggestions is non-empty —
  // we control `open` to avoid Popover's click-toggle interfering with the
  // keyboard-driven flow (suggestions appear while typing, not on click).
  // Popover's outside-pointerdown effect and Escape handler both fire when open.
  const hasSuggestions = suggestions.length > 0;

  return (
    <div className={`flex flex-col gap-2 ${className}`.trim()}>
      <Popover
        trigger={
          <Input
            testId="dirpicker-input"
            value={value}
            onValueChange={onChange}
            onKeyDown={onKeyDown}
            mono
            placeholder="/Volumes/data/ssrl/2026_04/…"
            aria-label="Data directory"
          />
        }
        label="Directory suggestions"
        className="w-full"
      >
        {hasSuggestions && (
          <ul
            data-testid="dirpicker-suggestions"
            role="listbox"
            aria-label="Directory suggestions"
            className="py-1"
          >
            {suggestions.map((s, i) => (
              <li key={s}>
                <button
                  type="button"
                  role="option"
                  aria-selected={i === active}
                  data-active={i === active ? "true" : undefined}
                  onClick={() => complete(s)}
                  className={
                    "flex w-full px-3 py-1.5 text-left font-mono text-sm transition-colors " +
                    (i === active ? "text-ink bg-paper-sunk" : "text-ink-soft hover:text-ink hover:bg-paper-sunk")
                  }
                >
                  {s}
                </button>
              </li>
            ))}
          </ul>
        )}
      </Popover>
      {validation && (
        <div
          data-testid="dirpicker-validation"
          data-ok={validation.ok ? "true" : "false"}
          className={"text-sm " + (validation.ok ? "text-ink-soft" : "text-error")}
        >
          {validation.ok
            ? `Matched ${validation.matched} of ${validation.scanned} files`
            : (validation.message ?? "No exposures found")}
        </div>
      )}
    </div>
  );
}
```

> `text-error` is a token utility (`--color-error` exists in `styles.css:28`; `Input` itself uses the `border-error` family). The suggestion list + validation line use only token colour utilities (`bg-plate`, `text-ink*`, `bg-paper-sunk`, `border-hair`) and named scale roles (`text-sm`) — no inline appearance literals, so the design guard passes. `rounded-sm` is the 5px step. **`onKeyDown` flows through `Input`'s `...rest` onto the inner `<input>`** (confirmed `Input.tsx:84-95`), so the test fires keydown on `getByTestId("dirpicker-input").querySelector("input")`, not the wrapper. **`Popover` provides the outside-pointerdown close, Escape close, and focus-return that the previous hand-rolled `<ul role=listbox>` lacked** (Popover.tsx:58-68 confirmed — the `useEffect` binds `pointerdown` to `document` only while open, and Escape fires `close()` which calls `triggerRef.current?.focus()`). Note: `Popover`'s `trigger` clone adds `onClick` (toggle) onto the Input — in jsdom this is harmless since the suggestions panel is conditionally rendered by `hasSuggestions`, but in production the Input receives a click-toggle in addition to its normal focus/change handlers. If this proves intrusive, the Popover open state can be lifted to a `useState(false)` driven by Input `onFocus`/`onBlur` — but test the simple composition first.

- [ ] **Step 4: Run it green** — `npm test -- test/DirectoryPickerField.test.tsx` → PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/components/DirectoryPickerField.tsx packages/HimalayaUI/frontend/test/DirectoryPickerField.test.tsx
git commit -m "feat(components): DirectoryPickerField (path input + suggestions + Tab-complete + validation)"
```

---

## Task 9: `ExperimentTopNav` (`src/print/shell/ExperimentTopNav.tsx`)

The `Experiments | Series` top nav used by `ExperimentShell` (which renders its own chrome outside `CorpusShell`). Router `<Link>`s with `pathname.startsWith` active detection, mirroring `CorpusTopbar`'s stage-tab pattern (tokens only — `text-ink`/`text-ink-soft`/`bg-paper-sunk`/`bg-print-accent`/`rounded`). The home phase bar is intentionally absent (Jonathan's round-3 note: drop the home phase bar).

**Files:**
- Create: `src/print/shell/ExperimentTopNav.tsx`
- Test: `test/ExperimentTopNav.test.tsx` (CREATE)

- [ ] **Step 1: Write the failing test**

```tsx
// test/ExperimentTopNav.test.tsx
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { ExperimentTopNav } from "../src/print/shell/ExperimentTopNav";

function at(path: string) {
  render(
    <MemoryRouter initialEntries={[path]}>
      <ExperimentTopNav />
    </MemoryRouter>,
  );
}

describe("ExperimentTopNav (Phase E1)", () => {
  it("links to /experiments and /series", () => {
    at("/experiments");
    expect(screen.getByTestId("topnav-experiments")).toHaveAttribute("href", "/experiments");
    expect(screen.getByTestId("topnav-series")).toHaveAttribute("href", "/series");
  });

  it("marks Experiments active on an experiment route", () => {
    at("/experiments/7/corpus");
    expect(screen.getByTestId("topnav-experiments")).toHaveAttribute("aria-current", "page");
    expect(screen.getByTestId("topnav-series")).not.toHaveAttribute("aria-current");
  });

  it("marks Series active on /series", () => {
    at("/series");
    expect(screen.getByTestId("topnav-series")).toHaveAttribute("aria-current", "page");
  });
});
```

- [ ] **Step 2: Run it red** — `npm test -- test/ExperimentTopNav.test.tsx` → FAIL (no component).

- [ ] **Step 3: Implement** — create `src/print/shell/ExperimentTopNav.tsx`:

```tsx
import { Link, useLocation } from "react-router-dom";

interface NavItem {
  id: "experiments" | "series";
  label: string;
  to: string;
}

const ITEMS: readonly NavItem[] = [
  { id: "experiments", label: "Experiments", to: "/experiments" },
  { id: "series", label: "Series", to: "/series" },
];

/**
 * ExperimentTopNav — the two top-level axes (spec §7: Experiments | Series).
 * Used by ExperimentShell, which renders its own chrome OUTSIDE CorpusShell so
 * the experiment header never stacks on CorpusTopbar. Router Links with
 * `startsWith` active detection (mirrors CorpusTopbar's stage tabs). The home
 * phase bar is intentionally dropped (round-3 note).
 */
export function ExperimentTopNav(): JSX.Element {
  const { pathname } = useLocation();
  return (
    <nav data-testid="experiment-top-nav" aria-label="Sections" className="flex gap-0.5">
      {ITEMS.map((it) => {
        const isActive = pathname.startsWith(it.to);
        return (
          <Link
            key={it.id}
            to={it.to}
            data-testid={`topnav-${it.id}`}
            aria-current={isActive ? "page" : undefined}
            className={
              "px-2.5 py-1.5 rounded text-xs font-semibold uppercase tracking-wide no-underline " +
              (isActive ? "text-ink bg-paper-sunk" : "text-ink-soft")
            }
          >
            <span
              aria-hidden="true"
              className="inline-block w-1 h-1 rounded-full mr-1.5 align-middle bg-print-accent"
            />
            {it.label}
          </Link>
        );
      })}
    </nav>
  );
}
```

> Mirrors CorpusTopbar's stage-tab token usage exactly (`text-ink`/`bg-paper-sunk`/`bg-print-accent`/`rounded` — all `@theme` tokens, not appearance literals), so the design guard passes for a shell file.

- [ ] **Step 4: Run it green** — `npm test -- test/ExperimentTopNav.test.tsx` → PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/shell/ExperimentTopNav.tsx packages/HimalayaUI/frontend/test/ExperimentTopNav.test.tsx
git commit -m "feat(shell): ExperimentTopNav (Experiments | Series axes)"
```

---

## Task 10: `ExperimentShell` (`src/print/shell/ExperimentShell.tsx`)

The `/experiments/:id` layout route. Renders its own top chrome (wordmark + `ExperimentTopNav` + experiment header with edit-in-place name/description, rescan/ingest status, `StatBar` ledger) and the `Corpus | Configuration` tab bar above an `<Outlet>`. Sits **outside** `CorpusShell` (spec §9.6). Pulls the experiment from `useExperiment(:id)` (called directly — no alias); rescan/ingest status reads `experiment.ingest_status` + `ingestInFlight[id]`. Name/description edit-in-place uses `Input variant="title"` wired to `updateExperiment` — for E1 the inline edit is a controlled local draft that calls `updateExperiment` on commit (`PATCH /api/experiments/:id` is widened by Task 1c — `name` and `description` are live-writable on the real backend; the test mocks the API). **Do NOT add a DirectoryPickerField or path autocomplete to this shell** — the data directory is fixed at creation; spec §13's directory-edit / FS-sandboxing open question is moot for v1.

**Files:**
- Create: `src/print/shell/ExperimentShell.tsx`
- Test: `test/ExperimentShell.test.tsx` (CREATE)

- [ ] **Step 1: Write the failing test**

```tsx
// test/ExperimentShell.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ExperimentShell } from "../src/print/shell/ExperimentShell";
import * as api from "../src/api";

function renderAt(path: string) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[path]}>
        <Routes>
          <Route path="/experiments/:id" element={<ExperimentShell />}>
            <Route path="corpus" element={<div>CORPUS BODY</div>} />
            <Route path="config" element={<div>CONFIG BODY</div>} />
          </Route>
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

const EXP: api.Experiment = {
  id: 7, name: "SSRL · 1p7m", description: null, path: "/d", data_dir: "/d/data", analysis_dir: "/d/an",
  manifest_path: null, created_at: "2026-04-12", q_units: "A^-1",
  beam_center_x: 1, beam_center_y: 1, pixel_size_um: 172, energy_kev: 9, flight_path_m: 1.8,
  energy_kev_source: "prp", flight_path_m_source: "setup", beam_center_x_source: "setup",
  beam_center_y_source: "setup", pixel_size_um_source: "prp", q_units_source: "prp",
  last_scanned_at: "2026-04-12T10:00:00", scan_signature: "sig", ingest_status: "complete",
  image_pattern: null, metadata_pattern: null, integration_pattern: null,
};

describe("ExperimentShell (Phase E1)", () => {
  beforeEach(() => { vi.spyOn(api, "getExperiment").mockResolvedValue(EXP); });
  afterEach(() => vi.restoreAllMocks());

  it("renders own chrome (top nav) + the experiment header name", async () => {
    renderAt("/experiments/7/corpus");
    expect(screen.getByTestId("experiment-top-nav")).toBeInTheDocument();
    // The header name is an Input variant='title' (edit-in-place). Its wrapper
    // div carries data-testid; the inner <input> carries the value attribute.
    const nameWrapper = await screen.findByTestId("experiment-header-name");
    expect(nameWrapper.querySelector("input")?.value).toBe("SSRL · 1p7m");
  });

  it("renders the Corpus | Configuration tab bar with the active tab", async () => {
    renderAt("/experiments/7/corpus");
    await screen.findByTestId("experiment-header-name");
    expect(screen.getByTestId("exp-tab-corpus")).toHaveAttribute("aria-current", "page");
    expect(screen.getByTestId("exp-tab-config")).not.toHaveAttribute("aria-current");
  });

  it("renders the child route via Outlet", async () => {
    renderAt("/experiments/7/config");
    expect(await screen.findByText("CONFIG BODY")).toBeInTheDocument();
    expect(screen.getByTestId("exp-tab-config")).toHaveAttribute("aria-current", "page");
  });

  it("does NOT render the corpus topbar (it lives outside CorpusShell)", async () => {
    renderAt("/experiments/7/corpus");
    await screen.findByTestId("experiment-header-name");
    expect(screen.queryByTestId("corpus-topbar")).toBeNull();
  });
});
```

- [ ] **Step 2: Run it red** — `npm test -- test/ExperimentShell.test.tsx` → FAIL (no component).

- [ ] **Step 3: Implement** — create `src/print/shell/ExperimentShell.tsx`:

```tsx
import { Link, Outlet, useLocation, useParams } from "react-router-dom";
import { useState } from "react";
import { useExperiment } from "../../queries";
import { useAppState } from "../../state";
import * as api from "../../api";
import { authOpts } from "../../lib/authOpts";
import { Wordmark } from "../ui/Wordmark";
import { Kicker } from "../ui/Kicker";
import { Input } from "../ui/Input";
import { StatBar } from "../ui/StatBar";
import type { StatBarStat } from "../ui/StatBar";
import { ProgressBar } from "../ui/ProgressBar";
import { PageFrame } from "../components/PageFrame";
import { ExperimentTopNav } from "./ExperimentTopNav";

interface TabDef {
  id: "corpus" | "config";
  label: string;
}
const TABS: readonly TabDef[] = [
  { id: "corpus", label: "Corpus" },
  { id: "config", label: "Configuration" },
];

/**
 * ExperimentShell — the /experiments/:id layout route. Renders its OWN chrome
 * (wordmark + ExperimentTopNav + experiment header + Corpus|Configuration tab
 * bar) outside CorpusShell so the header never stacks on CorpusTopbar (spec
 * §7/§9.6). Children (ExperimentCorpusPage / ConfigurationPage) render in the
 * Outlet. The grouping-review route reuses this shell too but hides the tab
 * bar's active state (E2 wires the banner → grouping surface).
 */
export function ExperimentShell(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  const expId = id ? Number(id) : 0;
  const { pathname } = useLocation();
  const exp = useExperiment(expId);
  const inFlight = useAppState((s) => s.ingestInFlight?.[expId]);
  const username = useAppState((s) => s.username);

  // Controlled draft for the edit-in-place name field. Initialized from the
  // loaded experiment; resets when the server data changes (e.g. after a
  // commit from another tab). The commit handler calls updateExperiment PATCH
  // (Task 1c widens the route to accept name as a plain-write field).
  const serverName = exp.data?.name ?? "";
  const [nameDraft, setNameDraft] = useState(serverName);
  // Sync draft when server data first loads or changes externally.
  // (useState initializer only runs once; use an effect for subsequent changes.)
  const [lastServerId, setLastServerId] = useState(expId);
  if (expId !== lastServerId) {
    setLastServerId(expId);
    setNameDraft(exp.data?.name ?? "");
  }

  const commitName = (): void => {
    const trimmed = nameDraft.trim();
    if (trimmed === serverName) return; // no-op
    // authOpts omits undefined keys (exactOptionalPropertyTypes: true).
    void api.updateExperiment(expId, { name: trimmed }, authOpts(username, undefined));
  };

  const name = serverName || `Experiment ${expId}`;
  const status = inFlight?.status ?? exp.data?.ingest_status ?? "idle";
  const isProcessing = status === "scanning" || status === "analyzing";

  const stats: StatBarStat[] = isProcessing
    ? [
        { key: "processed", caption: "Processed",
          value: inFlight ? `${inFlight.processed} / ~${inFlight.total}` : "—" },
        { key: "span", caption: "Span", value: "pending", pending: true },
      ]
    : [
        // E1 ledger reads the experiment detail. Real counts (loads/samples/
        // exposures/sessions) arrive from the extended GET /api/experiments/:id
        // roll-up (spec §9.2); until then show placeholders the rollup fills.
        { key: "loads", caption: "Loads", value: "—" },
        { key: "samples", caption: "Samples", value: "—" },
        { key: "exposures", caption: "Exposures", value: "—" },
        { key: "sessions", caption: "Sessions", value: "—" },
      ];

  return (
    <div
      data-testid="experiment-shell"
      className="h-full w-full flex flex-col min-h-0 bg-paper text-ink overflow-auto"
    >
      {/* Own top chrome (NOT CorpusTopbar). */}
      <header className="flex items-center gap-6 px-6 h-14 border-b border-hair shrink-0">
        <Link to="/experiments" className="rounded focus-visible:outline focus-visible:outline-2 focus-visible:outline-accent focus-visible:outline-offset-2">
          <Wordmark tail="SAXS">Himalaya</Wordmark>
        </Link>
        <ExperimentTopNav />
      </header>

      <PageFrame width="experiment" className="px-6 py-6 flex-1 min-h-0 flex flex-col">
        {/* Experiment header */}
        <div className="flex items-start justify-between gap-6">
          <div className="min-w-0">
            <Kicker>Experiment</Kicker>
            {/* Edit-in-place name: Input variant='title' controlled by a local
                draft; commits on blur / Enter via updateExperiment PATCH (Task
                1c widens the route to accept name as a plain-write field). The
                test spies on api.getExperiment so experiment-header-name is
                populated by the resolved exp.data.name. */}
            <Input
              variant="title"
              testId="experiment-header-name"
              value={nameDraft || name}
              onValueChange={setNameDraft}
              onBlur={commitName}
              onKeyDown={(e) => { if (e.key === "Enter") { commitName(); (e.target as HTMLElement).blur(); } }}
              aria-label="Experiment name"
            />
            {exp.data?.data_dir && (
              <p className="text-sm text-ink-soft font-mono mt-1 truncate">{exp.data.data_dir}</p>
            )}
          </div>
          <div data-testid="experiment-rescan-status" className="text-sm text-ink-soft shrink-0">
            {isProcessing
              ? (status === "scanning" ? "Scanning…" : "Analyzing…")
              : status === "failed" ? "Scan failed" : "Up to date"}
          </div>
        </div>

        {isProcessing && (
          <div className="mt-3">
            <ProgressBar
              value={inFlight ? inFlight.processed : 0}
              total={inFlight ? Math.max(inFlight.total, 1) : 1}
              label="Ingest progress"
            />
          </div>
        )}

        <StatBar aria-label="Experiment stats" stats={stats} className="my-5" />

        {/* Corpus | Configuration tab bar */}
        <nav data-testid="experiment-tab-bar" aria-label="Experiment views" className="flex gap-1 border-b border-hair">
          {TABS.map((t) => {
            const to = `/experiments/${expId}/${t.id}`;
            const isActive = pathname.startsWith(to);
            return (
              <Link
                key={t.id}
                to={to}
                data-testid={`exp-tab-${t.id}`}
                aria-current={isActive ? "page" : undefined}
                className={
                  "px-3 py-2 text-sm font-semibold no-underline -mb-px border-b-2 " +
                  (isActive ? "text-ink border-accent" : "text-ink-soft border-transparent")
                }
              >
                {t.label}
              </Link>
            );
          })}
        </nav>

        <div className="flex-1 min-h-0 pt-5">
          <Outlet />
        </div>
      </PageFrame>
    </div>
  );
}
```

> `ProgressBar`'s API (confirmed against `src/print/ui/ProgressBar.tsx`): props are `value: number`, `total: number` (REQUIRED — NOT `max`), and optional `label`. The `pct` is `value/total`. All colours are `@theme` tokens; `border-b-2 border-accent` is a 2px tab underline (not a side-stripe — the guard only flags `border-l/r` > 1px), so the design guard passes. `Input variant="title"` is confirmed in `Input.tsx` as a valid variant (renders a styled headline-scale input). The wrapper div carries `data-testid`; `onBlur`/`onKeyDown` flow through `...rest` onto the inner `<input>` — the same `...rest` spread that `DirectoryPickerField` relies on. The `updateExperiment` call on commit threads `username` from Zustand (same pattern as other PATCH callers).

- [ ] **Step 4: Run it green** — `npm test -- test/ExperimentShell.test.tsx` → PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/shell/ExperimentShell.tsx packages/HimalayaUI/frontend/test/ExperimentShell.test.tsx
git commit -m "feat(shell): ExperimentShell layout (own chrome + header + StatBar + Corpus|Config tabs)"
```

---

## Task 11: `ExperimentsHomePage` (`src/print/pages/ExperimentsHomePage.tsx`)

`/experiments`: the gallery of experiment cards (status dot · serif name · `date range · directory` · sample/exposure counts · "N need grouping review") + a `+ New experiment` CTA, or an `EmptyState` when there are none. Selecting a card persists `activeExperimentId` to Zustand and routes to that experiment's corpus.

**Files:**
- Create: `src/print/pages/ExperimentsHomePage.tsx`
- Test: `test/ExperimentsHomePage.test.tsx` (CREATE)

- [ ] **Step 1: Write the failing test**

```tsx
// test/ExperimentsHomePage.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ExperimentsHomePage } from "../src/print/pages/ExperimentsHomePage";
import * as api from "../src/api";

const navigate = vi.fn();
vi.mock("react-router-dom", async (orig) => {
  const m = await orig<typeof import("react-router-dom")>();
  return { ...m, useNavigate: () => navigate };
});

function renderPage() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter><ExperimentsHomePage /></MemoryRouter>
    </QueryClientProvider>,
  );
}

const EXP = (over: Partial<api.Experiment>): api.Experiment => ({
  id: 7, name: "SSRL · 1p7m", description: null, path: "/d", data_dir: "/d/data", analysis_dir: "/d/an",
  manifest_path: null, created_at: "2026-04-12", q_units: "A^-1", beam_center_x: 1,
  beam_center_y: 1, pixel_size_um: 172, energy_kev: 9, flight_path_m: 1.8,
  energy_kev_source: "prp", flight_path_m_source: "setup", beam_center_x_source: "setup",
  beam_center_y_source: "setup", pixel_size_um_source: "prp", q_units_source: "prp",
  last_scanned_at: null, scan_signature: null, ingest_status: "complete",
  image_pattern: null, metadata_pattern: null, integration_pattern: null, ...over,
});

describe("ExperimentsHomePage (Phase E1)", () => {
  beforeEach(() => { navigate.mockClear(); });
  afterEach(() => vi.restoreAllMocks());

  it("lists experiment cards", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([EXP({}), EXP({ id: 8, name: "JC plate" })]);
    renderPage();
    expect(await screen.findByText("SSRL · 1p7m")).toBeInTheDocument();
    expect(screen.getByText("JC plate")).toBeInTheDocument();
  });

  it("routes to the corpus tab on card click", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([EXP({})]);
    renderPage();
    fireEvent.click(await screen.findByTestId("experiment-card-7"));
    expect(navigate).toHaveBeenCalledWith("/experiments/7/corpus");
  });

  it("routes to /experiments/new on the New CTA", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([EXP({})]);
    renderPage();
    fireEvent.click(await screen.findByTestId("new-experiment-cta"));
    expect(navigate).toHaveBeenCalledWith("/experiments/new");
  });

  it("shows an empty state when there are no experiments", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
    renderPage();
    expect(await screen.findByText(/no experiments yet/i)).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run it red** — `npm test -- test/ExperimentsHomePage.test.tsx` → FAIL (no page).

- [ ] **Step 3: Implement** — create `src/print/pages/ExperimentsHomePage.tsx`:

```tsx
import { useNavigate } from "react-router-dom";
import { useExperiments } from "../../queries";
import { useAppState } from "../../state";
import { PageFrame } from "../components/PageFrame";
import { Card } from "../ui/Card";
import { Button } from "../ui/Button";
import { Kicker } from "../ui/Kicker";
import { Dot } from "../ui/Dot";
import { EmptyState } from "../ui/EmptyState";
import type { DotTone } from "../ui/Dot";
import type { Experiment } from "../../api";

// Dot tones are accent | success | muted | neutral (Dot.tsx — there is NO
// `ok`/`danger` tone). failed → muted (a quiet hollow ring, NOT a red alarm —
// the row's "Scan failed" text carries the severity), scanning/analyzing →
// accent, complete → success.
function statusTone(s: Experiment["ingest_status"]): DotTone {
  if (s === "failed") return "muted";
  if (s === "scanning" || s === "analyzing") return "accent";
  if (s === "complete") return "success";
  return "neutral";
}

/**
 * ExperimentsHomePage — /experiments gallery (spec §7/§8.7). Cards carry a
 * status dot, serif name, `date range · directory` meta, counts, and a
 * "N need grouping review" hint; selecting one persists activeExperimentId and
 * routes to its corpus. Empty → EmptyState + New CTA.
 */
export function ExperimentsHomePage(): JSX.Element {
  const navigate = useNavigate();
  const setActiveExperiment = useAppState((s) => s.setActiveExperiment);
  const experiments = useExperiments();

  const open = (id: number): void => {
    setActiveExperiment(id);
    navigate(`/experiments/${id}/corpus`);
  };

  const list = experiments.data ?? [];

  return (
    <PageFrame width="home" className="px-6 py-8">
      <div className="flex items-start justify-between gap-6 mb-6">
        <div>
          <Kicker>Experiments</Kicker>
          <h1 className="text-display text-ink">Your beamtimes</h1>
        </div>
        <Button
          variant="accent"
          data-testid="new-experiment-cta"
          onClick={() => navigate("/experiments/new")}
        >
          + New experiment
        </Button>
      </div>

      {experiments.isSuccess && list.length === 0 ? (
        <EmptyState
          title="No experiments yet"
          body="Point Himalaya at a directory of exposures and it scans, groups them into samples, and derives the geometry — no manifest."
          action={
            <Button variant="accent" onClick={() => navigate("/experiments/new")}>
              + New experiment
            </Button>
          }
        />
      ) : (
        <ul className="flex flex-col gap-3">
          {list.map((e) => (
            <Card
              as="li"
              key={e.id}
              interactive
              padding="md"
              data-testid={`experiment-card-${e.id}`}
              onClick={() => open(e.id)}
            >
              <div className="flex items-center gap-3">
                <Dot tone={statusTone(e.ingest_status)} />
                <span className="text-headline text-ink">{e.name ?? `Experiment ${e.id}`}</span>
                <span className="ml-auto font-mono text-sm text-ink-soft">{e.data_dir}</span>
              </div>
            </Card>
          ))}
        </ul>
      )}
    </PageFrame>
  );
}
```

> Confirmed against the live primitives: `Button` has `variant="accent"` (`Button.tsx` `ButtonVariant` = solid|accent|success|ghost|danger|outline|ghostInverse); `Dot`'s `tone` ∈ `accent|success|muted|neutral` (NO `ok`/`danger` — mapped above). `Button` accepts `data-testid` via `...props` (`ButtonHTMLAttributes`) and `Card` accepts `onClick`/`data-testid` via `...rest` (`ElementProps[T]`); `Card`'s `as="li"`/`interactive`/`padding="md"` are real own-props. The counts/"N need grouping review" hint come from the loads roll-up (`useLoads`, samples whose `flag` is non-null) once the page wires it — E1 shows name + dir + status.

- [ ] **Step 4: Run it green** — `npm test -- test/ExperimentsHomePage.test.tsx` → PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/pages/ExperimentsHomePage.tsx packages/HimalayaUI/frontend/test/ExperimentsHomePage.test.tsx
git commit -m "feat(pages): ExperimentsHomePage gallery + empty state + New CTA"
```

---

## Task 12: `NewExperimentPage` (`src/print/pages/NewExperimentPage.tsx`)

`/experiments/new`: wraps `DirectoryPickerField`, fetches suggestions (`suggestPaths`) + validation (`validatePath`) as the path changes, shows an optional name field + a "File patterns" advanced disclosure, and on submit calls `createExperiment` → navigates to the new experiment's corpus (where the live-ingest unfolds). For E1 the suggestion/validation fetching can be debounced-on-change via local effects; the test mocks the API.

**Files:**
- Create: `src/print/pages/NewExperimentPage.tsx`
- Test: `test/NewExperimentPage.test.tsx` (CREATE)

- [ ] **Step 1: Write the failing test**

```tsx
// test/NewExperimentPage.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { NewExperimentPage } from "../src/print/pages/NewExperimentPage";
import * as api from "../src/api";

const navigate = vi.fn();
vi.mock("react-router-dom", async (orig) => {
  const m = await orig<typeof import("react-router-dom")>();
  return { ...m, useNavigate: () => navigate };
});

function renderPage() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter><NewExperimentPage /></MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("NewExperimentPage (Phase E1)", () => {
  beforeEach(() => { navigate.mockClear(); });
  afterEach(() => vi.restoreAllMocks());

  it("renders the directory picker", () => {
    renderPage();
    expect(screen.getByTestId("dirpicker-input")).toBeInTheDocument();
  });

  it("creates the experiment and routes to its corpus on submit", async () => {
    vi.spyOn(api, "validatePath").mockResolvedValue({ ok: true, matched: 682, scanned: 700, message: null });
    vi.spyOn(api, "createExperiment").mockResolvedValue({ id: 9 } as api.Experiment);
    renderPage();
    fireEvent.change(screen.getByTestId("dirpicker-input").querySelector("input")!, {
      target: { value: "/Volumes/data/ssrl/2026_04/1p7m" },
    });
    await waitFor(() => expect(screen.getByTestId("create-submit")).not.toBeDisabled());
    fireEvent.click(screen.getByTestId("create-submit"));
    await waitFor(() => expect(api.createExperiment).toHaveBeenCalled());
    expect(navigate).toHaveBeenCalledWith("/experiments/9/corpus");
  });

  it("keeps submit disabled until validation is ok", () => {
    renderPage();
    expect(screen.getByTestId("create-submit")).toBeDisabled();
  });
});
```

- [ ] **Step 2: Run it red** — `npm test -- test/NewExperimentPage.test.tsx` → FAIL (no page).

- [ ] **Step 3: Implement** — create `src/print/pages/NewExperimentPage.tsx`:

```tsx
import { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import { useAppState } from "../../state";
import * as api from "../../api";
import { authOpts } from "../../lib/authOpts";
import { getClientId } from "../../lib/clientId";
import { PageFrame } from "../components/PageFrame";
import { DirectoryPickerField } from "../components/DirectoryPickerField";
import { Button } from "../ui/Button";
import { Kicker } from "../ui/Kicker";
import { Input } from "../ui/Input";
import type { ValidatePathResponse } from "../../api";

const CLIENT_ID = getClientId();

/**
 * NewExperimentPage — /experiments/new (spec §8.7). Directory picker +
 * suggestions + validation + optional name + a File-patterns advanced
 * disclosure; submit creates the experiment and routes to its corpus, where the
 * live-ingest unfold runs. No file browser.
 */
export function NewExperimentPage(): JSX.Element {
  const navigate = useNavigate();
  const username = useAppState((s) => s.username);
  const setActiveExperiment = useAppState((s) => s.setActiveExperiment);

  const [path, setPath] = useState("");
  const [name, setName] = useState("");
  const [suggestions, setSuggestions] = useState<string[]>([]);
  const [validation, setValidation] = useState<ValidatePathResponse | null>(null);
  const [submitting, setSubmitting] = useState(false);

  // Debounced suggestion + validation fetch on path change.
  useEffect(() => {
    if (path.trim() === "") {
      setSuggestions([]);
      setValidation(null);
      return;
    }
    let live = true;
    const t = setTimeout(() => {
      void api.suggestPaths(path).then((r) => { if (live) setSuggestions(r.suggestions); }).catch(() => {});
      void api.validatePath(path).then((r) => { if (live) setValidation(r); }).catch(() => {});
    }, 200);
    return () => { live = false; clearTimeout(t); };
  }, [path]);

  const canSubmit = validation?.ok === true && !submitting;

  const submit = async (): Promise<void> => {
    if (!canSubmit) return;
    setSubmitting(true);
    try {
      const exp = await api.createExperiment(
        { path, ...(name.trim() ? { name: name.trim() } : {}) },
        // authOpts(username, clientId, clientOpId?) — clientId is POSITIONAL
        // (lib/authOpts.ts). Thread CLIENT_ID; no per-op idempotency key here
        // (create is a one-shot, not a queue mutation).
        authOpts(username, CLIENT_ID),
      );
      setActiveExperiment(exp.id);
      navigate(`/experiments/${exp.id}/corpus`);
    } finally {
      setSubmitting(false);
    }
  };

  return (
    <PageFrame width="home" className="px-6 py-8">
      <button
        type="button"
        onClick={() => navigate("/experiments")}
        className="text-sm text-ink-soft hover:text-ink mb-4"
      >
        ← Experiments
      </button>
      <Kicker>New experiment</Kicker>
      <h1 className="text-display text-ink mb-1">Point at a beamtime</h1>
      <p className="text-base text-ink-soft mb-6">
        Choose the directory of exposures. Himalaya reads the PRP and setup files,
        groups the frames into samples, and derives the geometry.
      </p>

      <div className="flex flex-col gap-5 max-w-[680px]">
        <div>
          <label htmlFor="dirpicker" className="block text-xs font-bold uppercase tracking-wide text-ink-faint mb-1.5">
            Data directory
          </label>
          <DirectoryPickerField
            value={path}
            onChange={setPath}
            suggestions={suggestions}
            validation={validation}
          />
        </div>

        <div>
          <label htmlFor="exp-name" className="block text-xs font-bold uppercase tracking-wide text-ink-faint mb-1.5">
            Experiment name <span className="font-normal normal-case text-ink-soft">(optional)</span>
          </label>
          <Input
            testId="exp-name"
            value={name}
            onValueChange={setName}
            placeholder="Defaults to the directory name"
            aria-label="Experiment name"
          />
        </div>

        <div className="flex items-center gap-3">
          <Button variant="ghost" onClick={() => navigate("/experiments")}>Cancel</Button>
          <Button
            variant="accent"
            data-testid="create-submit"
            disabled={!canSubmit}
            onClick={() => { void submit(); }}
          >
            Scan and create
          </Button>
        </div>
      </div>
    </PageFrame>
  );
}
```

> `authOpts` lives in `src/lib/authOpts.ts` (NOT `queries.ts` — confirmed; `queries.ts` itself imports it from there). Signature is `authOpts(username, clientId, clientOpId?)` with clientId POSITIONAL — pass `CLIENT_ID` from `getClientId()`. `Button`'s `ghost`/`accent` variants confirmed against `Button.tsx`. The File-patterns advanced disclosure (image/metadata/integration patterns) is shown in the prototype; for E1 it can be added as a collapsed `<details>` slotting three `Input`s threaded into `createExperiment`'s `patterns` — include it only if straightforward, else defer the patterns plumbing to E2 (the create body already types `patterns?`).

- [ ] **Step 4: Run it green** — `npm test -- test/NewExperimentPage.test.tsx` → PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/pages/NewExperimentPage.tsx packages/HimalayaUI/frontend/test/NewExperimentPage.test.tsx
git commit -m "feat(pages): NewExperimentPage directory picker + create flow"
```

---

## Task 13: `ExperimentCorpusPage` (`src/print/pages/ExperimentCorpusPage.tsx`)

The Corpus tab body: reuses the shipped contact-sheet `SheetTable` scoped to the experiment, with a sticky **"N samples need grouping review →"** banner above it that routes to the grouping surface (`/experiments/:id/grouping`). The grouping surface internals are **E2's `GroupingReviewPage`** — this page only renders the banner + the table and links to it. The live-ingest unfold (E2's `LiveIngestUnfold`) replaces the table while `ingestInFlight[id]` is active; for E1, render a labelled placeholder where `LiveIngestUnfold` will slot.

**Files:**
- Create: `src/print/pages/ExperimentCorpusPage.tsx`
- Test: `test/ExperimentCorpusPage.test.tsx` (CREATE)

- [ ] **Step 1: Write the failing test**

```tsx
// test/ExperimentCorpusPage.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ExperimentCorpusPage } from "../src/print/pages/ExperimentCorpusPage";
import { useAppState } from "../src/state";
import * as api from "../src/api";
import type { Load } from "../src/api";

// Minimal load stub with a flagged merge discrepancy.
const LOAD_WITH_FLAG: Load = {
  load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null,
  frame_count: 4, note: null,
  samples: [
    { sample_id: 9, name: "HA85", slot_index: 1, grouping_source: "computed",
      name_source: "computed", merged_into_id: null,
      flag: { kind: "merge", merge_with_sample_id: 10, merge_with_label: "HA86" },
      exposures: [] },
    { sample_id: 10, name: "HA86", slot_index: 2, grouping_source: "computed",
      name_source: "computed", merged_into_id: null, flag: null, exposures: [] },
    { sample_id: 11, name: "HA87", slot_index: 3, grouping_source: "computed",
      name_source: "computed", merged_into_id: null,
      flag: { kind: "split", split_at_index: 2, jump_from: 12.4, jump_to: 48.1 },
      exposures: [] },
    { sample_id: 12, name: "HA88", slot_index: 4, grouping_source: "computed",
      name_source: "computed", merged_into_id: null, flag: null, exposures: [] },
  ],
};

function renderAt(loads: Load[], processing = false) {
  if (processing) useAppState.setState({ ingestInFlight: { 7: { processed: 10, total: 100, status: "scanning" } } });
  else useAppState.setState({ ingestInFlight: null });
  vi.spyOn(api, "listLoads").mockResolvedValue(loads);
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/experiments/7/corpus"]}>
        <Routes>
          <Route path="/experiments/:id/corpus" element={<ExperimentCorpusPage />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("ExperimentCorpusPage (Phase E1)", () => {
  beforeEach(() => vi.restoreAllMocks());
  afterEach(() => { useAppState.setState({ ingestInFlight: null }); vi.restoreAllMocks(); });

  it("shows the grouping-review banner with the count + link when reviews are pending", async () => {
    renderAt([LOAD_WITH_FLAG]);
    // 2 flagged samples: LOAD_WITH_FLAG.samples[0] (merge) + [2] (split)
    const banner = await screen.findByTestId("grouping-review-banner");
    expect(banner).toHaveTextContent("2");
    expect(screen.getByTestId("grouping-review-link")).toHaveAttribute("href", "/experiments/7/grouping");
  });

  it("hides the banner when nothing needs review", async () => {
    // Load with no flagged samples.
    const cleanLoad: Load = { ...LOAD_WITH_FLAG,
      samples: LOAD_WITH_FLAG.samples.map((s) => ({ ...s, flag: null })) };
    renderAt([cleanLoad]);
    await waitFor(() => expect(api.listLoads).toHaveBeenCalled());
    expect(screen.queryByTestId("grouping-review-banner")).toBeNull();
  });

  it("renders the live-ingest placeholder while processing", async () => {
    renderAt([], true);
    expect(screen.getByTestId("live-ingest-slot")).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run it red** — `npm test -- test/ExperimentCorpusPage.test.tsx` → FAIL (no page).

- [ ] **Step 3: Implement** — create `src/print/pages/ExperimentCorpusPage.tsx`:

```tsx
import { Link, useParams } from "react-router-dom";
import { useAppState } from "../../state";
import { useLoads } from "../../queries";
import { Badge } from "../ui/Badge";

/**
 * ExperimentCorpusPage — the Corpus tab body. Reuses the shipped SheetTable
 * contact sheet scoped to the experiment (E2 wires the table), with a sticky
 * grouping-review banner above it linking to E2's GroupingReviewPage
 * (/experiments/:id/grouping). The live-ingest unfold (E2 LiveIngestUnfold)
 * replaces the table while ingestInFlight[id] is active.
 */
export function ExperimentCorpusPage(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  const expId = id ? Number(id) : 0;
  const inFlight = useAppState((s) => s.ingestInFlight?.[expId]);
  const loads = useLoads(expId);

  // Review count: LoadSamples across all loads whose `flag` is non-null
  // (a flagged merge/split discrepancy). Derived from useLoads — tests mock
  // api.listLoads via vi.spyOn, same pattern as Tasks 11–14.
  const reviewCount = (loads.data ?? []).reduce(
    (n, l) => n + l.samples.filter((s) => s.flag !== null).length,
    0,
  );
  const processing = inFlight?.status === "scanning" || inFlight?.status === "analyzing";

  return (
    <div className="flex flex-col gap-4">
      {reviewCount > 0 && (
        <div
          data-testid="grouping-review-banner"
          className="sticky top-0 z-10 flex items-center gap-3 rounded-sm border border-hair-strong bg-paper-sunk px-4 py-2.5"
        >
          <span className="text-sm text-ink">
            {reviewCount} {reviewCount === 1 ? "sample needs" : "samples need"} grouping review
            <Badge>{reviewCount}</Badge>
          </span>
          <Link
            to={`/experiments/${expId}/grouping`}
            data-testid="grouping-review-link"
            className="ml-auto text-sm font-semibold text-accent hover:underline"
          >
            Review grouping →
          </Link>
        </div>
      )}

      {processing ? (
        // E2 LiveIngestUnfold slots here (ProgressBar + StatBar counters +
        // skeleton rows driven by ingestInFlight). E1 renders the labelled slot.
        <div data-testid="live-ingest-slot" className="text-sm text-ink-soft">
          Processing exposures…
        </div>
      ) : (
        // E2 mounts the scoped SheetTable here. E1 renders the labelled slot so
        // the page is assemblable without the corpus query wiring.
        <div data-testid="corpus-sheet-slot" className="text-sm text-ink-soft">
          {/* SheetTable scoped to experiment {expId} — wired in E2 */}
        </div>
      )}
    </div>
  );
}
```

> `NoticePill` has only `new`/`draft` tones (`NoticePillTone` in `NoticePill.tsx`) — there is NO `warning` tone, so it cannot carry the review count. Use `Badge` (the inline mono count primitive) inside the banner span instead; the banner itself is a token-styled div (`bg-paper-sunk` + `border-hair-strong`). The scoped `SheetTable` reuse + the corpus-samples query plumbing belong to E2 (it owns the grouping-review data plane); E1 lands the banner (count derived from `useLoads`) + the slots. All colours are tokens; `rounded-sm` is the 5px step. The test uses `vi.spyOn(api, "listLoads")` — the same QueryClient wrapper pattern used by Tasks 11–14 — rather than a test-only prop, keeping `ExperimentCorpusPage`'s public API free of production seams.

- [ ] **Step 4: Run it green** — `npm test -- test/ExperimentCorpusPage.test.tsx` → PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/pages/ExperimentCorpusPage.tsx packages/HimalayaUI/frontend/test/ExperimentCorpusPage.test.tsx
git commit -m "feat(pages): ExperimentCorpusPage banner + live-ingest/sheet slots"
```

---

## Task 14: `ConfigurationPage` shell (`src/print/pages/ConfigurationPage.tsx`)

The Configuration tab body shell: editable description region, a **Geometry** ledger region, and a **Sources** region. The internals — `GeometryLedger`, `AcquisitionTimeline`, `SourcesCard` — are **E2 components**; this page lays out the three regions and slots them as labelled placeholders, reading the experiment via `useExperiment` (called directly).

> **E2 ownership of SourcesCard:** E2 builds `SourcesCard` with the following field layout (DECISION):
> - `image_pattern`, `metadata_pattern`, `integration_pattern` — **EDITABLE** edit-in-place rows wired to `updateExperiment` (patch triggers a rescan server-side via `scan_signature` invalidation).
> - `data_dir`, `analysis_dir` — **READ-ONLY** display rows. The directory is fixed at experiment creation. **Do NOT add a DirectoryPickerField, path autocomplete, or validate-path backend call here** — spec §13's directory-edit / FS-sandboxing open question is moot for v1.
> - Acquisition timeline: read-only / derived from the loads roll-up.

**Files:**
- Create: `src/print/pages/ConfigurationPage.tsx`
- Test: `test/ConfigurationPage.test.tsx` (CREATE)

- [ ] **Step 1: Write the failing test**

```tsx
// test/ConfigurationPage.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ConfigurationPage } from "../src/print/pages/ConfigurationPage";
import * as api from "../src/api";

const EXP: api.Experiment = {
  id: 7, name: "SSRL", description: null, path: "/d", data_dir: "/d/data", analysis_dir: "/d/an",
  manifest_path: null, created_at: "2026-04-12", q_units: "A^-1", beam_center_x: 421,
  beam_center_y: 836, pixel_size_um: 172, energy_kev: 9, flight_path_m: 1.8095,
  energy_kev_source: "prp", flight_path_m_source: "setup", beam_center_x_source: "setup",
  beam_center_y_source: "setup", pixel_size_um_source: "prp", q_units_source: "prp",
  last_scanned_at: null, scan_signature: null, ingest_status: "complete",
  image_pattern: null, metadata_pattern: null, integration_pattern: null,
};

function renderPage() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/experiments/7/config"]}>
        <Routes>
          <Route path="/experiments/:id/config" element={<ConfigurationPage />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("ConfigurationPage (Phase E1 shell)", () => {
  beforeEach(() => { vi.spyOn(api, "getExperiment").mockResolvedValue(EXP); });
  afterEach(() => vi.restoreAllMocks());

  it("renders the Geometry and Sources regions", async () => {
    renderPage();
    expect(await screen.findByTestId("config-geometry-region")).toBeInTheDocument();
    expect(screen.getByTestId("config-sources-region")).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run it red** — `npm test -- test/ConfigurationPage.test.tsx` → FAIL (no page).

- [ ] **Step 3: Implement** — create `src/print/pages/ConfigurationPage.tsx`:

```tsx
import { useParams } from "react-router-dom";
import { useExperiment } from "../../queries";
import { Card } from "../ui/Card";
import { Kicker } from "../ui/Kicker";

/**
 * ConfigurationPage — the Configuration tab body shell (spec §8.1). Lays out
 * three regions: editable description, the Geometry ledger, and Sources. The
 * region INTERNALS are E2 components — GeometryLedger (per-field value +
 * provenance chip + Override + discrepancy banner), AcquisitionTimeline, and
 * SourcesCard. E1 lays out + slots them as labelled placeholders.
 */
export function ConfigurationPage(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  const expId = id ? Number(id) : 0;
  const exp = useExperiment(expId);

  return (
    <div className="flex flex-col gap-6">
      <Card padding="md" data-testid="config-geometry-region">
        <Kicker>Geometry</Kicker>
        {/* E2 GeometryLedger: per-field value + prp/setup/user chip + Override
            + the multi-setup discrepancy banner. E1 shows the derived values. */}
        <dl className="mt-3 grid grid-cols-[auto_1fr] gap-x-6 gap-y-1 text-sm">
          <dt className="text-ink-soft">Energy</dt>
          <dd className="font-mono text-ink">{exp.data?.energy_kev ?? "—"} keV</dd>
          <dt className="text-ink-soft">Flight path</dt>
          <dd className="font-mono text-ink">{exp.data?.flight_path_m ?? "—"} m</dd>
          <dt className="text-ink-soft">Beam center</dt>
          <dd className="font-mono text-ink">
            {exp.data ? `${exp.data.beam_center_x ?? "—"}, ${exp.data.beam_center_y ?? "—"}` : "—"}
          </dd>
        </dl>
      </Card>

      <Card padding="md" data-testid="config-sources-region">
        <Kicker>Sources</Kicker>
        {/* E2 SourcesCard: editable pattern rows + read-only directory rows.
            DECISION (E1→E2 contract):
            - image_pattern / metadata_pattern / integration_pattern: EDITABLE
              edit-in-place rows wired to updateExperiment (rescan-on-change).
            - data_dir / analysis_dir: READ-ONLY display rows (directory is
              fixed at creation; no DirectoryPickerField / path-validation here).
            E1 shows the static read-only dirs as a placeholder. */}
        <dl className="mt-3 grid grid-cols-[auto_1fr] gap-x-6 gap-y-1 text-sm">
          <dt className="text-ink-soft">Data</dt>
          <dd className="font-mono text-ink truncate">{exp.data?.data_dir ?? "—"}</dd>
          <dt className="text-ink-soft">Analysis</dt>
          <dd className="font-mono text-ink truncate">{exp.data?.analysis_dir ?? "—"}</dd>
        </dl>
      </Card>
    </div>
  );
}
```

- [ ] **Step 4: Run it green** — `npm test -- test/ConfigurationPage.test.tsx` → PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/pages/ConfigurationPage.tsx packages/HimalayaUI/frontend/test/ConfigurationPage.test.tsx
git commit -m "feat(pages): ConfigurationPage shell (Geometry + Sources regions)"
```

---

## Task 15: `AppRoutes` — the `/experiments/*` routing tree + `/` redirect

Add `/experiments` (home), `/experiments/new`, and an `/experiments/:id` layout route (`ExperimentShell`) with children `corpus` / `grouping` / `config`, **outside** `CorpusShell`. Redirect `/` → `/experiments` (was `/samples`). The `grouping` child route element is E2's `GroupingReviewPage`; for E1, route it to a labelled placeholder element and cross-reference E2.

**Files:**
- Modify: `src/print/shell/AppRoutes.tsx` (imports + the `<Routes>` table)
- Test: `test/AppRoutes.ingestion.test.tsx` (CREATE)

- [ ] **Step 1: Write the failing test**

```tsx
// test/AppRoutes.ingestion.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { AppRoutes } from "../src/print/shell/AppRoutes";
import * as api from "../src/api";

function renderAt(path: string) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[path]}><AppRoutes /></MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("AppRoutes ingestion tree (Phase E1)", () => {
  beforeEach(() => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
    vi.spyOn(api, "getExperiment").mockResolvedValue({ id: 7, name: "X", ingest_status: "complete" } as api.Experiment);
  });
  afterEach(() => vi.restoreAllMocks());

  it("/experiments renders the home gallery", async () => {
    renderAt("/experiments");
    expect(await screen.findByText("Your beamtimes")).toBeInTheDocument();
  });

  it("/experiments/new renders the directory picker", () => {
    renderAt("/experiments/new");
    expect(screen.getByTestId("dirpicker-input")).toBeInTheDocument();
  });

  it("/experiments/:id/corpus mounts ExperimentShell chrome + corpus body", async () => {
    renderAt("/experiments/7/corpus");
    expect(await screen.findByTestId("experiment-shell")).toBeInTheDocument();
    expect(screen.getByTestId("experiment-top-nav")).toBeInTheDocument();
    // ExperimentShell is OUTSIDE CorpusShell → no corpus topbar.
    expect(screen.queryByTestId("corpus-topbar")).toBeNull();
  });

  it("/experiments/:id/grouping mounts the grouping route element", async () => {
    renderAt("/experiments/7/grouping");
    expect(await screen.findByTestId("grouping-review-placeholder")).toBeInTheDocument();
  });

  it("/ redirects to /experiments", async () => {
    renderAt("/");
    expect(await screen.findByText("Your beamtimes")).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run it red** — `npm test -- test/AppRoutes.ingestion.test.tsx` → FAIL (routes not registered; `/` still goes to `/samples`).

- [ ] **Step 3: Implement**

Add imports at the top of `src/print/shell/AppRoutes.tsx`:

```tsx
import { ExperimentShell } from "./ExperimentShell";
import { ExperimentsHomePage } from "../pages/ExperimentsHomePage";
import { NewExperimentPage } from "../pages/NewExperimentPage";
import { ExperimentCorpusPage } from "../pages/ExperimentCorpusPage";
import { ConfigurationPage } from "../pages/ConfigurationPage";
```

Add the `/experiments/*` routes **before** the closing `</Routes>`, OUTSIDE the `CorpusShell` layout route (place them alongside the other outside-shell redirects, after line 127):

```tsx
      {/* Ingestion redesign (spec §7/§9.6): the experiments tree sits OUTSIDE
          CorpusShell so ExperimentShell's own chrome never stacks on
          CorpusTopbar. */}
      <Route path="/experiments" element={<ExperimentsHomePage />} />
      <Route path="/experiments/new" element={<NewExperimentPage />} />
      <Route path="/experiments/:id" element={<ExperimentShell />}>
        <Route index element={<Navigate to="corpus" replace />} />
        <Route path="corpus" element={<ExperimentCorpusPage />} />
        <Route path="config" element={<ConfigurationPage />} />
        {/* E2 GroupingReviewPage mounts here. E1 routes to a placeholder. */}
        <Route
          path="grouping"
          element={<div data-testid="grouping-review-placeholder">Grouping review (E2)</div>}
        />
      </Route>
```

Change the `/` redirect (line 116) from `/samples` to `/experiments`:

```tsx
      <Route path="/" element={<Navigate to="/experiments" replace />} />
```

> Leave `/samples`, `/series`, and the loupe/focus routes under `CorpusShell` unchanged — the corpus pages still exist; only the entry default moves. `/series` stays as the existing folio surface (the Series axis ships as-is per spec §7). Update any e2e spec that asserted `/`→`/samples` separately (spec §9.6 — note as a follow-up; e2e is Task 22's concern only if a spec breaks the build, which it does not since `npm run build` doesn't run Playwright).

- [ ] **Step 4: Run it green** — `npm test -- test/AppRoutes.ingestion.test.tsx` → PASS. Then run the existing `npm test -- test/AppRoutes.test.tsx` and fix any assertion that depended on `/`→`/samples` (update it to `/experiments`).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/shell/AppRoutes.tsx packages/HimalayaUI/frontend/test/AppRoutes.ingestion.test.tsx packages/HimalayaUI/frontend/test/AppRoutes.test.tsx
git commit -m "feat(shell): /experiments routing tree outside CorpusShell; / -> /experiments"
```

---

## Task 16: `CorpusTopbar` — add the Experiments stage

Per spec §9.6: add "Experiments" to `CorpusTopbar`'s `STAGES` array so the non-experiment corpus surfaces (`/samples`, `/series`) still expose a link back to the experiments home. `pathname.startsWith` active detection already handles it.

**Files:**
- Modify: `src/print/shell/CorpusTopbar.tsx` (`Stage` type ~16-21; `STAGES` ~27-30)
- Test: `test/CorpusTopbar.ingestion.test.tsx` (CREATE)

- [ ] **Step 1: Write the failing test**

```tsx
// test/CorpusTopbar.ingestion.test.tsx
import { describe, it, expect, vi, afterEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { CorpusTopbar } from "../src/print/shell/CorpusTopbar";
import * as api from "../src/api";

function at(path: string) {
  vi.spyOn(api, "listExperiments").mockResolvedValue([]);
  vi.spyOn(api, "listCorpusSamples").mockResolvedValue([]);
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[path]}><CorpusTopbar /></MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("CorpusTopbar Experiments stage (Phase E1)", () => {
  afterEach(() => vi.restoreAllMocks());

  it("includes an Experiments stage tab linking to /experiments", () => {
    at("/samples");
    expect(screen.getByTestId("stage-tab-experiments")).toHaveAttribute("href", "/experiments");
  });
});
```

- [ ] **Step 2: Run it red** — `npm test -- test/CorpusTopbar.ingestion.test.tsx` → FAIL (no `experiments` stage).

- [ ] **Step 3: Implement**

Widen the `Stage` id union (`src/print/shell/CorpusTopbar.tsx:17`):

```tsx
interface Stage {
  id: "experiments" | "samples" | "series";
  label: string;
  to: string;
}
```

Prepend the Experiments stage to `STAGES` (line 27-30):

```tsx
const STAGES: readonly Stage[] = [
  { id: "experiments", label: "Experiments", to: "/experiments" },
  { id: "samples", label: "Samples", to: "/samples" },
  { id: "series", label: "Series", to: "/series" },
];
```

> No other change — the `stageTabs` map already uses `pathname.startsWith(s.to)` and `stage-tab-${s.id}` testids, so the new tab gets active detection + testid for free.

- [ ] **Step 4: Run it green** — `npm test -- test/CorpusTopbar.ingestion.test.tsx` → PASS. Run `npm test -- test/CorpusTopbar` (any existing topbar specs) and fix count-based assertions if a spec hard-coded "two stage tabs".

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/shell/CorpusTopbar.tsx packages/HimalayaUI/frontend/test/CorpusTopbar.ingestion.test.tsx
git commit -m "feat(shell): add Experiments stage tab to CorpusTopbar"
```

---

## Task 17: ~~`NavModal` — drop the home phase bar~~ — DROPPED (no such element exists)

**DROPPED.** Verified live (2026-06-18): there is NO phase bar / phase-distribution strip in `NavModal.tsx` or anywhere under `src/print/shell/`. `git grep -n "PhaseStrip\|phase-bar\|phaseDistribution\|nav-home-phase" src/print/shell/` returns nothing; `NavModal.tsx` references `phase` only via the `display_name`/`name` search predicates (handled by the Task 1b sweep). The `PhaseStrip` primitive exists in `src/print/ui/PhaseStrip.tsx` but is not consumed by any shell file. Writing a test that asserts a never-rendered element is absent would be vacuous (it passes trivially against today's tree), so this task is removed rather than shipped.

> If the round-3 "drop the home phase bar" note referred to a phase bar that lived on the OLD home/Index surface, that surface was already retired in the greenfield cutover — nothing to remove here. If a future surface (e.g. an ExperimentsHomePage phase summary) is built and then needs the bar dropped, that's a fresh task, not this one. NavModal's `display_name`→`name` reads are converted in Task 1b.

---

## Task 18: `ingest_*` cache arms in `applyRemoteToCache`

Per spec §9.3/§9.6: add the four `ingest_*` cache arms in `applyRemoteToCache` that discriminate on the top-level `kind`, read `payload.experiment_id` (NEVER `remote.entity_id`), and invalidate the experiment's loads/samples caches — placed **before** the `default:` arm so they never reach the poisoning peaks/indices invalidation. The arms do **cache invalidation only**; the `ingestInFlight` store writes live in the separate `App.tsx` listener (Task 20), keeping `applyRemoteToCache` pure.

Also **export a shared `invalidateIngestFrameCache(qc, expId, isComplete)` helper** (~6 lines) so Task 20's App listener can call it instead of duplicating the `invalidateQueries` calls inline. This eliminates the drift risk the self-review flags (open question #4). The helper is pure (no Zustand); store writes stay in the App listener.

> **`SseEvent.kind` is ALREADY `string`** (confirmed `types.ts:148` — `kind: string`, not a literal union). So there is NOTHING to extend on the kind union — the new arms just add `case "ingest_progress":` etc. to the switch. `payload` is typed `unknown` (`types.ts:155`); the arms read it via the existing `payload?.experiment_id as number | undefined` cast idiom (same as every other arm). The ingest frame is **payload-wrapped** — `{ kind, payload: { experiment_id, processed, total } }` riding the "curation" event — so read `payload.experiment_id`, NOT `remote.entity_id`.

**Files:**
- Modify: `src/lib/queue/applyRemoteToCache.ts` (new arms before `default:`)
- Test: `test/applyRemoteToCache.ingest.test.ts` (CREATE)

- [ ] **Step 1: Confirm the shapes** — `SseEvent.kind: string` + `payload: unknown` (`types.ts:146-157`); the `default:` arm fires `peaks(id)`/`indices(id)` invalidations off `remote.entity_id` (lines 334-337). No `types.ts` edit is needed.

- [ ] **Step 2: Write the failing test**

```ts
// test/applyRemoteToCache.ingest.test.ts
import { describe, it, expect, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../src/lib/queue/applyRemoteToCache";
import { queryKeys } from "../src/queries";
import type { SseEvent } from "../src/lib/queue/types";

function ingestFrame(kind: string, experiment_id: number, extra: Record<string, unknown> = {}): SseEvent {
  // broadcast_progress! (events.jl:1158-1186) emits kind + payload at top-level;
  // it does NOT emit a top-level entity_id. The payload sub-object carries
  // experiment_id, processed, total. entity_id is -1 sentinel here (not emitted
  // by the real backend for ingest frames; the arms read payload.experiment_id).
  return {
    kind,
    entity_id: -1,
    entity_type: "experiment",
    payload: { experiment_id, ...extra },
  } as unknown as SseEvent;
}

describe("ingest_* cache arms (Phase E1)", () => {
  it("ingest_progress invalidates the experiment loads + samples, never peaks/indices", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(ingestFrame("ingest_progress", 7, { processed: 100, total: 680 }), qc);
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.loads(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.samples(7)));
    // the default arm's peaks/indices invalidation must NOT fire
    expect(keys).not.toContain(JSON.stringify(queryKeys.peaks(7)));
    expect(keys).not.toContain(JSON.stringify(queryKeys.indices(7)));
  });

  it("ingest_complete invalidates loads + samples + the experiment detail", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(ingestFrame("ingest_complete", 7), qc);
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.experiment(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.loads(7)));
  });

  it("ingest_started / ingest_failed do not throw and do not poison peaks/indices", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(ingestFrame("ingest_started", 7), qc);
    applyRemoteToCache(ingestFrame("ingest_failed", 7), qc);
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).not.toContain(JSON.stringify(queryKeys.peaks(7)));
  });
});
```

- [ ] **Step 3: Run it red** — `npm test -- test/applyRemoteToCache.ingest.test.ts` → FAIL (kinds hit the `default:` arm → peaks/indices invalidations fire on `entity_id`).

- [ ] **Step 4: (no types.ts change needed)** — `SseEvent.kind` is already `string` (confirmed). The new `case` labels below compile against it directly. Skip straight to adding the arms.

- [ ] **Step 5: Add the `invalidateIngestFrameCache` helper and the `ingest_*` arms** in `applyRemoteToCache.ts`.

First, export the shared helper **near the top of the file** (before the `applyRemoteToCache` function), so both the switch arms and the App listener can call it without duplicating the `invalidateQueries` calls:

```ts
/**
 * Shared cache-invalidation helper for ingest-frame side-effects. Exported so
 * `App.tsx`'s SSE listener can call the same invalidations without duplicating
 * them inline (the listener still owns the Zustand store write; this helper is
 * pure cache-only). `isComplete=true` also refetches the experiment detail row
 * (so `ingest_status` transitions from "analyzing" → "complete").
 */
export function invalidateIngestFrameCache(
  qc: QueryClient,
  expId: number,
  isComplete: boolean,
): void {
  qc.invalidateQueries({ queryKey: queryKeys.loads(expId) });
  qc.invalidateQueries({ queryKey: queryKeys.samples(expId) });
  if (isComplete) qc.invalidateQueries({ queryKey: queryKeys.experiment(expId) });
}
```

Then add the `ingest_*` arms **before** the `default:` case (line 334), calling the helper:

```ts
    case "ingest_started":
    case "ingest_progress":
    case "ingest_failed": {
      // Broadcast-only progress (spec §9.3): rides the curation channel, carries
      // a positive experiment_id in the payload (NEVER a sentinel entity_id).
      // Read experiment_id from the payload, not remote.entity_id. Cache effect
      // is invalidation-only (the ingestInFlight store write lives in the
      // separate App.tsx listener — applyRemoteToCache stays pure).
      const expId = payload?.experiment_id as number | undefined;
      if (expId !== undefined) invalidateIngestFrameCache(qc, expId, false);
      break;
    }
    case "ingest_complete": {
      // Authoritative terminal frame (the 64-slot channel may drop progress
      // frames at 680-exposure scale; treat complete as the source of truth).
      const expId = payload?.experiment_id as number | undefined;
      if (expId !== undefined) invalidateIngestFrameCache(qc, expId, true);
      break;
    }
```

- [ ] **Step 6: Run it green** — `npm test -- test/applyRemoteToCache.ingest.test.ts` → PASS.

- [ ] **Step 7: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/applyRemoteToCache.ts packages/HimalayaUI/frontend/test/applyRemoteToCache.ingest.test.ts
git commit -m "feat(queue): ingest_* cache arms before default (payload.experiment_id)"
```

---

## Task 19: SSE receive path for the structural event kinds

Per the canonical contract (spec §9.3/§9.6) the Phase-D structural grouping edits emit these payload-wrapped frames (read top-level `kind` + the `payload` fields below; NEVER `remote.entity_id` for the experiment):

- `sample_renamed { name, experiment_id }` — splice the new `name` into the sample cache (entity_id is the sample id) + invalidate the experiment + corpus sample listings.
- `exposure_moved { sample_id, from_sample_id, experiment_id }` — one exposure left `from_sample_id` for `sample_id` (the destination). Invalidate BOTH samples' exposure lists + the loads roll-up.
- `sample_created { experiment_id }` — a brand-new sample row (e.g. the destination of a split). **Invalidate-only** (loads + both sample listings). **MUST precede `default:`** so its entity_id (a sample id) does not poison `peaks(id)`/`indices(id)`.
- `sample_split { new_sample_id, exposure_ids, experiment_id }` — invalidate-only (loads + sample listings); ids aren't worth a surgical splice (the `series_created` precedent).
- `grouping_flag_dismissed { flag_kind, merge_with_sample_id?, experiment_id }` — a flag was cleared (the "merge" gesture is modeled as dismissing the merge flag, not a `sample_merged` event). Invalidate loads + sample listings so the banner count + the fold re-derive.

> **There is NO `sample_merged` event** (per the contract — drop any such arm; merging is `grouping_flag_dismissed` + the move of exposures). The structural MUTATIONS (mutators + `useQueueMutation` wiring) are **E2**; this task wires only the cache RECEIVE path so foreign-tab edits and own-op echoes reconcile. `SseEvent.kind` is already `string`, so no union edit.

**Files:**
- Modify: `src/lib/queue/applyRemoteToCache.ts` (new arms before `default:`)
- Test: `test/applyRemoteToCache.structural.test.ts` (CREATE)

- [ ] **Step 1: Write the failing test**

```ts
// test/applyRemoteToCache.structural.test.ts
import { describe, it, expect, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../src/lib/queue/applyRemoteToCache";
import { queryKeys } from "../src/queries";
import type { SseEvent } from "../src/lib/queue/types";
import type { Sample } from "../src/api";

function frame(kind: string, entity_id: number, payload: Record<string, unknown>): SseEvent {
  return { kind, entity_id, entity_type: "sample", payload } as unknown as SseEvent;
}

describe("structural grouping SSE receive (Phase E1)", () => {
  it("sample_renamed splices the new name into the sample cache + invalidates listings", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    qc.setQueryData<Sample>(queryKeys.sample(5), {
      id: 5, experiment_id: 7, name: "Old", notes: null, tags: [],
    });
    applyRemoteToCache(frame("sample_renamed", 5, { name: "HA85 (S01P15)", experiment_id: 7 }), qc);
    expect(qc.getQueryData<Sample>(queryKeys.sample(5))?.name).toBe("HA85 (S01P15)");
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.samples(7)));
  });

  it("exposure_moved invalidates both samples' exposures + the loads roll-up", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    // Contract payload: sample_id = destination, from_sample_id = source.
    applyRemoteToCache(
      frame("exposure_moved", 99, { sample_id: 2, from_sample_id: 1, experiment_id: 7 }),
      qc,
    );
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.exposures(1)));
    expect(keys).toContain(JSON.stringify(queryKeys.exposures(2)));
    expect(keys).toContain(JSON.stringify(queryKeys.loads(7)));
  });

  it("sample_created invalidates loads + sample listings, NEVER peaks/indices (precedes default)", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(frame("sample_created", 42, { experiment_id: 7 }), qc);
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.loads(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.samples(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.corpusSamples));
    // entity_id 42 is a sample id — the default arm would poison these:
    expect(keys).not.toContain(JSON.stringify(queryKeys.peaks(42)));
    expect(keys).not.toContain(JSON.stringify(queryKeys.indices(42)));
  });

  it("sample_split / grouping_flag_dismissed refetch loads + sample listings (invalidate-only)", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(frame("sample_split", 2, { new_sample_id: 9, exposure_ids: [1, 2], experiment_id: 7 }), qc);
    applyRemoteToCache(frame("grouping_flag_dismissed", 2, { flag_kind: "merge", merge_with_sample_id: 4, experiment_id: 7 }), qc);
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.loads(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.samples(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.corpusSamples));
  });
});
```

- [ ] **Step 2: Run it red** — `npm test -- test/applyRemoteToCache.structural.test.ts` → FAIL (kinds hit `default:`; `sample_renamed` splice doesn't happen).

- [ ] **Step 3: Implement** — add these arms in `applyRemoteToCache.ts` before `default:` (after the ingest arms from Task 18):

```ts
    case "sample_renamed": {
      // Single-entity, payload-derivable: splice the renamed label into the
      // sample cache + refresh the corpus/experiment listings. entity_id is the
      // sample id; the new label rides payload.name.
      const newName = payload?.name as string | undefined;
      if (newName !== undefined) {
        qc.setQueryData<Sample>(queryKeys.sample(id), (old) =>
          old ? { ...old, name: newName } : old);
      }
      qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
      const expId = payload?.experiment_id as number | undefined;
      if (expId !== undefined) qc.invalidateQueries({ queryKey: queryKeys.samples(expId) });
      break;
    }
    case "exposure_moved": {
      // One exposure left from_sample_id for sample_id (the destination).
      // Invalidate both sides' exposure lists + the loads roll-up (which
      // re-derives the grouping view).
      const dest = payload?.sample_id as number | undefined;
      const from = payload?.from_sample_id as number | undefined;
      if (from !== undefined) qc.invalidateQueries({ queryKey: queryKeys.exposures(from) });
      if (dest !== undefined) qc.invalidateQueries({ queryKey: queryKeys.exposures(dest) });
      const expId = payload?.experiment_id as number | undefined;
      if (expId !== undefined) qc.invalidateQueries({ queryKey: queryKeys.loads(expId) });
      break;
    }
    case "sample_created":
    case "sample_split":
    case "grouping_flag_dismissed": {
      // Orchestrations that create sample rows / change grouping whose new ids
      // aren't worth a surgical splice (the series_created precedent). Refetch
      // the loads roll-up + both sample listings. MUST be before `default:` so
      // entity_id (a sample id) never poisons peaks(id)/indices(id). There is
      // NO sample_merged event — merging is grouping_flag_dismissed + the move.
      const expId = payload?.experiment_id as number | undefined;
      if (expId !== undefined) {
        qc.invalidateQueries({ queryKey: queryKeys.loads(expId) });
        qc.invalidateQueries({ queryKey: queryKeys.samples(expId) });
      }
      qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
      break;
    }
```

- [ ] **Step 4: Run it green** — `npm test -- test/applyRemoteToCache.structural.test.ts` → PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/applyRemoteToCache.ts packages/HimalayaUI/frontend/test/applyRemoteToCache.structural.test.ts
git commit -m "feat(queue): structural grouping SSE receive (rename splice; move/created/split/flag-dismiss invalidate)"
```

---

## Task 20: `App.tsx` — `ingestInFlight` store write + short-circuit before `handleRemoteEvent`

Per spec §9.3/§9.6: the `ingest_*` store writes go in the `App.tsx` SSE listener (Zustand is NOT imported in `applyRemoteToCache`, which stays pure). On an `ingest_*` frame:

1. write/clear the `ingestInFlight` store entry, then
2. call `invalidateIngestFrameCache(qc, expId, isComplete)` (exported from `applyRemoteToCache.ts` in Task 18 — NOT an inline duplicate), then
3. **`return` — do NOT fall through to `handleRemoteEvent`**.

The short-circuit matters: `handleRemoteEvent` runs the full queue reconciliation (own-op deferred match → otherwise rollback-and-rerun every pending mutation). An `ingest_*` frame is broadcast-only — it is no tab's own op — so feeding it to `handleRemoteEvent` would needlessly roll back and re-run the user's in-flight peak/tag edits on every progress tick (hundreds of frames at 680-exposure scale). Returning early avoids that.

> **On Task 18 / Task 20 coordination:** Task 20 calls `invalidateIngestFrameCache` (the exported helper from Task 18) — there is NO duplication of `invalidateQueries` calls between the two. The App listener owns the store write (Zustand stays out of `applyRemoteToCache`); the cache invalidation is a single function defined once in `applyRemoteToCache.ts` and called from both the switch arm (the path for any non-App caller) and the App listener. The STRUCTURAL frames (Task 19) are NOT short-circuited — they DO flow through `handleRemoteEvent` (they can be own-ops in E2), so their cache effect stays in `applyRemoteToCache`.

**Files:**
- Modify: `src/print/App.tsx` (the SSE `useEffect` ~31-42)
- Test: `test/App.ingestListener.test.tsx` (CREATE)

- [ ] **Step 1: Write the failing test** (driving the store via a dispatched MessageEvent through a mocked EventSource):

```tsx
// test/App.ingestListener.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { MemoryRouter } from "react-router-dom";
import { PrintApp } from "../src/print/App";
import { useAppState } from "../src/state";
import * as api from "../src/api";

// Minimal EventSource stub that captures listeners so the test can fire frames.
class FakeES {
  static last: FakeES | null = null;
  listeners: Record<string, ((e: MessageEvent) => void)[]> = {};
  url: string;
  constructor(url: string) { this.url = url; FakeES.last = this; }
  addEventListener(t: string, fn: (e: MessageEvent) => void) {
    (this.listeners[t] ??= []).push(fn);
  }
  removeEventListener() {}
  close() {}
  emit(type: string, data: unknown) {
    for (const fn of this.listeners[type] ?? []) fn({ data: JSON.stringify(data) } as MessageEvent);
  }
}

describe("App ingestInFlight SSE listener (Phase E1)", () => {
  beforeEach(() => {
    (globalThis as unknown as { EventSource: unknown }).EventSource = FakeES;
    useAppState.setState({ ingestInFlight: null });
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
  });
  afterEach(() => { useAppState.setState({ ingestInFlight: null }); vi.restoreAllMocks(); });

  function mount() {
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={["/experiments"]}><PrintApp /></MemoryRouter>
      </QueryClientProvider>,
    );
  }

  it("writes ingestInFlight on ingest_progress", () => {
    mount();
    FakeES.last!.emit("curation", {
      kind: "ingest_progress", entity_id: 7, entity_type: "experiment",
      payload: { experiment_id: 7, processed: 300, total: 680 },
    });
    expect(useAppState.getState().ingestInFlight?.[7]).toEqual({
      processed: 300, total: 680, status: "scanning",
    });
  });

  it("clears ingestInFlight on ingest_complete", () => {
    mount();
    FakeES.last!.emit("curation", { kind: "ingest_started", entity_id: 7, entity_type: "experiment", payload: { experiment_id: 7, processed: 0, total: 680 } });
    FakeES.last!.emit("curation", { kind: "ingest_complete", entity_id: 7, entity_type: "experiment", payload: { experiment_id: 7, processed: 680, total: 680 } });
    expect(useAppState.getState().ingestInFlight).toBeNull();
  });

  it("an ingest_* frame is short-circuited (does NOT reach handleRemoteEvent)", async () => {
    // Spy on the queue reconciler: an ingest frame must NOT trigger it (it is
    // broadcast-only, never an own-op; running it would roll back pending edits).
    const rc = await import("../src/lib/queue/replayCoordinator");
    const spy = vi.spyOn(rc, "handleRemoteEvent");
    mount();
    FakeES.last!.emit("curation", {
      kind: "ingest_progress", entity_id: 7, entity_type: "experiment",
      payload: { experiment_id: 7, processed: 300, total: 680 },
    });
    expect(spy).not.toHaveBeenCalled();
    // A NON-ingest frame still flows through (control).
    FakeES.last!.emit("curation", { kind: "peak_added", entity_id: 1, entity_type: "exposure", payload: {} });
    expect(spy).toHaveBeenCalledTimes(1);
  });
});
```

- [ ] **Step 2: Run it red** — `npm test -- test/App.ingestListener.test.tsx` → FAIL (no store write happens).

- [ ] **Step 3: Implement** — in `src/print/App.tsx`, extend the SSE `useEffect`. Import the store, the helper, and the two setters:

```tsx
import { useAppState } from "../state";
import { invalidateIngestFrameCache } from "../lib/queue/applyRemoteToCache";
```

Inside `PrintApp`, before the SSE effect:

```tsx
  const setIngestProgress = useAppState((s) => s.setIngestProgress);
  const clearIngestProgress = useAppState((s) => s.clearIngestProgress);
```

Replace the SSE listener body (lines 33-40) to short-circuit `ingest_*` (store write + shared cache helper + RETURN), then run the queue path for everything else:

```tsx
    es.addEventListener("curation", (e) => {
      try {
        const parsed = JSON.parse((e as MessageEvent).data as string) as SseEvent;
        // Ingest progress is broadcast-only (spec §9.3): never an own-op. Handle
        // it HERE and RETURN — do NOT feed it to handleRemoteEvent, whose queue
        // reconciliation would roll back + re-run every pending edit on each
        // progress tick. The store write lives here (Zustand is NOT imported in
        // applyRemoteToCache, which stays pure); the cache invalidation delegates
        // to the shared invalidateIngestFrameCache helper (defined + tested in
        // applyRemoteToCache.ts, Task 18) — no duplication.
        if (
          parsed.kind === "ingest_started" || parsed.kind === "ingest_progress" ||
          parsed.kind === "ingest_complete" || parsed.kind === "ingest_failed"
        ) {
          const p = (parsed as { payload?: { experiment_id?: number; processed?: number; total?: number } }).payload;
          const expId = p?.experiment_id;
          if (expId !== undefined) {
            if (parsed.kind === "ingest_started" || parsed.kind === "ingest_progress") {
              // E1 maps both pre-terminal frames to "scanning" (advisory header
              // label); refine to "analyzing" when the backend frame carries a
              // phase discriminator.
              setIngestProgress(expId, {
                processed: p?.processed ?? 0,
                total: p?.total ?? 0,
                status: "scanning",
              });
            } else {
              // Terminal (complete/failed): drop the in-flight entry; the
              // experiment's ingest_status (refetched below) is the resting truth.
              clearIngestProgress(expId);
            }
            // Delegate cache invalidation to the shared helper (defined in
            // applyRemoteToCache.ts) — the single source of truth for which
            // query keys the ingest frames affect.
            invalidateIngestFrameCache(qc, expId, parsed.kind === "ingest_complete");
          }
          return; // do NOT run the queue reconciler for a broadcast-only frame
        }
        handleRemoteEvent(parsed, qc, mc);
      } catch {
        // malformed frame, ignore
      }
    });
```

Add `setIngestProgress, clearIngestProgress` to the effect's dependency array (both stable selector results, so the EventSource lifetime is unchanged):

```tsx
  }, [qc, mc, setIngestProgress, clearIngestProgress]);
```

> Structural frames (Task 19: `sample_renamed`/`exposure_moved`/`sample_created`/`sample_split`/`grouping_flag_dismissed`) are NOT short-circuited — they fall through to `handleRemoteEvent` (they can be own-ops in E2), and their cache effect stays in `applyRemoteToCache`.

- [ ] **Step 4: Run it green** — `npm test -- test/App.ingestListener.test.tsx` → PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/App.tsx packages/HimalayaUI/frontend/test/App.ingestListener.test.tsx
git commit -m "feat(app): separate ingestInFlight SSE listener (store write outside pure cache path)"
```

---

## Task 21: ~~`persistence.ts` — `SCHEMA_VERSION` 4 → 5~~ — FOLDED INTO TASK 1b

**FOLDED into Task 1b.** The queue-op `SCHEMA_VERSION` 4→5 bump (so pre-deploy `update_sample` ops carrying the old `display_name` key drop on rehydrate) is part of the same atomic `display_name → name` collapse, so it lands in Task 1b alongside the consumer sweep (`test/displayNameCollapse.test.ts` asserts `SCHEMA_VERSION === 5`). It is the queue-op schema version, distinct from the Zustand `persist` `version` (which stays 5, `state.ts:499`, and is NOT touched). Nothing to do here.

---

## Task 22: Full-suite + build gate

Run the whole frontend test suite and the production build (which runs `tsc --noEmit`, the `lint:design` appearance guard, and `vite build`). Fix any regression the new work introduced (most likely: a stage-count assertion, the `/`→`/samples` assertion, or the `display_name` label sweep typechecking against `Sample`).

**Files:**
- No new files (fixes only)

- [ ] **Step 1: Run the full Vitest suite**

`npm test` (from `packages/HimalayaUI/frontend/`). Expected: green. Triage failures (the known ones, verified against live tests 2026-06-18):
- `test/AppRoutes.test.tsx` — the "bare / redirects to the corpus contact sheet (/samples)" spec (asserts `samples-page` after `renderRoutes("/")`) BREAKS once `/`→`/experiments` (Task 15). Update it to expect the experiments home ("Your beamtimes" / the experiments-home testid). Note `renderRoutes("/")` also appears at ~lines 221/228 in unrelated stale-path specs — re-check those land where intended after the redirect change. The `/index*` redirects still target `/samples` (unchanged).
- `test/CorpusTopbar.test.tsx:70` — "renders **two** stage-tabs (Samples, Series)" BREAKS once the Experiments tab is prepended (Task 16). Update it to three tabs (and add a `stage-tab-experiments` assertion). `CorpusTopbarFocus.test.tsx` only asserts `stage-tab-index` is absent (still true).
- The `display_name`/`name` sweep is already done in **Task 1b** (a sequenced task, not a Task-22 chore). If `git grep -n display_name packages/HimalayaUI/frontend/src` still returns a non-test hit here, a file was missed in Task 1b — fix it.

- [ ] **Step 2: Run the build gate**

`npm run build`. This runs `tsc --noEmit` + `lint:design` + `vite build`. Expected: pass. Triage:
- `lint:design` failure → an inline appearance utility leaked into a consumer (not a `print/ui/**` file). Move the appearance into the relevant primitive or swap to a token utility. The new primitives (`Dropdown`, `StatBar`) are in `print/ui/**` (exempt); the pages/components use only token utilities + named scale roles by construction.
- `tsc` failure on `display_name` / nullable `name` → a Task 1b straggler; finish the sweep.

- [ ] **Step 3: Commit any fixes**

```bash
git add packages/HimalayaUI/frontend/src packages/HimalayaUI/frontend/test
git commit -m "test/build: green full suite + production build for ingestion E1 shell"
```

---

## Self-Review

**What I verified against live source (read in full; every API below is confirmed, not assumed):**
- `api.ts` — `Experiment` (8-22, no `*_source`/scan cols yet → Task 1 adds them), `Sample` (31-38, carries BOTH `name: string|null` AND `display_name: string|null` → Task 1 collapses to non-null `name`), `updateExperiment`'s `Record<string,never>` stub (123-127 → Task 2 widens to `ExperimentPatch` — **NOT** `ExperimentGeometryPatch`; the canonical name is `ExperimentPatch`, E2 imports it), `updateSample`'s `display_name` patch key (134), `request<T>(method, path, body, opts)`.
- `lib/authOpts.ts` — **`authOpts` lives HERE, not in `queries.ts`**; signature `authOpts(username, clientId, clientOpId?)` with clientId POSITIONAL (Fix #4 applied: Task 12 imports from `../../lib/authOpts` and threads `CLIENT_ID`).
- `state.ts` — `partialize` (500-508), Zustand persist `version: 5` (499, NOT bumped), `LS_KEY = "himalaya-ui:state"` (24, = the persist `name` → the Task 3 test's literal key is correct), the `set((s)=>…)` idiom.
- `queries.ts` — `queryKeys.samples` = `["experiment", id ?? "none", "samples"]` (51-52); **so `loads` MUST be `["experiment", id ?? "none", "loads"]`** (Fix #7 applied — the plan's earlier `["experiment", 7, "loads"]` was right only for id=7 but the SIGNATURE must be `id ?? "none"`-tolerant). `experiment(id)` = `["experiment", id]` (50, no "config" segment → Fix #10: `experimentConfig` dropped as a dead key, `useExperimentDetail` alias dropped — callers use `useExperiment` directly). `useExperiment` exists (112, `enabled: id>0`).
- `ProgressBar.tsx` — props are `value` + **`total` (REQUIRED, NOT `max`)** + optional `label` (Fix #1 applied, Task 10).
- `Dot.tsx` — `DotTone = accent|success|muted|neutral` (NO `ok`/`danger`) (Fix #2 applied: failed→muted, scanning/analyzing→accent, complete→success, returns `DotTone`).
- `NoticePill.tsx` — `NoticePillTone = new|draft` ONLY (NO `warning`) (Fix #3 applied: Task 13 uses `Badge` for the review count instead).
- `Menu.tsx` (`open`/`options`/`onSelect`/`onClose`/`activeValue`/`className`, owner returns trigger focus), `Input.tsx` (`variant="title"`, `testId`, `value`/`onValueChange`, `...rest` carries `onKeyDown`/`aria-label`; note the wrapper div holds `data-testid`, the inner `<input>` gets the keydown — the tests' bubbling/`querySelector("input")` both work), `Card.tsx` (`as`/`interactive`/`padding` own-props; `onClick`/`data-testid` via `...rest`), `Button.tsx` (`accent`/`ghost` real; `data-testid` via `...props`), `EmptyState.tsx` (`title`/`body`/`action`/`as`), `Badge.tsx` (inline mono count), `PageFrame.tsx`, `ui/index.ts` — all confirmed.
- `applyRemoteToCache.ts` — `payload = remote.payload as Record<string,unknown>|undefined` (106), `default:` fires `peaks(id)`/`indices(id)` off `entity_id` (334-337); the existing `update_sample` arm (314-318) does `{...old, ...payload}` (handles a `name` patch unchanged). New `ingest_*`/structural arms placed before `default:`, reading `payload.experiment_id`.
- `lib/queue/types.ts` — **`SseEvent.kind` is already `string`** (148), `payload: unknown` (155). Fix #10 applied: NO kind-union edit; the `case` labels compile directly (Task 18 Step 4 reduced to a no-op note; Task 18 no longer touches `types.ts`).
- `App.tsx` — single `addEventListener("curation", …)` (33-40); `handleRemoteEvent` imported from `replayCoordinator` (8); Zustand NOT imported. Fix #9 applied: Task 20 short-circuits `ingest_*` (store write + `invalidateIngestFrameCache` call + `return`) BEFORE `handleRemoteEvent`, with a spy test proving the reconciler isn't called for an ingest frame. The cache invalidation delegates to the shared helper (no duplication).
- `CorpusTopbar.tsx` (`STAGES` two entries + `stage-tab-${id}` + `startsWith`), `AppRoutes.tsx` (`/`→`/samples` at 116, outside-shell redirects). Live-test breakage confirmed: `CorpusTopbar.test.tsx:70` asserts "two stage-tabs"; `AppRoutes.test.tsx` asserts `/`→`samples-page` — both updated in Task 22.
- `persistence.ts` — `SCHEMA_VERSION = 4` (14) → bumped to 5 in Task 1b (Task 21 folded).
- **Fix #5 (Task 17):** `git grep -n "PhaseStrip\|phase-bar\|phaseDistribution\|nav-home-phase" src/print/shell/` returns NOTHING; NavModal references `phase` only via `display_name`/`name` search predicates. **Task 17 DROPPED** (a removal test for a never-rendered element is vacuous). The `PhaseStrip` primitive exists in `src/print/ui/` but no shell file consumes it.
- **Fix #6 (Task 1b):** the full `display_name → name` consumer list was captured via `git grep -n display_name src/` (12 non-test source files) and promoted to a single dedicated sequenced task (1b) that converts every consumer + bumps `SCHEMA_VERSION`, leaving `tsc` green at its commit (Task 1 + 1b are an intentional red→green pair).

**Canonical contract (E1 OWNS; E2 consumes; pinned spec §8.8):**
- `api.ts`: `GroupingFlag` (merge/split/null), `LoadExposure`, `LoadSample` (the (load,slot) coord + `flag`), nested `Load { load_id, …, samples: LoadSample[] }` — **Fix #7 applied: the plan's earlier FLAT `Load` is replaced by this nested shape.** `Sample.name` non-null. `updateExperiment` patch widened to `ExperimentPatch` (canonical name; includes name, description, geometry ×6, file patterns ×3 — **E2 imports `ExperimentPatch` from `api.ts`, must NOT redefine it**).
- `queryKeys.loads(id)` = `["experiment", id ?? "none", "loads"]`; `useLoads(id)` gated `id>0`.
- Ingest SSE frame = payload-wrapped `{ kind, payload:{experiment_id, processed, total} }`; arms read top-level `kind` + `payload.experiment_id`.
- Structural payloads: `exposure_moved {sample_id, from_sample_id, experiment_id}`, `sample_renamed {name, experiment_id}`, `sample_created {experiment_id}`, `sample_split {new_sample_id, exposure_ids, experiment_id}`, `grouping_flag_dismissed {flag_kind, merge_with_sample_id?, experiment_id}`. **NO `sample_merged`** (Fix #8 applied: dropped; the prior `to_sample_id` field corrected to `sample_id`; `sample_created` arm added before `default:`).

**Honest open questions / call-outs (decisions a reviewer should confirm):**
1. **`/` → `/experiments` is more disruptive than one assertion.** `AppRoutes.test.tsx` has ≥1 spec hard-asserting `/`→`samples-page`, plus two more `renderRoutes("/")` call sites in stale-path specs. Task 22 flags them; a reviewer should confirm the experiments-home renders a stable testid for those updated assertions (Task 11 should expose one — see #2).
2. **`ExperimentsHomePage` testid for "Your beamtimes".** The AppRoutes/Home tests assert on the literal text "Your beamtimes". If that copy changes, add an explicit page testid (e.g. `experiments-home`) and assert on it instead — text-coupling is brittle. (Left as-is to match the existing AppRoutes test idiom, which asserts `samples-page` by testid; consider adding `data-testid="experiments-home"` to the Home PageFrame.)
3. **Task 13 test uses `vi.spyOn(api, "listLoads")` — no test-only prop** (`__testReviewCount` removed). The `ExperimentCorpusPage` component has no production seam; the review count derives purely from `useLoads(expId)` which the test controls via the spy. The `retry: false` QueryClient ensures the mock resolves cleanly without retries.
4. **App-listener / `applyRemoteToCache` cache invalidation: single shared helper.** Task 18 exports `invalidateIngestFrameCache(qc, expId, isComplete)`; Task 20's App listener calls it. There is no duplication — the function is defined once and called from two sites. The store write stays exclusively in the App listener.
5. **`/series` axis** stays the existing folio surface (spec §7). `ExperimentTopNav`'s Series link points at `/series` under `CorpusShell`, so clicking it from inside an experiment leaves `ExperimentShell` — intended IA (Series is cross-cutting).
6. **No-placeholder check:** every referenced type (`IngestProgress`, nested `Load`/`LoadSample`/`LoadExposure`/`GroupingFlag`, `ExperimentPatch` [canonical — E2 imports this, not `ExperimentGeometryPatch`], `GeometrySource`, `IngestStatus`, `ValidatePathResponse`, `PathSuggestResponse`, `DropdownProps`, `StatBarStat`, `DirectoryPickerFieldProps`) is defined in some task; every E2-deferred slot names the E2 component + renders a labelled placeholder; every step shows the command + expected red/green.
