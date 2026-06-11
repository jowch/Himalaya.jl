# Editable Sample Tags Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make sample tags editable — a per-sample "Manage tags" modal with corpus-aware key/value suggestion dropdowns, a backend `edit_tag` path, and a scoping upsert — closing the LO-TAGDUP P1 (a duplicate manual/scoping tag that can't be individually removed because the loupe resolves removal by key+value, not id).

**Architecture:** Five slices, smallest-correct-first. Slice 1 (id-exact delete) is the standalone literal P1 fix. Slices 2-4 add the edit endpoint, counts, and the modal. Slice 5 (scoping upsert, the risk locus) goes last. Backend mirrors the existing `add_tag`/`remove_tag` route + event machinery; frontend mirrors the existing `addCorpusSampleTagMutator`/`removeCorpusSampleTagMutator` queue chain; the modal reuses the custom-index `ModalShell`.

**Tech Stack:** Julia/Oxygen.jl/SQLite backend; React/TypeScript/TanStack Query/Zustand frontend with the project mutation queue; Vitest + Playwright; the project's `with_idempotency`/`apply_event!`/SSE machinery.

**Spec:** `docs/superpowers/specs/2026-06-11-editable-sample-tags-design.md` (read it first — it carries the rationale for every decision).

**Gate (run from `packages/HimalayaUI/frontend` for frontend slices):** `npm run lint:design` → `npx tsc --noEmit -p tsconfig.build.json` → `npx vitest run <touched files>` → (full `npx vitest run` + `npm run e2e` + `npm run build` before declaring a slice done). Backend: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'` (slow; capture once).

---

## Review-driven corrections (3 reviewers, folded 2026-06-11 — READ BEFORE EXECUTING)

All three plan reviewers returned READY-WITH-FIXES (no architectural blocker). Apply these alongside each task:

**Exact test-file paths** (the plan's generic placeholders resolved):
- T3 + T10 backend tests → `packages/HimalayaUI/test/test_routes_samples.jl` (existing tag-route tests at ~3-54; batch tests at ~114-296). Bare-`SQLite.DB()` harness: `db = SQLite.DB(); HimalayaUI.create_schema!(db); …; with_test_server(db) do port, base`. Every new `@testset` opens its own `with_test_server`. Build headers inline `["Content-Type"=>"application/json","X-Username"=>"alice"]` — **there is no `ct` variable** (the plan's test snippets use `ct`; define it or inline). Seed samples/tags via `HimalayaUI.create_sample!` etc.; all `$sid/$T/$A/$B/$S/$X` placeholders must be real seeded ids.
- T6 backend test → `packages/HimalayaUI/test/test_picker_routes.jl` (corpus test ~150-193). This file uses the `mktempdir()` + `_setup_analyzed_exposure(tmp)` + `with_test_server(ctx.db)` harness — match the file's local convention, not the bare-DB one.
- T1 → `test/print-pages/loupeAdapters.test.ts` — **rewrite/remove the existing `describe("toLoupeTags / findSampleTagId", …)` block (~125-138)**: deleting `findSampleTagId` without fixing that block breaks the build.
- T2 → `test/print-pages/LoupePage.test.tsx` (already mocks `useRemoveCorpusSampleTag` → `removeTagMutate`; reuse it as the spy).
- T5 → new `test/queue/editTag.test.ts` **and** add an `edit_tag` arm assertion to `test/queue/mutatorRegistry.test.ts` (both `resolveMutator` sampleId-only and `resolveMutatorForEvent`).
- T6 frontend → `test/proposeOrdering.test.ts`.

**T2 delete-wiring is deeper than written (BIG):** the loupe tags render through `LoupeSidePanel` → `<TagList>` → `<TagPill>` → `<Chip>`, NOT a per-pill loop in `LoupeSidePanel`. `TagList`'s `onRemove(t)` passes a `Tag`, and `Chip`'s removable variant hardcodes `aria-label="Remove"`. To make removal id-exact AND name the tag, T2 must: add a `LoupeTag`-aware path through `TagList`/`TagPill` that carries `id` and a computed per-tag remove label, and thread a `removeLabel` prop into `TagPill`→`Chip` (Chip/TagPill/TagList are in `src/print/ui/**`, design-guard-exempt). Sequence: T1 (adapter) → this threading → id-exact `onRemove(id)`.

**Backend code fixes:**
- T10: **hoist the `exp_id` SELECT to the top of the `for t in entries` loop body** (all three branches use it; today it's computed mid-loop after the INSERT). Ensure all three branches fall through to ONE `push!(created, …)` using the resolved `tag_id`/`value` so the batch response contract (one entry per input row) is preserved — read the current `created` builder first.
- T6: the experiment-scoped sibling query is JOINed and aliased — write `SELECT t.key, t.value, COUNT(DISTINCT t.sample_id) AS count … WHERE s.experiment_id = ? GROUP BY t.key, t.value`. Add `:count => Int(r.count)` to BOTH routes' JSON comprehensions.
- T3: a spec wrinkle — spec §Testing's "round-trip parity" line is **superseded** by §Architecture-#2 (sample_tags is a base table, not a rebuilt view); the plan's user_actions test is correct. **Strengthen** it to also assert no view-write + the broadcast fires (not just `+1` row).

**Test-coverage fixes:**
- **T5 contract test (false-green risk — the documented SSE-wins-partial trap):** the single `onMutate` test is insufficient. Add own-op SSE-confirmation (`synthesizeFromSse` resolves a synthesized partial) and foreign-event replay assertions, and treat it as the **pair** of the backend `apply_event!` test (six-layer rule).
- **T7 (TagSuggest) and T8 (ManageTagsModal) are too coarse** — split each into 2-3 red→green→commit sub-tasks: T7a "combobox + listbox ARIA + keyboard nav", T7b "counts + create-as-typed"; T8a "shell + rows + edit/delete wiring", T8b "add-row + duplicate-key rejection + focus-restore-to-trigger + `lib/announce`". Add the currently-missing tests for `onRemove`, focus-restore-on-close (ModalShell does NOT provide it — the modal self-manages it), the live-announce text, and the create-as-typed keyboard path.
- T10: add a batch idempotency-replay test (same `client_op_id` → no new rows/events) and confirm the duplicate-`sample_id`-in-batch guard still holds.

**Design-guard:** `TagSuggest.tsx` is in `src/print/ui/**` (guard-exempt — may author appearance). `ManageTagsModal.tsx` is in `src/print/components/` (**scanned** — token-only classes; no `text-[…]`/`rounded-[…]`/raw color, mirror `CustomIndexModal.tsx`).

**Naming:** the spec says `useEditSampleTag`/`editSampleTagMutator`; the plan (authoritative) uses `useEditCorpusSampleTag`/`editCorpusSampleTagMutator` consistently. T9 must instantiate `const editTag = useEditCorpusSampleTag(sample.id)` (the plan implies but doesn't write it). `synthesizeFromSse`'s hardcoded `source:"manual"` is a faint-marker-only inaccuracy for a scoping-edited tag — acceptable, comment it.

**Confirmed correct (no change):** `findSampleTagId` has no caller outside the loupe (safe to delete); `resolveMutator` tri-scope routes a sampleId-only `edit_tag` to the corpus mutator; `applyRemoteToCache` `case "edit_tag":` added to the `add_tag`/`remove_tag` block inherits the right 3 invalidations; all Task-5 imports are real; `proposeOrdering` genuinely ignores `count`; `announce`, `ModalShell`/`ModalHead`/`ModalFooter`, `TagEditor.commit()` all exist as used.

## File structure

**Slice 1 — id-exact delete (frontend only):**
- Modify `packages/HimalayaUI/frontend/src/print/pages/loupeAdapters.ts` — `toLoupeTags` emits a loupe-local `LoupeTag` view-model carrying `id`+`source`; delete `findSampleTagId`.
- Modify `packages/HimalayaUI/frontend/src/print/components/LoupeSidePanel.tsx` + `LoupePage.tsx` — pills carry their `id`; `onRemoveTag` receives the id.
- Tests: `test/print-pages/loupeAdapters.test.ts`, `test/LoupePage.*` (existing loupe tests).

**Slice 2 — edit_tag path (backend + frontend queue):**
- Modify `packages/HimalayaUI/src/routes_samples.jl` — `PATCH /api/samples/{id}/tags/{tag_id}`.
- Modify `packages/HimalayaUI/src/events.jl` — `edit_tag` no-op dispatcher guard.
- Modify frontend `src/api.ts` (`editSampleTag`), `src/queries.ts` (`useEditSampleTag`), `src/lib/queue/types.ts` (`OpKind`), `src/lib/queue/mutators/trivial.ts` (`editCorpusSampleTagMutator`), `src/lib/queue/mutatorRegistry.ts` (`edit_tag` arms), `src/lib/queue/applyRemoteToCache.ts` (`edit_tag` case).
- Tests: `packages/HimalayaUI/test/test_routes_samples.jl` (or the nearest tag-route test), an `edit_tag` `apply_event!` test, `frontend/test/queue/*` contract test.

**Slice 3 — counted corpus tags:**
- Modify `packages/HimalayaUI/src/routes_picker.jl` — counted `GET /api/sample-tags` + experiment sibling.
- Modify frontend `src/api.ts` (`SampleTagPair` gains `count`).
- Tests: backend picker-route test; `proposeOrdering` UNCHANGED (a test pinning that it ignores `count`).

**Slice 4 — TagSuggest + ManageTagsModal:**
- Create `frontend/src/print/ui/TagSuggest.tsx` (+ `.stories.tsx`, `test/print-ui/TagSuggest.test.tsx`).
- Create `frontend/src/print/components/ManageTagsModal.tsx` (+ `test/print-components/ManageTagsModal.test.tsx`).
- Modify `src/print/ui/TagEditor.tsx` + `TagList.tsx` — thread existing-keys for inline duplicate-key rejection.
- Modify `LoupeSidePanel.tsx`/`LoupePage.tsx` — "Manage" affordance + modal wiring + `lib/announce`.
- E2E `e2e/loupe-tags.spec.ts` (mocked).

**Slice 5 — scoping upsert:**
- Modify `packages/HimalayaUI/src/routes_samples.jl` — batch upsert.
- Tests: batch upsert (insert/update/no-op/idempotency) in the batch-route test file.

---

## Slice 1 — id-exact inline delete (the standalone LO-TAGDUP correctness fix)

### Task 1: `LoupeTag` view-model carries the tag id

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/print/pages/loupeAdapters.ts`
- Test: `packages/HimalayaUI/frontend/test/print-pages/loupeAdapters.test.ts`

- [ ] **Step 1: Write the failing test**

In `test/print-pages/loupeAdapters.test.ts`, add:

```ts
import { describe, it, expect } from "vitest";
import { toLoupeTags } from "../../src/print/pages/loupeAdapters";
import type { SampleTag } from "../../src/api";

describe("toLoupeTags (LO-TAGDUP)", () => {
  it("carries the tag id and source so byte-identical duplicates are distinguishable", () => {
    const tags: SampleTag[] = [
      { id: 3, key: "dose", value: "10", source: "manual" },
      { id: 7, key: "dose", value: "10", source: "scoping" },
    ];
    const out = toLoupeTags(tags);
    expect(out).toEqual([
      { id: 3, key: "dose", value: "10", source: "manual" },
      { id: 7, key: "dose", value: "10", source: "scoping" },
    ]);
    // The two pills are now distinguishable by id even though key+value match.
    expect(out[0]!.id).not.toBe(out[1]!.id);
  });

  it("omits an empty value (keeps the exactOptionalPropertyTypes contract)", () => {
    const out = toLoupeTags([{ id: 1, key: "k", value: "", source: "manual" }]);
    expect(out[0]).toEqual({ id: 1, key: "k", value: undefined, source: "manual" });
  });
});
```

- [ ] **Step 2: Run it, verify it fails**

Run: `npx vitest run test/print-pages/loupeAdapters.test.ts`
Expected: FAIL — `toLoupeTags` currently returns `{key, value?}` without `id`/`source`.

- [ ] **Step 3: Implement**

In `loupeAdapters.ts`, replace `toLoupeTags` and **delete** `findSampleTagId`:

```ts
import type { Tag } from "../ui";

/** Loupe-local tag view-model: the ui `Tag` ({key, value?}) plus the backend
 *  identity (`id`) and provenance (`source`). The `id` is what makes two
 *  byte-identical pills (manual dose=10 + scoping dose=10) individually
 *  addressable for removal/edit — `ui/tag.ts` is deliberately NOT widened. */
export interface LoupeTag extends Tag {
  id: number;
  source: string;
}

/** SampleTag[] → LoupeTag[]; omit empty value (exactOptionalPropertyTypes). */
export function toLoupeTags(tags: SampleTag[]): LoupeTag[] {
  return tags.map((t) =>
    t.value
      ? { id: t.id, key: t.key, value: t.value, source: t.source }
      : { id: t.id, key: t.key, value: undefined, source: t.source },
  );
}
```

(Remove the old `findSampleTagId` function entirely — Task 2 removes its last caller.)

- [ ] **Step 4: Run it, verify it passes**

Run: `npx vitest run test/print-pages/loupeAdapters.test.ts`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/print/pages/loupeAdapters.ts packages/HimalayaUI/frontend/test/print-pages/loupeAdapters.test.ts
git commit -m "LO-TAGDUP: loupe tag view-model carries id+source (delete findSampleTagId)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 2: remove a tag by its exact id

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/print/pages/LoupePage.tsx` (the `handleRemoveTag` callback + the `toLoupeTags(sample.tags)` pass-down)
- Modify: `packages/HimalayaUI/frontend/src/print/components/LoupeSidePanel.tsx` (the tags section: each pill's remove fires with its `id`)
- Test: existing loupe page test (`test/LoupePage.test.tsx` or the nearest) + a new assertion

- [ ] **Step 1: Write the failing test**

Add to the loupe page test (the file that mounts `LoupePage` with a mocked sample carrying two same-key tags):

```tsx
it("removing the second of two byte-identical dose pills deletes tag id 7, not id 3 (LO-TAGDUP)", async () => {
  // sample seeded with tags [{id:3,key:dose,value:10,source:manual},
  //                          {id:7,key:dose,value:10,source:scoping}]
  renderLoupe(/* sample with the duplicate */);
  const removes = screen.getAllByRole("button", { name: /remove dose/i });
  fireEvent.click(removes[1]!); // the scoping twin
  expect(removeTagSpy).toHaveBeenCalledWith(7);
  expect(removeTagSpy).not.toHaveBeenCalledWith(3);
});
```

(Wire `removeTagSpy` to the `useRemoveCorpusSampleTag().mutate` mock the test already stubs.)

- [ ] **Step 2: Run it, verify it fails**

Run: `npx vitest run test/LoupePage.test.tsx`
Expected: FAIL — today `handleRemoveTag` re-derives the id via `findSampleTagId(sample.tags, t)`, which returns the FIRST key+value match (id 3) for both pills.

- [ ] **Step 3: Implement**

In `LoupePage.tsx`, replace the `handleRemoveTag` callback so it takes an id directly:

```tsx
const handleRemoveTag = useCallback((tagId: number) => {
  removeTag.mutate(tagId);
}, [removeTag]);
```

Pass `toLoupeTags(sample.tags)` to the side panel (it already does — now the rows carry `id`). In `LoupeSidePanel.tsx`, change the tags-section prop type to `LoupeTag[]` and make each pill's remove call `onRemoveTag(tag.id)` and its `aria-label` name the tag: `aria-label={\`Remove ${tag.key}${tag.value ? ` ${tag.value}` : ""}\`}`. Remove the `findSampleTagId` import.

- [ ] **Step 4: Run it, verify it passes**

Run: `npx vitest run test/LoupePage.test.tsx test/print-pages/loupeAdapters.test.ts`
Expected: PASS. Also run `npx tsc --noEmit -p tsconfig.build.json` (verifies no remaining `findSampleTagId` reference).

- [ ] **Step 5: Render-verify + commit**

Render-verify at :5182 on the sample-3 loupe: the two `dose=10` pills are individually removable (delete the scoping one, the manual one remains). Then:

```bash
git add packages/HimalayaUI/frontend/src/print/pages/LoupePage.tsx packages/HimalayaUI/frontend/src/print/components/LoupeSidePanel.tsx packages/HimalayaUI/frontend/test/LoupePage.test.tsx
git commit -m "LO-TAGDUP: remove a sample tag by exact id, not key+value match

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

> **Slice 1 closes the LO-TAGDUP P1 literally.** Re-run the full gate (`lint:design` + `tsc` + `vitest` + `e2e` + `build`) before moving on. The remaining slices deliver the full "editable tags" feature.

---

## Slice 2 — the `edit_tag` path (backend route + event + frontend queue chain)

### Task 3: backend `PATCH /api/samples/{id}/tags/{tag_id}`

**Files:**
- Modify: `packages/HimalayaUI/src/routes_samples.jl` (add after the `@delete .../tags/{tag_id}` handler, ~line 195)
- Test: `packages/HimalayaUI/test/` — the tag-route test file (grep for `"/api/samples/"` + `tags` in `test/`)

- [ ] **Step 1: Write the failing test**

Add to the tag-route test:

```julia
@testset "PATCH /api/samples/:id/tags/:tag_id edits in place" begin
    # seed a sample with one tag id=T
    r = HTTP.request("PATCH", "$base/api/samples/$sid/tags/$T",
        ["Content-Type"=>"application/json"],
        JSON3.write(Dict(:value => "12")); status_exception=false)
    @test r.status == 200
    body = JSON3.read(String(r.body))
    @test body.id == T && body.value == "12"
    # the row updated in place, id stable
    rows = DBInterface.execute(db, "SELECT value FROM sample_tags WHERE id=?", [T]) |> Tables.rowtable
    @test rows[1].value == "12"
end

@testset "PATCH rejects a duplicate-key collision with 409" begin
    # sample has dose=10 (id A) and temp=25 (id B); editing B's key to "dose" collides
    r = HTTP.request("PATCH", "$base/api/samples/$sid/tags/$B",
        ["Content-Type"=>"application/json"],
        JSON3.write(Dict(:key => "dose")); status_exception=false)
    @test r.status == 409
end

@testset "PATCH a non-matching (tag_id,sample_id) is 404" begin
    r = HTTP.request("PATCH", "$base/api/samples/$sid/tags/999999",
        ["Content-Type"=>"application/json"],
        JSON3.write(Dict(:value => "x")); status_exception=false)
    @test r.status == 404
end
```

- [ ] **Step 2: Run it, verify it fails**

Run: `julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")' 2>&1 | tee /tmp/jl.out` then grep the testset. Expected: FAIL (route not defined → 404 on the edit itself / no 409).

- [ ] **Step 3: Implement the route**

In `routes_samples.jl`, after the `@delete .../tags/{tag_id}` handler:

```julia
    @patch "/api/samples/{id}/tags/{tag_id}" function(req::HTTP.Request, id::Int, tag_id::Int)
        db   = current_db()
        body = json(req)
        # Partial update: key and/or value, each a string if present.
        has_key = haskey(body, :key); has_val = haskey(body, :value)
        (has_key || has_val) ||
            return HTTP.Response(400, ["Content-Type"=>"application/json"],
                JSON3.write(Dict(:error => "provide at least one of: key, value")))
        if (has_key && !(body.key isa AbstractString)) ||
           (has_val && !(body.value isa AbstractString))
            return HTTP.Response(400, ["Content-Type"=>"application/json"],
                JSON3.write(Dict(:error => "key and value must be strings")))
        end
        return with_idempotency(db, req) do
            cur = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, key, value, source FROM sample_tags
                 WHERE id = ? AND sample_id = ?", [tag_id, id]))
            isempty(cur) && return HTTP.Response(404, ["Content-Type"=>"application/json"],
                JSON3.write(Dict(:error => "tag not found")))
            new_key = has_key ? String(body.key) : String(cur[1].key)
            new_val = has_val ? String(body.value) : String(cur[1].value)
            # Single-valued-key rule: a key-edit must not collide with another
            # tag on the same sample (id<>? excludes the row being edited).
            if has_key
                coll = Tables.rowtable(DBInterface.execute(db,
                    "SELECT 1 FROM sample_tags
                     WHERE sample_id = ? AND key = ? AND id <> ?",
                    [id, new_key, tag_id]))
                isempty(coll) || return HTTP.Response(409, ["Content-Type"=>"application/json"],
                    JSON3.write(Dict(:error => "sample already has a '$new_key' tag")))
            end
            srows = Tables.rowtable(DBInterface.execute(db,
                "SELECT experiment_id FROM samples WHERE id = ?", [id]))
            exp_id = (isempty(srows) || ismissing(srows[1].experiment_id)) ?
                     nothing : Int(srows[1].experiment_id)
            DBInterface.execute(db,
                "UPDATE sample_tags SET key = ?, value = ? WHERE id = ? AND sample_id = ?",
                [new_key, new_val, tag_id, id])
            result = apply_event!(InTransaction(), db, req;
                kind = "edit_tag",
                entity_type = "sample", entity_id = id,
                payload = Dict(:tag_id => tag_id, :key => new_key, :value => new_val,
                               :experiment_id => exp_id))
            _enqueue_broadcast_from_result!(result, "edit_tag", "sample", id)
            HTTP.Response(200, ["Content-Type"=>"application/json"],
                JSON3.write(Dict(:id => tag_id, :key => new_key, :value => new_val,
                                 :source => String(cur[1].source))))
        end
    end
```

- [ ] **Step 4: Add the `edit_tag` dispatcher guard**

In `events.jl`, alongside the `add_tag`/`remove_tag` no-op returns (~line 406):

```julia
    kind == "edit_tag" && return nothing
```

- [ ] **Step 5: Add the `apply_event!` coverage test**

```julia
@testset "edit_tag is a non-view log event (user_actions row, no view write)" begin
    n0 = only(DBInterface.execute(db, "SELECT COUNT(*) c FROM user_actions") |> Tables.rowtable).c
    HTTP.request("PATCH", "$base/api/samples/$sid/tags/$T",
        ["Content-Type"=>"application/json"], JSON3.write(Dict(:value=>"13")))
    n1 = only(DBInterface.execute(db, "SELECT COUNT(*) c FROM user_actions") |> Tables.rowtable).c
    @test n1 == n0 + 1
end
```

- [ ] **Step 6: Run, verify pass, commit**

Run the suite, grep the new testsets → PASS.

```bash
git add packages/HimalayaUI/src/routes_samples.jl packages/HimalayaUI/src/events.jl packages/HimalayaUI/test/<tag-route-test>.jl
git commit -m "edit_tag: PATCH .../tags/:id route + non-view event (single-valued-key 409)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 4: frontend `api.editSampleTag` + `OpKind`

**Files:**
- Modify: `src/api.ts`, `src/lib/queue/types.ts`
- Test: covered by Task 5's queue test (no standalone test for the thin fetcher)

- [ ] **Step 1: add the fetcher** to `src/api.ts` next to `removeSampleTag`:

```ts
export const editSampleTag = (id: number, tag_id: number, patch: { key?: string; value?: string }, opts?: AuthOpts) =>
  request<SampleTag>("PATCH", `/api/samples/${id}/tags/${tag_id}`, patch, opts);
```

- [ ] **Step 2: add `edit_tag` to `OpKind`** in `src/lib/queue/types.ts` (the union at line 41):

```ts
  | "add_tag" | "remove_tag" | "edit_tag"
```

- [ ] **Step 3: typecheck + commit**

Run: `npx tsc --noEmit -p tsconfig.build.json` → clean.

```bash
git add packages/HimalayaUI/frontend/src/api.ts packages/HimalayaUI/frontend/src/lib/queue/types.ts
git commit -m "edit_tag: api.editSampleTag fetcher + OpKind member

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 5: `editCorpusSampleTagMutator` + registry + cache + `useEditSampleTag`

**Files:**
- Modify: `src/lib/queue/mutators/trivial.ts`, `src/lib/queue/mutatorRegistry.ts`, `src/lib/queue/applyRemoteToCache.ts`, `src/queries.ts`
- Test: `frontend/test/queue/editTag.test.ts` (new) — paired with the backend per the contract-testing rule

- [ ] **Step 1: Write the failing queue test**

```ts
import { describe, it, expect } from "vitest";
import { editCorpusSampleTagMutator } from "../../src/lib/queue/mutators/trivial";
import { QueryClient } from "@tanstack/react-query";
import { queryKeys } from "../../src/queries";

describe("editCorpusSampleTagMutator", () => {
  it("optimistically patches the tag in place by id (NOT key+value)", () => {
    const qc = new QueryClient();
    qc.setQueryData(queryKeys.corpusSamples, [
      { id: 1, tags: [
        { id: 3, key: "dose", value: "10", source: "manual" },
        { id: 7, key: "dose", value: "10", source: "scoping" },
      ] },
    ]);
    editCorpusSampleTagMutator.onMutate(
      { sampleId: 1, tagId: 7, key: "dose", value: "12" } as any, qc,
    );
    const tags = (qc.getQueryData(queryKeys.corpusSamples) as any)[0].tags;
    expect(tags.find((t: any) => t.id === 7).value).toBe("12"); // edited
    expect(tags.find((t: any) => t.id === 3).value).toBe("10"); // untouched
  });
});
```

- [ ] **Step 2: Run, verify it fails** (`editCorpusSampleTagMutator` undefined).

- [ ] **Step 3: Add the mutator** in `trivial.ts` after `removeCorpusSampleTagMutator`:

```ts
export type EditCorpusSampleTagInput = { tagId: number; key: string; value: string };
type EditCorpusSampleTagScope = BaseScope & { sampleId: number };

export const editCorpusSampleTagMutator: Mutator<
  EditCorpusSampleTagInput, EditCorpusSampleTagScope, SampleTag
> = {
  kind: "edit_tag",
  onMutate: (p, qc): RollbackContext => {
    const key = queryKeys.corpusSamples;
    const prev = qc.getQueryData<CorpusSample[]>(key);
    if (prev) {
      qc.setQueryData<CorpusSample[]>(key, prev.map((s) =>
        s.id === p.sampleId
          ? { ...s, tags: s.tags.map((t) =>
              t.id === p.tagId ? { ...t, key: p.key, value: p.value } : t) }
          : s,
      ));
    }
    return { restore: () => { if (prev !== undefined) qc.setQueryData(key, prev); } };
  },
  request: (p) => api.editSampleTag(p.sampleId, p.tagId, { key: p.key, value: p.value }, buildAuthOpts(p)),
  // id is stable across the edit, so the optimistic row already matches the
  // server row — no placeholder reconciliation needed (unlike add).
  onSuccess: () => {},
  // CORRECTED SIGNATURE (matches types.ts:215 + addSampleTagMutator.synthesizeFromSse,
  // trivial.ts:159): two args (remote, base); payload lives at remote.payload;
  // guard the missing field and spread ...base.
  synthesizeFromSse: (remote, base) => {
    const p = remote.payload as Record<string, unknown>;
    if (p.tag_id === undefined) return undefined;
    return { ...base, id: p.tag_id as number, key: p.key as string,
             value: p.value as string, source: "manual" };
  },
};
```

- [ ] **Step 4: Register it** in `mutatorRegistry.ts` — add an `edit_tag` arm to `resolveMutator` (mirror the corpus `add_tag` arm; it routes on `sampleId`-only scope) AND to `resolveMutatorForEvent` (so foreign `edit_tag` frames replay). Import `editCorpusSampleTagMutator`.

```ts
    case "edit_tag":
      return editCorpusSampleTagMutator;
```

- [ ] **Step 5: Cache invalidation** — in `applyRemoteToCache.ts`, add `case "edit_tag":` to the SAME block as `add_tag`/`remove_tag` (line 357), so it falls into the identical sample-vs-exposure invalidation (corpusSamples + corpusSampleTags + corpusPickerSamples). One line:

```ts
    case "add_tag":
    case "remove_tag":
    case "edit_tag": {
```

- [ ] **Step 6: the hook** in `src/queries.ts`. `useQueueMutation` takes **positional** `(mutator, scope)` args (NOT a `{kind}` object) — mirror `useRemoveCorpusSampleTag` exactly (read it at `src/queries.ts:605`), pulling `username`/`clientId` from `useAppState` the same way it does:

```ts
export function useEditCorpusSampleTag(sampleId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(editCorpusSampleTagMutator, { sampleId, username, clientId: CLIENT_ID });
  return { ...inner, mutate: (a: { tagId: number; key: string; value: string }) => inner.mutate(a) };
}
```

(Use whatever `username`/`clientId` source `useRemoveCorpusSampleTag` uses — copy its exact scope construction so routing through `resolveMutator` is identical.)

- [ ] **Step 7: Run the queue test + tsc + commit**

Run: `npx vitest run test/queue/editTag.test.ts` → PASS; `npx tsc --noEmit -p tsconfig.build.json` → clean.

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts packages/HimalayaUI/frontend/src/lib/queue/mutatorRegistry.ts packages/HimalayaUI/frontend/src/lib/queue/applyRemoteToCache.ts packages/HimalayaUI/frontend/src/queries.ts packages/HimalayaUI/frontend/test/queue/editTag.test.ts
git commit -m "edit_tag: queue mutator + registry arms + cache invalidation + useEditCorpusSampleTag

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Slice 3 — counted `GET /api/sample-tags`

### Task 6: count per (key,value); proposeOrdering stays unchanged

**Files:**
- Modify: `packages/HimalayaUI/src/routes_picker.jl` (both `/api/sample-tags` and `/api/experiments/{eid}/sample-tags`)
- Modify: `frontend/src/api.ts` (`SampleTagPair` gains `count?: number`)
- Test: backend picker test asserts the count; `frontend/test/.../proposeOrdering.test.ts` pins that ranking is unchanged

- [ ] **Step 1: Write the failing backend test**

```julia
@testset "GET /api/sample-tags returns per-(key,value) sample count" begin
    body = JSON3.read(String(HTTP.get("$base/api/sample-tags").body))
    row = first(filter(p -> p.key == "dose" && p.value == "10", body))
    @test haskey(row, :count) && row.count >= 1
end
```

- [ ] **Step 2: Run, verify fail** (no `count` key).

- [ ] **Step 3: Implement** — in `routes_picker.jl`, change the corpus query from `SELECT DISTINCT key, value FROM sample_tags ...` to:

```julia
"SELECT key, value, COUNT(DISTINCT sample_id) AS count
 FROM sample_tags GROUP BY key, value ORDER BY key, value"
```

Apply the identical change to the experiment-scoped sibling (add a `WHERE experiment_id`/join as that route already does). Add `count` to the JSON row.

- [ ] **Step 4: Frontend type** — in `src/api.ts`, `SampleTagPair` gains `count?: number;`.

- [ ] **Step 5: Pin proposeOrdering is unchanged** — add to its test:

```ts
it("ranks ordering keys by number of distinct values, ignoring the sample count field", () => {
  // corpus: keyA has 3 distinct values (low total count), keyB has 1 value (high count)
  const pairs = [
    { key: "keyA", value: "1", count: 1 }, { key: "keyA", value: "2", count: 1 },
    { key: "keyA", value: "3", count: 1 }, { key: "keyB", value: "x", count: 99 },
  ] as any;
  expect(proposeOrdering(pairs, samples).orderingKey).toBe("keyA"); // distinct-value winner, not count
});
```

- [ ] **Step 6: Run both suites, commit**

```bash
git add packages/HimalayaUI/src/routes_picker.jl packages/HimalayaUI/frontend/src/api.ts packages/HimalayaUI/frontend/test/<proposeOrdering>.test.ts packages/HimalayaUI/test/<picker>.jl
git commit -m "edit_tag: counted /api/sample-tags (modal ranking); proposeOrdering unchanged

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Slice 4 — `TagSuggest` + `ManageTagsModal`

### Task 7: `TagSuggest` combobox primitive

**Files:**
- Create: `frontend/src/print/ui/TagSuggest.tsx`, `frontend/src/print/ui/TagSuggest.stories.tsx`
- Modify: `frontend/src/print/ui/index.ts` (export)
- Test: `frontend/test/print-ui/TagSuggest.test.tsx`

- [ ] **Step 1: Write failing tests** (WAI-ARIA combobox semantics + counts + create-as-typed):

```tsx
import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { TagSuggest } from "../../src/print/ui/TagSuggest";

const opts = [{ label: "dose", count: 132 }, { label: "lipid_ratio", count: 88 }];

describe("<TagSuggest>", () => {
  it("is a combobox; opening shows a listbox of options with counts", () => {
    render(<TagSuggest mode="key" value="" options={opts} onCommit={() => {}} label="Key" />);
    const input = screen.getByRole("combobox", { name: "Key" });
    fireEvent.focus(input);
    expect(screen.getByRole("listbox")).toBeInTheDocument();
    expect(screen.getByRole("option", { name: /dose/ })).toHaveTextContent("132");
  });

  it("offers a create-as-typed option for a novel value and commits it", () => {
    const onCommit = vi.fn();
    render(<TagSuggest mode="value" value="" options={opts} onCommit={onCommit} label="Value" />);
    const input = screen.getByRole("combobox", { name: "Value" });
    fireEvent.change(input, { target: { value: "newval" } });
    fireEvent.click(screen.getByRole("option", { name: /use .newval./i }));
    expect(onCommit).toHaveBeenCalledWith("newval");
  });

  it("ArrowDown sets aria-activedescendant; Enter commits the active option", () => {
    const onCommit = vi.fn();
    render(<TagSuggest mode="key" value="" options={opts} onCommit={onCommit} label="Key" />);
    const input = screen.getByRole("combobox", { name: "Key" });
    fireEvent.focus(input);
    fireEvent.keyDown(input, { key: "ArrowDown" });
    expect(input).toHaveAttribute("aria-activedescendant");
    fireEvent.keyDown(input, { key: "Enter" });
    expect(onCommit).toHaveBeenCalledWith("dose");
  });
});
```

- [ ] **Step 2: Run, verify fail.**

- [ ] **Step 3: Implement `TagSuggest.tsx`** — a controlled combobox. Input `role="combobox"` with `aria-expanded`, `aria-controls`, `aria-autocomplete="list"`, `aria-activedescendant`; a `role="listbox"` of `role="option"` rows (filtered by the typed text, each showing `label` + faint `count`, low-count de-emphasized) plus a trailing create-as-typed option when the text matches no option. Down/Up move active index, Enter commits the active option (or the typed text), Escape collapses. Appearance (panel, option highlight, count styling) uses the house token classes that `Menu`/`Popover` already use — closed-look, only placement `className` exposed. Export from `ui/index.ts`. (Write the full component; ~120 lines. Reference `src/print/ui/Menu.tsx` for the listbox option styling vocabulary.)

- [ ] **Step 4: Run, verify pass; design-guard.**

Run: `npx vitest run test/print-ui/TagSuggest.test.tsx` → PASS; `npm run lint:design` → clean.

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/print/ui/TagSuggest.tsx packages/HimalayaUI/frontend/src/print/ui/TagSuggest.stories.tsx packages/HimalayaUI/frontend/src/print/ui/index.ts packages/HimalayaUI/frontend/test/print-ui/TagSuggest.test.tsx
git commit -m "TagSuggest: corpus-aware key/value combobox with counts (WAI-ARIA)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 8: `ManageTagsModal`

**Files:**
- Create: `frontend/src/print/components/ManageTagsModal.tsx`
- Test: `frontend/test/print-components/ManageTagsModal.test.tsx`

- [ ] **Step 1: Write failing tests** — rows render with ids; editing a value calls `onEdit(tagId, key, value)`; deleting calls `onRemove(tagId)`; adding calls `onAdd(key, value)`; adding a key the sample already has is rejected inline (no `onAdd`); the modal is `role=dialog` `aria-modal` labelled by its `h2`.

```tsx
it("editing a value commits the edit by tag id", () => {
  const onEdit = vi.fn();
  render(<ManageTagsModal open sampleName="LL37 Only"
    tags={[{ id: 3, key: "dose", value: "10", source: "manual" }]}
    keyOptions={[]} valueOptionsFor={() => []}
    onEdit={onEdit} onAdd={() => {}} onRemove={() => {}} onClose={() => {}} />);
  // open the value combobox for the dose row, pick/type 12
  // ... commit ...
  expect(onEdit).toHaveBeenCalledWith(3, "dose", "12");
});

it("rejects adding a key the sample already has (single-valued rule)", () => {
  const onAdd = vi.fn();
  render(<ManageTagsModal open sampleName="x"
    tags={[{ id: 3, key: "dose", value: "10", source: "manual" }]}
    keyOptions={[]} valueOptionsFor={() => []}
    onEdit={() => {}} onAdd={onAdd} onRemove={() => {}} onClose={() => {}} />);
  // type "dose" into the add-row key + a value, attempt commit
  // ... expect an alert and NO onAdd ...
  expect(onAdd).not.toHaveBeenCalled();
});
```

- [ ] **Step 2: Run, verify fail.**

- [ ] **Step 3: Implement `ManageTagsModal.tsx`** on `ModalShell size="md"` + `ModalHead` (`Kicker "Sample tags"` / serif `h2` = `sampleName`) + body of rows (`[key TagSuggest][value TagSuggest][IconButton dismiss aria-label="Remove ${key} ${value}"]`) + an add row + `ModalFooter` (note "Tags apply to the whole sample. Saved as you go." + `Done`). Props: `open`, `sampleName`, `tags: LoupeTag[]`, `keyOptions: {label,count}[]`, `valueOptionsFor: (key) => {label,count}[]`, `onEdit(tagId,key,value)`, `onAdd(key,value)`, `onRemove(tagId)`, `onClose`. Add-row + key-edit reject a key already in `tags` (inline `role=alert`, cleared on edit). Capture the opener ref via a `triggerRef` prop and restore focus on close. Use `lib/announce` for add/edit/remove. (Full component; reuse `CustomIndexModal.tsx` as the structural template.)

- [ ] **Step 4: Run, verify pass; design-guard; commit.**

```bash
git add packages/HimalayaUI/frontend/src/print/components/ManageTagsModal.tsx packages/HimalayaUI/frontend/test/print-components/ManageTagsModal.test.tsx
git commit -m "ManageTagsModal: per-sample tag editor on the custom-index shell

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 9: wire the modal into the loupe + inline duplicate-key rejection

**Files:**
- Modify: `frontend/src/print/components/LoupeSidePanel.tsx` + `LoupePage.tsx` (the "Manage" affordance, modal mount, `useEditCorpusSampleTag` + `useCorpusSampleTags` wiring)
- Modify: `frontend/src/print/ui/TagEditor.tsx` + `TagList.tsx` (thread the sample's existing keys; reject a duplicate-key inline `+ tag`)
- Test: extend the loupe page test + a `TagEditor` test
- E2E: `frontend/e2e/loupe-tags.spec.ts`

- [ ] **Step 1: Failing tests** — (a) `TagEditor` with an `existingKeys={["dose"]}` prop rejects committing key "dose" (no `onAdd`, shows the alert); (b) the loupe renders a "Manage" button that opens the modal; (c) an E2E that opens Manage, edits a value, and deletes the redundant duplicate row (assert rendered semantics, not classes).

- [ ] **Step 2: Run, verify fail.**

- [ ] **Step 3: Implement** — add `existingKeys?: string[]` to `TagEditor`/`TagList`; in `TagEditor.commit()`, after the empty-key check, reject when `existingKeys?.includes(key)` (set the same `aria-invalid` + alert it already uses for empty keys). In `LoupeSidePanel`, render a "Manage" `Button`/`IconButton` by the "Sample tags" kicker; in `LoupePage`, mount `ManageTagsModal` with `tags={toLoupeTags(sample.tags)}`, `keyOptions`/`valueOptionsFor` derived from `useCorpusSampleTags()` (keys ranked by summed count; values filtered by key and ranked by count), `onEdit={(tagId,key,value)=>editTag.mutate({tagId,key,value})}`, `onAdd`/`onRemove` to the existing corpus add/remove hooks, and pass the sample's current keys to the inline `TagList`.

- [ ] **Step 4: Run unit + e2e; render-verify at :5182** (open Manage on sample 3, edit a value, see the count dropdown, delete the duplicate). 

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/frontend/src/print/components/LoupeSidePanel.tsx packages/HimalayaUI/frontend/src/print/pages/LoupePage.tsx packages/HimalayaUI/frontend/src/print/ui/TagEditor.tsx packages/HimalayaUI/frontend/src/print/ui/TagList.tsx packages/HimalayaUI/frontend/test/<touched> packages/HimalayaUI/frontend/e2e/loupe-tags.spec.ts
git commit -m "Loupe: Manage-tags affordance + modal wiring + inline duplicate-key rejection

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Slice 5 — scoping batch upsert (the risk locus, last)

### Task 10: batch upserts on (sample_id, key)

**Files:**
- Modify: `packages/HimalayaUI/src/routes_samples.jl` (the `@post "/api/samples/tags/batch"` body loop)
- Test: the batch-route test file

- [ ] **Step 1: Write failing tests**

```julia
@testset "batch upsert: existing key updates in place, new key inserts, unchanged is a no-op" begin
    # sample S already has dose=10 (id X). Batch {key:dose, tags:[{sample_id:S, value:12}]}.
    n_events0 = only(DBInterface.execute(db,"SELECT COUNT(*) c FROM user_actions") |> Tables.rowtable).c
    HTTP.request("POST", "$base/api/samples/tags/batch", ct,
        JSON3.write(Dict(:key=>"dose", :source=>"scoping", :tags=>[Dict(:sample_id=>S, :value=>"12")])))
    rows = DBInterface.execute(db,"SELECT id,value FROM sample_tags WHERE sample_id=? AND key='dose'", [S]) |> Tables.rowtable
    @test length(rows) == 1 && rows[1].id == X && rows[1].value == "12"   # updated in place, no twin
end

@testset "batch upsert: unchanged value writes nothing and emits no event" begin
    # sample S has dose=12 (from above). Re-scope with the SAME value.
    n0 = only(DBInterface.execute(db,"SELECT COUNT(*) c FROM user_actions") |> Tables.rowtable).c
    HTTP.request("POST", "$base/api/samples/tags/batch", ct,
        JSON3.write(Dict(:key=>"dose", :source=>"scoping", :tags=>[Dict(:sample_id=>S, :value=>"12")])))
    n1 = only(DBInterface.execute(db,"SELECT COUNT(*) c FROM user_actions") |> Tables.rowtable).c
    @test n1 == n0   # idempotent: no new row, no new event
end
```

- [ ] **Step 2: Run, verify fail** (today it inserts a twin).

- [ ] **Step 3: Implement** — inside the `for t in entries` loop, before the INSERT, SELECT the existing row and branch:

```julia
                existing = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, value FROM sample_tags
                     WHERE sample_id = ? AND key = ? ORDER BY id LIMIT 1", [sample_id, key]))
                if isempty(existing)
                    res = DBInterface.execute(db,
                        "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?,?,?,?)",
                        [sample_id, key, value, source])
                    tag_id = Int(DBInterface.lastrowid(res))
                    result = apply_event!(InTransaction(), db, req; kind="add_tag",
                        entity_type="sample", entity_id=sample_id,
                        payload=Dict(:key=>key, :value=>value, :tag_id=>tag_id, :experiment_id=>exp_id))
                    _enqueue_broadcast_from_result!(result, "add_tag", "sample", sample_id)
                elseif String(existing[1].value) != value
                    tag_id = Int(existing[1].id)
                    DBInterface.execute(db,
                        "UPDATE sample_tags SET value = ? WHERE id = ?", [value, tag_id])
                    result = apply_event!(InTransaction(), db, req; kind="edit_tag",
                        entity_type="sample", entity_id=sample_id,
                        payload=Dict(:tag_id=>tag_id, :key=>key, :value=>value, :experiment_id=>exp_id))
                    _enqueue_broadcast_from_result!(result, "edit_tag", "sample", sample_id)
                    tag_id  # for the `created` push
                else
                    tag_id = Int(existing[1].id)   # unchanged: no write, no event
                end
```

(Adjust the `created` push to use the resolved `tag_id`. The pre-existing duplicate-`sample_id`-in-batch guard stays — it still prevents two entries racing on the same `(sample_id, key)` SELECT.)

- [ ] **Step 4: Run, verify pass.** Confirm the existing scoping E2E (`series-scoping.spec.ts`) still passes (no twin now).

- [ ] **Step 5: Commit.**

```bash
git add packages/HimalayaUI/src/routes_samples.jl packages/HimalayaUI/test/<batch-test>.jl
git commit -m "Scoping: batch upserts on (sample_id,key) — no more duplicate tags

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Final verification (after all slices)

- [ ] Full gate green: `lint:design` + `tsc(build)` + full `vitest` + `e2e` + `build`; Julia suite green.
- [ ] Render-verify at :5182: open the loupe on the duplicate sample, Manage → edit a value (count dropdown scoped to the key) → delete the redundant twin → re-scope a sample and confirm no new twin is born.
- [ ] Mark **LO-TAGDUP done** in `docs/superpowers/notes/loop-backlog.md` with the slice SHAs; re-score the Loupe surface.

## Self-review notes (writing-plans)

- **Spec coverage:** every spec section maps to a task — PATCH route (T3), edit_tag event (T3), upsert (T10), counted endpoint (T6), TagSuggest (T7), ManageTagsModal (T8), loupe integration + inline reject (T9), id-exact delete (T1-2), the four queue files (T4-5), a11y (T7-8). Non-goals (Layer 2, value typing, did-you-mean) are excluded.
- **Single-valued enforcement** appears in all write paths: modal/inline (T8/T9), PATCH 409 (T3), scoping upsert (T10).
- **Build order** is smallest-correct-first; the risk locus (upsert) is last behind the most tests.
- **Open implementation detail for the executor:** the exact test-file paths in `packages/HimalayaUI/test/` and a few existing frontend test files are named generically (grep for the nearest existing tag-route / loupe test and extend it rather than create a parallel file).
