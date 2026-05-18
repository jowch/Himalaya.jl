# Corpus Sample-Tag Mutations Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make sample-tag add/remove work from the corpus contact sheet by wiring a corpus-scoped path through all six layers of the mutation queue.

**Architecture:** Two new dedicated mutators (`addCorpusSampleTagMutator`, `removeCorpusSampleTagMutator`) patch the corpus query key `queryKeys.corpusSamples`. `resolveMutator` becomes a true tri-way split (exposure / experiment-sample / corpus-sample). `applyRemoteToCache` invalidates the corpus key on any sample-tag SSE frame. The single-tag route gains an optional `source` parameter. `resolveMutatorForEvent` stays deliberately 2-arm — verified in the design spec.

**Tech Stack:** Julia (Oxygen.jl REST, SQLite, stdlib `Test`); TypeScript/React (TanStack Query, Zustand, Vitest).

**Spec:** `docs/superpowers/specs/2026-05-18-corpus-sample-tag-mutations-design.md`

---

## Background the engineer needs

- **Mutation queue.** Every user gesture that writes data goes through a "mutator": `onMutate` writes an optimistic cache patch and returns a `restore` rollback closure; `request` calls the HTTP API; `onSuccess` replaces the optimistic placeholder with the server's canonical row. Cross-tab sync is handled separately by `applyRemoteToCache` (foreign SSE events). See `docs/mutation-queue.md`.
- **The corpus query.** `queryKeys.corpusSamples` is the static key `["corpus","samples"]` (defined `queries.ts:103`). `useCorpusSamples()` fetches `GET /api/samples` into a `CorpusSample[]` cache. `CorpusSample extends Sample` with an extra `q_units: string` field (`api.ts:40`).
- **Six-layer contract testing.** A tag mutation touches route emit, SSE frame, `applyRemoteToCache`, the cache row shape, `onMutate`, and `onSuccess`. A fix at one layer needs a test at every layer the bug class can occur. See `docs/contract-testing.md`.
- **Why two new mutators, not editing the existing ones.** The existing `addSampleTagMutator` / `removeSampleTagMutator` are experiment-scoped (they patch `queryKeys.samples(experimentId)`). They stay untouched. The corpus mutators are separate so the existing path keeps its behaviour and its passing tests. See the spec's "Scope decision".
- **The tri-scope discriminator order matters.** A corpus-tag op carries `{sampleId}` only — no `experimentId`. An exposure-tag op carries `{sampleId, exposureId}` — also no `experimentId`. So `resolveMutator` must test `exposureId` **first**, then `experimentId`, then fall through to corpus. Getting the order wrong sends corpus tags to the exposure mutator.

## Running tests

**Frontend** — from `packages/HimalayaUI/frontend/`:
```bash
npm test -- test/queue/<file>.test.ts      # one-shot (the --run flag is added by a hook)
npx tsc --noEmit                           # typecheck only
npm run build                              # tsc --noEmit + vite build
```

**Backend (targeted single file)** — from the worktree root:
```bash
julia --project=packages/HimalayaUI -e '
using HimalayaUI, Test, HTTP, JSON3, SQLite, DBInterface, Tables
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_routes_samples.jl")
'
```
`test_http.jl` defines the `with_test_server` helper and must be `include`d first (it runs a small fast testset of its own — that is expected and harmless).

## File map

| File | Change |
|---|---|
| `packages/HimalayaUI/src/routes_samples.jl` | Modify `POST /api/samples/{id}/tags` — optional `source` param |
| `packages/HimalayaUI/test/test_routes_samples.jl` | New testset for `source` parameterization |
| `packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts` | Add `addCorpusSampleTagMutator` + `removeCorpusSampleTagMutator` |
| `packages/HimalayaUI/frontend/src/lib/queue/mutatorRegistry.ts` | Tri-scope `resolveMutator`; comment `resolveMutatorForEvent` |
| `packages/HimalayaUI/frontend/src/lib/queue/applyRemoteToCache.ts` | `add_tag`/`remove_tag` branch invalidates the corpus key |
| `packages/HimalayaUI/frontend/src/queries.ts` | Add `useAddCorpusSampleTag` / `useRemoveCorpusSampleTag` hooks |
| `packages/HimalayaUI/frontend/test/queue/cache-shape.test.ts` | Corpus `onSuccess` row-shape test |
| `packages/HimalayaUI/frontend/test/queue/rollbackSymmetry.test.ts` | Corpus `onMutate` ↔ rollback test |
| `packages/HimalayaUI/frontend/test/queue/authHeaders.test.ts` | Corpus mutator header-propagation test |
| `packages/HimalayaUI/frontend/test/queue/mutatorRegistry.test.ts` | Tri-scope resolver test |
| `packages/HimalayaUI/frontend/test/queue/sseEventPayload.contract.test.ts` | Corpus-key invalidation test |
| `packages/HimalayaUI/frontend/test/queue/mutatorOnSseWins.test.ts` | Cross-mutator SSE-wins handoff test |

`test_route_response_shapes.jl` needs **no change** — it already pins `POST /api/samples/:id/tags → {id,key,value,source}` (`:448`); the key set is unchanged by this work.

---

## Task 1: Route — parameterize tag `source`

**Files:**
- Modify: `packages/HimalayaUI/src/routes_samples.jl:123-160` (the `POST /api/samples/{id}/tags` handler)
- Test: `packages/HimalayaUI/test/test_routes_samples.jl`

- [ ] **Step 1: Write the failing test**

Append this `@testset` to `packages/HimalayaUI/test/test_routes_samples.jl` (after the existing `"POST /api/samples/:id/tags is idempotent under retry"` testset):

```julia
@testset "POST /api/samples/:id/tags — source defaults to manual, accepts an explicit value, rejects non-string" begin
    db = SQLite.DB()
    HimalayaUI.create_schema!(db)
    exp_id = HimalayaUI.init_experiment!(db; path="/tsrc", data_dir="/tsrc/d",
                                             analysis_dir="/tsrc/a")
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1", display_name="UX1")

    with_test_server(db) do port, base
        hdrs = ["Content-Type" => "application/json", "X-Username" => "alice"]

        # No source field → defaults to 'manual'.
        r = HTTP.post("$base/api/samples/$s_id/tags";
            body = JSON3.write(Dict(:key => "lipid", :value => "DOPC")),
            headers = hdrs)
        @test r.status == 201
        @test JSON3.read(String(r.body)).source == "manual"

        # Explicit source is honored end to end.
        r = HTTP.post("$base/api/samples/$s_id/tags";
            body = JSON3.write(Dict(:key => "temp", :value => "25C",
                                    :source => "scoping")),
            headers = hdrs)
        @test r.status == 201
        @test JSON3.read(String(r.body)).source == "scoping"

        # Non-string source → 400, nothing written.
        r = HTTP.post("$base/api/samples/$s_id/tags";
            body = JSON3.write(Dict(:key => "x", :value => "y", :source => 42)),
            headers = hdrs)
        @test r.status == 400

        # Exactly the two successful inserts persisted, with the right source.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT key, source FROM sample_tags WHERE sample_id = ? ORDER BY id",
            [s_id]))
        @test length(rows) == 2
        @test rows[1].source == "manual"
        @test rows[2].source == "scoping"
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
julia --project=packages/HimalayaUI -e '
using HimalayaUI, Test, HTTP, JSON3, SQLite, DBInterface, Tables
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_routes_samples.jl")
'
```
Expected: FAIL — the explicit-source assertion gets `"manual"` (route hard-codes it) and the non-string-source POST returns `201` instead of `400`.

- [ ] **Step 3: Modify the route handler**

In `packages/HimalayaUI/src/routes_samples.jl`, replace the body of the `@post "/api/samples/{id}/tags"` handler. The current handler reads:

```julia
    @post "/api/samples/{id}/tags" function(req::HTTP.Request, id::Int)
        db    = current_db()
        body  = json(req)
        if !haskey(body, :key) || !haskey(body, :value)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "missing required fields: key, value")))
        end
        key   = String(body.key)
        value = String(body.value)
        return with_idempotency(db, req) do
            res   = DBInterface.execute(db,
                "INSERT INTO sample_tags (sample_id, key, value, source)
                 VALUES (?, ?, ?, 'manual')",
                [id, key, value])
```

Replace exactly those lines with:

```julia
    @post "/api/samples/{id}/tags" function(req::HTTP.Request, id::Int)
        db    = current_db()
        body  = json(req)
        if !haskey(body, :key) || !haskey(body, :value)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "missing required fields: key, value")))
        end
        # `source` is optional and defaults to 'manual'. Validated BEFORE the
        # transaction opens so a malformed request is a clean 400 (mirrors the
        # batch route's validation discipline). The series-scoping path passes
        # an explicit non-'manual' source.
        if haskey(body, :source) && !(body.source isa AbstractString)
            return HTTP.Response(400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict(:error => "source must be a string")))
        end
        key    = String(body.key)
        value  = String(body.value)
        source = haskey(body, :source) ? String(body.source) : "manual"
        return with_idempotency(db, req) do
            res   = DBInterface.execute(db,
                "INSERT INTO sample_tags (sample_id, key, value, source)
                 VALUES (?, ?, ?, ?)",
                [id, key, value, source])
```

Then, further down in the same handler, change the 201 response — the current line:

```julia
            HTTP.Response(201, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:id => tag_id,
                                 :key => key, :value => value, :source => "manual")))
```

becomes:

```julia
            HTTP.Response(201, ["Content-Type" => "application/json"],
                JSON3.write(Dict(:id => tag_id,
                                 :key => key, :value => value, :source => source)))
```

Leave the `apply_event!` call and its payload (`:key, :value, :tag_id, :experiment_id`) **unchanged** — the event payload deliberately does not carry `source` (see spec "Out of scope").

- [ ] **Step 4: Run the test to verify it passes**

```bash
julia --project=packages/HimalayaUI -e '
using HimalayaUI, Test, HTTP, JSON3, SQLite, DBInterface, Tables
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_routes_samples.jl")
'
```
Expected: PASS — all `samples routes` testsets green, including the new one.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/src/routes_samples.jl packages/HimalayaUI/test/test_routes_samples.jl
git commit -m "$(cat <<'EOF'
Parameterize tag source on the single sample-tag route (#159)

POST /api/samples/{id}/tags now accepts an optional `source` body
field, defaulting to 'manual', validated string-if-present. Brings the
single-tag route to parity with the batch route. The add_tag event
payload is unchanged.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Corpus mutators + cache-shape test

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts`
- Test: `packages/HimalayaUI/frontend/test/queue/cache-shape.test.ts`

- [ ] **Step 1: Write the failing test**

In `packages/HimalayaUI/frontend/test/queue/cache-shape.test.ts`, add `addCorpusSampleTagMutator` to the import from `../../src/lib/queue/mutators/trivial` (the existing import block already pulls `addSampleTagMutator` from that module — add the new name to it).

Then add this test inside the `describe("Cache-shape integrity ...")` block, right after the existing `"addSampleTag inserts a SampleTag with exactly 4 keys ..."` test:

```ts
  it("addCorpusSampleTag inserts a SampleTag with exactly 4 keys into the corpus cache", async () => {
    // Corpus contact sheet path: the op carries no experimentId; the mutator
    // patches queryKeys.corpusSamples, not queryKeys.samples(experimentId).
    const initialSample = {
      id: 10, experiment_id: 1, display_name: "D1", name: "n", notes: null,
      q_units: "nm^-1", tags: [],
    };
    qc.setQueryData(queryKeys.corpusSamples, [initialSample]);
    mockFetchOnce({
      id: 50, sample_id: 10, key: "color", value: "red", source: "manual",
    }, 201);
    await runMutator(qc, addCorpusSampleTagMutator, {
      kind: "add_tag",
      clientOpId: "op-corpus-shape-1",
      sampleId: 10, username: "alice", clientId: "tab-1",
      key: "color", value: "red",
      payload: { key: "color", value: "red" },
    });
    const list = qc.getQueryData<{ tags: unknown[]; q_units: string }[]>(
      queryKeys.corpusSamples);
    const tag = list![0]!.tags[0];
    assertKeys(tag, SAMPLE_TAG_KEYS, "addCorpusSampleTag cache row");
    // The mutator mutates `tags` in place — q_units survives untouched.
    expect(list![0]!.q_units).toBe("nm^-1");
  });
```

- [ ] **Step 2: Run the test to verify it fails**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- test/queue/cache-shape.test.ts`
Expected: FAIL — `addCorpusSampleTagMutator` is not exported (import error / undefined).

- [ ] **Step 3: Add the corpus mutators**

In `packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts`:

First, add `CorpusSample` to the type import from `../../../api` (the existing import already lists `Sample, SampleTag, Exposure, ExposureTag, SampleMessage, ComparisonMessage, AuthOpts` — add `CorpusSample`).

Then add this section immediately after the `removeExposureTagMutator` definition (after the `// remove_tag (exposure)` block, before `// post_message`):

```ts
// ---------------------------------------------------------------------------
// add_tag / remove_tag (corpus sample) — #159
// ---------------------------------------------------------------------------
// The corpus contact sheet reads queryKeys.corpusSamples (a CorpusSample[]),
// not the per-experiment queryKeys.samples(experimentId). These mutators are
// separate from addSampleTagMutator/removeSampleTagMutator so the existing
// experiment-scoped path is untouched. The op scope carries `sampleId` only —
// no `experimentId` — which is how resolveMutator routes here. They define no
// `synthesizeFromSse`: the own-op SSE-wins response shape is produced by
// addSampleTagMutator.synthesizeFromSse via resolveMutatorForEvent (shared
// SampleTag shape) — see the comment in mutatorRegistry.ts.

export type AddCorpusSampleTagInput = { key: string; value: string };
type AddCorpusSampleTagScope = BaseScope & { sampleId: number };

export const addCorpusSampleTagMutator: Mutator<
  AddCorpusSampleTagInput, AddCorpusSampleTagScope, SampleTag
> = {
  kind: "add_tag",
  onMutate: (p, qc): RollbackContext => {
    const key = queryKeys.corpusSamples;
    const prev = qc.getQueryData<CorpusSample[]>(key);
    const placeholderId = nextOptimisticId();
    if (prev) {
      qc.setQueryData<CorpusSample[]>(key, prev.map((s) =>
        s.id === p.sampleId
          ? { ...s, tags: [...s.tags, {
              id: placeholderId, key: p.key, value: p.value, source: "manual",
            }] }
          : s,
      ));
    }
    return {
      restore: () => {
        if (prev !== undefined) qc.setQueryData(key, prev);
      },
    };
  },
  request: (p) => api.addSampleTag(p.sampleId, p.key, p.value, buildAuthOpts(p)),
  onSuccess: (p, response, qc) => {
    // Route emits {id, sample_id, key, value, source}; SampleTag omits
    // sample_id. Strip it so the cached row matches the type.
    const tag: SampleTag = {
      id: response.id, key: response.key, value: response.value,
      source: response.source,
    };
    const key = queryKeys.corpusSamples;
    qc.setQueryData<CorpusSample[]>(key, (list) => {
      if (!list) return list;
      return list.map((s) =>
        s.id !== p.sampleId
          ? s
          : {
              ...s,
              tags: replacePlaceholder(
                s.tags,
                tag,
                (t) => t.key === p.key && t.value === p.value,
              ),
            },
      );
    });
  },
};

export type RemoveCorpusSampleTagInput = { tagId: number };
type RemoveCorpusSampleTagScope = BaseScope & { sampleId: number };

export const removeCorpusSampleTagMutator: Mutator<
  RemoveCorpusSampleTagInput, RemoveCorpusSampleTagScope, void
> = {
  kind: "remove_tag",
  onMutate: (p, qc): RollbackContext => {
    const key = queryKeys.corpusSamples;
    const prev = qc.getQueryData<CorpusSample[]>(key);
    if (prev) {
      qc.setQueryData<CorpusSample[]>(key, prev.map((s) =>
        s.id === p.sampleId
          ? { ...s, tags: s.tags.filter((t) => t.id !== p.tagId) }
          : s,
      ));
    }
    return {
      restore: () => {
        if (prev !== undefined) qc.setQueryData(key, prev);
      },
    };
  },
  request: (p) => api.removeSampleTag(p.sampleId, p.tagId, buildAuthOpts(p)),
  onSuccess: () => {},
  // 404 = the tag is already gone → desired end state.
  treats404AsSuccess: true,
};
```

- [ ] **Step 4: Run the test to verify it passes**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- test/queue/cache-shape.test.ts`
Expected: PASS — all cache-shape tests green, including `addCorpusSampleTag`.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts \
        packages/HimalayaUI/frontend/test/queue/cache-shape.test.ts
git commit -m "$(cat <<'EOF'
Add corpus sample-tag mutators (#159)

addCorpusSampleTagMutator / removeCorpusSampleTagMutator patch the
corpus query key (queryKeys.corpusSamples). Scope carries sampleId
only — no experimentId — keeping them distinct from the untouched
experiment-scoped mutators.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Rollback-symmetry test for the corpus mutators

**Files:**
- Test: `packages/HimalayaUI/frontend/test/queue/rollbackSymmetry.test.ts`

- [ ] **Step 1: Write the failing test**

In `packages/HimalayaUI/frontend/test/queue/rollbackSymmetry.test.ts`, add `addCorpusSampleTagMutator, removeCorpusSampleTagMutator` to the import from `../../src/lib/queue/mutators/trivial` (the block that already imports `addSampleTagMutator, removeSampleTagMutator`).

Then add these two entries to the `SPECS` array, right after the existing `"removeSampleTag"` spec:

```ts
  {
    name: "addCorpusSampleTag",
    keys: [queryKeys.corpusSamples],
    seed: (qc) => qc.setQueryData(queryKeys.corpusSamples, [SAMPLE]),
    run: (qc) => {
      const ctx = addCorpusSampleTagMutator.onMutate({
        kind: "add_tag", clientOpId: "op",
        payload: { key: "x", value: "y" },
        sampleId: 10, username: "alice", clientId: "tab",
        key: "x", value: "y",
      } as any, qc);
      ctx.restore();
    },
  },
  {
    name: "removeCorpusSampleTag",
    keys: [queryKeys.corpusSamples],
    seed: (qc) => qc.setQueryData(queryKeys.corpusSamples, [SAMPLE]),
    run: (qc) => {
      const ctx = removeCorpusSampleTagMutator.onMutate({
        kind: "remove_tag", clientOpId: "op",
        payload: { tagId: 1 },
        sampleId: 10, username: "alice", clientId: "tab",
        tagId: 1,
      } as any, qc);
      ctx.restore();
    },
  },
```

(The existing `SAMPLE` fixture — `{ id: 10, ..., tags: [{ id: 1, key: "k", value: "v", source: "manual" }] }` — is reused; `removeCorpusSampleTag` with `tagId: 1` removes that seeded tag, and `restore()` must put it back.)

- [ ] **Step 2: Run the test to verify it fails**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- test/queue/rollbackSymmetry.test.ts`
Expected: FAIL — `addCorpusSampleTagMutator` / `removeCorpusSampleTagMutator` import fails (this task runs before Task 2 only if reordered; if Task 2 is already done the test instead fails only if the mutators' `restore` is wrong — it should be correct, so confirm a real run-then-pass cycle by temporarily breaking nothing). If Task 2 is complete, this is a test-only addition that should pass immediately; treat Step 2 as "run, observe, the two new specs are present and green".

- [ ] **Step 3: (no implementation needed)**

This is a test-only task — the mutators from Task 2 already satisfy it. If Step 2 shows a failure, it is a real rollback bug in the Task 2 mutators; fix the `restore` closure there.

- [ ] **Step 4: Run the test to verify it passes**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- test/queue/rollbackSymmetry.test.ts`
Expected: PASS — `addCorpusSampleTag` and `removeCorpusSampleTag` specs green.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/test/queue/rollbackSymmetry.test.ts
git commit -m "$(cat <<'EOF'
Pin rollback symmetry for the corpus sample-tag mutators (#159)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Auth-header propagation test for the corpus mutators

**Files:**
- Test: `packages/HimalayaUI/frontend/test/queue/authHeaders.test.ts`

- [ ] **Step 1: Write the failing test**

In `packages/HimalayaUI/frontend/test/queue/authHeaders.test.ts`, add `addCorpusSampleTagMutator, removeCorpusSampleTagMutator` to the import from `../../src/lib/queue/mutators/trivial`.

Then add these two entries to the `SPECS` array, after the existing `"removeExposureTag"` spec:

```ts
  {
    name: "addCorpusSampleTag",
    run: (qc) => addCorpusSampleTagMutator.request(
      { ...FLAT_BASE, kind: "add_tag",
        payload: { key: "k", value: "v" },
        sampleId: 10, key: "k", value: "v" } as any,
      new AbortController().signal),
  },
  {
    name: "removeCorpusSampleTag",
    run: (qc) => removeCorpusSampleTagMutator.request(
      { ...FLAT_BASE, kind: "remove_tag",
        payload: { tagId: 1 },
        sampleId: 10, tagId: 1 } as any,
      new AbortController().signal),
  },
```

- [ ] **Step 2: Run the test to verify it fails**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- test/queue/authHeaders.test.ts`
Expected: FAIL before Task 2 (import error). With Task 2 complete this is a test-only addition — run it and confirm the two new specs appear and pass.

- [ ] **Step 3: (no implementation needed)**

Test-only task — the corpus mutators reuse `buildAuthOpts(p)` exactly as the experiment mutators do, so the headers propagate. If Step 2 fails on a header assertion, the bug is a missing `buildAuthOpts` in the Task 2 mutator's `request`; fix it there.

- [ ] **Step 4: Run the test to verify it passes**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- test/queue/authHeaders.test.ts`
Expected: PASS — `addCorpusSampleTag` and `removeCorpusSampleTag` send all three headers.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/test/queue/authHeaders.test.ts
git commit -m "$(cat <<'EOF'
Pin auth-header propagation for the corpus sample-tag mutators (#159)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Corpus tag hooks in `queries.ts`

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/queries.ts`

- [ ] **Step 1: Add the mutator imports**

In `packages/HimalayaUI/frontend/src/queries.ts`, the import block ending at line 19 (`} from "./lib/queue/mutators/trivial";`) already imports the trivial mutators. Add `addCorpusSampleTagMutator` and `removeCorpusSampleTagMutator` to that import list.

- [ ] **Step 2: Add the two hooks**

In `queries.ts`, immediately after the existing `useRemoveSampleTag` function, add:

```ts
/**
 * Corpus contact-sheet sample tagging (#159). Scope carries `sampleId` only —
 * NO `experimentId` — which is how resolveMutator routes the op to the
 * corpus mutator rather than the experiment-scoped or exposure mutator.
 */
export function useAddCorpusSampleTag(sampleId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(
    addCorpusSampleTagMutator,
    { sampleId, username, clientId: CLIENT_ID },
  );
}

export function useRemoveCorpusSampleTag(sampleId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    removeCorpusSampleTagMutator,
    { sampleId, username, clientId: CLIENT_ID },
  );
  return {
    ...inner,
    mutate: (tagId: number) => inner.mutate({ tagId }),
  };
}
```

- [ ] **Step 3: Typecheck**

Run (from `packages/HimalayaUI/frontend/`): `npx tsc --noEmit`
Expected: no errors. (These hooks have no dedicated unit test — they are thin wrappers; `tsc` is the contract here, and the contact-sheet UI in #160 will exercise them.)

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/src/queries.ts
git commit -m "$(cat <<'EOF'
Add useAddCorpusSampleTag / useRemoveCorpusSampleTag hooks (#159)

Scope carries sampleId only — no experimentId — so resolveMutator
routes the op to the corpus mutator.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Tri-scope `resolveMutator` + `resolveMutatorForEvent` comment

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutatorRegistry.ts`
- Test: `packages/HimalayaUI/frontend/test/queue/mutatorRegistry.test.ts`

- [ ] **Step 1: Write the failing test**

In `packages/HimalayaUI/frontend/test/queue/mutatorRegistry.test.ts`, add `addCorpusSampleTagMutator, removeCorpusSampleTagMutator` to the import from `../../src/lib/queue/mutators/trivial`.

Add this test inside `describe("resolveMutator", ...)`, after the existing `"dispatches remove_tag the same way"` test:

```ts
  it("dispatches a corpus sample-tag op (no experimentId, no exposureId) to the corpus mutator", () => {
    expect(
      resolveMutator({
        kind: "add_tag",
        payload: { sampleId: 10, key: "k", value: "v" },
      }),
    ).toBe(addCorpusSampleTagMutator);
    expect(
      resolveMutator({
        kind: "remove_tag",
        payload: { sampleId: 10, tagId: 7 },
      }),
    ).toBe(removeCorpusSampleTagMutator);
  });
```

Then **replace** the existing `"tolerates undefined/null payload by treating add_tag as exposure-scoped"` test with this (the tri-scope fallthrough changes the degenerate-payload target from exposure to corpus):

```ts
  it("tolerates undefined/null payload by treating add_tag as corpus-scoped", () => {
    // No exposureId and no experimentId → corpus fallthrough. Must not throw.
    expect(
      resolveMutator({ kind: "add_tag", payload: undefined }),
    ).toBe(addCorpusSampleTagMutator);
    expect(resolveMutator({ kind: "add_tag", payload: null })).toBe(
      addCorpusSampleTagMutator,
    );
  });
```

> **Do NOT** add the corpus mutators to the `describe("resolveMutator ↔ resolveMutatorForEvent consistency")` `cases` array at the bottom of this file. That array round-trips every mutator through `resolveMutatorForEvent`, which deliberately stays 2-arm and **never** returns a corpus mutator (it returns `addSampleTagMutator` for `add_tag`/`sample`). A corpus row there would assert `resolveMutatorForEvent("add_tag","sample") === addCorpusSampleTagMutator` and fail. The corpus mutators are intentionally absent from that cross-check — leave the array as it is.

- [ ] **Step 2: Run the test to verify it fails**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- test/queue/mutatorRegistry.test.ts`
Expected: FAIL — a corpus op currently routes to `addExposureTagMutator` (no `experimentId` → exposure branch), so both the new test and the rewritten degenerate-payload test fail.

- [ ] **Step 3: Make `resolveMutator` tri-scope**

In `packages/HimalayaUI/frontend/src/lib/queue/mutatorRegistry.ts`:

Add the two corpus mutators to the import from `./mutators/trivial` (the block that imports `addSampleTagMutator`, `removeSampleTagMutator`, etc.).

Replace the `add_tag` and `remove_tag` cases in `resolveMutator` — current:

```ts
    case "add_tag":
      return p?.experimentId !== undefined
        ? addSampleTagMutator
        : addExposureTagMutator;
    case "remove_tag":
      return p?.experimentId !== undefined
        ? removeSampleTagMutator
        : removeExposureTagMutator;
```

with:

```ts
    case "add_tag":
      // Tri-scope (#159). Peel exposure off FIRST: an exposure-tag op carries
      // {sampleId, exposureId} and — like a corpus-tag op — has no
      // experimentId, so the two collide unless exposureId is tested first.
      if (p?.exposureId !== undefined) return addExposureTagMutator;
      return p?.experimentId !== undefined
        ? addSampleTagMutator
        : addCorpusSampleTagMutator;
    case "remove_tag":
      if (p?.exposureId !== undefined) return removeExposureTagMutator;
      return p?.experimentId !== undefined
        ? removeSampleTagMutator
        : removeCorpusSampleTagMutator;
```

(The local payload cast at the top of `resolveMutator` already declares `exposureId?: number` — no type change is needed.)

- [ ] **Step 4: Add the `resolveMutatorForEvent` comment**

In the same file, in `resolveMutatorForEvent`, the `add_tag` / `remove_tag` cases currently read:

```ts
    case "add_tag":
      return entityType === "sample" ? addSampleTagMutator : addExposureTagMutator;
    case "remove_tag":
      return entityType === "sample" ? removeSampleTagMutator : removeExposureTagMutator;
```

Replace with the same logic plus a documenting comment:

```ts
    // add_tag / remove_tag stay 2-arm here — deliberately NOT tri-scope like
    // resolveMutator above. Verified (#159): a corpus-originated and an
    // experiment-originated sample tag are byte-identical on the SSE wire
    // (same entity_type="sample", same payload — the route always resolves
    // experiment_id from the sample row). There is no corpus-vs-experiment
    // discriminator to key on, and none is needed: this function only picks a
    // synthesizeFromSse for the own-op SSE-confirmation *response shape*.
    // addSampleTagMutator.synthesizeFromSse yields a SampleTag, which the
    // corpus mutator's own onSuccess consumes; the cache patch is done by the
    // pending mutation's own onSuccess (the corpus mutator). A literal third
    // arm would be unreachable dead code. This is a conscious, plan-aware
    // deviation from the literal "convert both to tri-scope" wording of #159
    // / master plan §11. See
    // docs/superpowers/specs/2026-05-18-corpus-sample-tag-mutations-design.md.
    case "add_tag":
      return entityType === "sample" ? addSampleTagMutator : addExposureTagMutator;
    case "remove_tag":
      return entityType === "sample" ? removeSampleTagMutator : removeExposureTagMutator;
```

- [ ] **Step 5: Run the test to verify it passes**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- test/queue/mutatorRegistry.test.ts`
Expected: PASS — all `resolveMutator`, `resolveMutatorForEvent`, and consistency-cross-check tests green. (The existing `add_tag (experimentId)` → `addSampleTagMutator` and `add_tag (exposureId)` → `addExposureTagMutator` assertions still hold; the consistency block calls only `resolveMutatorForEvent`, which is unchanged in behaviour.)

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutatorRegistry.ts \
        packages/HimalayaUI/frontend/test/queue/mutatorRegistry.test.ts
git commit -m "$(cat <<'EOF'
Make resolveMutator tri-scope for corpus sample tags (#159)

resolveMutator now routes add_tag/remove_tag three ways: exposureId
present -> exposure mutator; experimentId present -> experiment-sample
mutator; neither -> corpus-sample mutator. exposureId is tested first
because corpus and exposure ops both lack experimentId.

resolveMutatorForEvent stays 2-arm with a documenting comment — a
corpus and an experiment sample tag are byte-identical on the SSE
wire, so a literal third arm would be dead code.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: `applyRemoteToCache` corpus-key invalidation

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/applyRemoteToCache.ts`
- Test: `packages/HimalayaUI/frontend/test/queue/sseEventPayload.contract.test.ts`

- [ ] **Step 1: Write the failing test**

In `packages/HimalayaUI/frontend/test/queue/sseEventPayload.contract.test.ts`, **replace** the existing `"add_tag (sample) invalidates the experiment's samples list (uses experiment_id)"` test with:

```ts
  it("add_tag (sample) invalidates BOTH the experiment samples list and the corpus list", () => {
    const invalidated: unknown[] = [];
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidated.push(arg.queryKey); return Promise.resolve();
    }) as typeof qc.invalidateQueries;
    // Mirrors routes_samples.jl POST: payload always includes experiment_id.
    // A sample tag lives in two cached projections — the per-experiment list
    // and the corpus contact-sheet list — so both keys must invalidate.
    const evt: SseEvent = {
      id: 99, kind: "add_tag", entity_type: "sample", entity_id: 10,
      payload: { key: "k", value: "v", tag_id: 50, experiment_id: 1 },
    };
    applyRemoteToCache(evt, qc);
    expect(invalidated).toEqual([queryKeys.samples(1), queryKeys.corpusSamples]);
  });
```

And **replace** the existing `"remove_tag (sample) invalidates the experiment's samples list"` test with:

```ts
  it("remove_tag (sample) invalidates BOTH the experiment samples list and the corpus list", () => {
    const invalidated: unknown[] = [];
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidated.push(arg.queryKey); return Promise.resolve();
    }) as typeof qc.invalidateQueries;
    const evt: SseEvent = {
      id: 99, kind: "remove_tag", entity_type: "sample", entity_id: 10,
      payload: { tag_id: 50, experiment_id: 1 },
    };
    applyRemoteToCache(evt, qc);
    expect(invalidated).toEqual([queryKeys.samples(1), queryKeys.corpusSamples]);
  });
```

Leave the two `add_tag (exposure)` / `remove_tag (exposure)` tests unchanged — the exposure branch still invalidates exactly one key.

- [ ] **Step 2: Run the test to verify it fails**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- test/queue/sseEventPayload.contract.test.ts`
Expected: FAIL — the current branch invalidates only `queryKeys.samples(1)`, so `invalidated` is `[["experiment",1,"samples"]]`, not the expected two-element array.

- [ ] **Step 3: Make the `add_tag`/`remove_tag` branch tri-scope**

In `packages/HimalayaUI/frontend/src/lib/queue/applyRemoteToCache.ts`, replace the `add_tag` / `remove_tag` case — current:

```ts
    case "add_tag":
    case "remove_tag": {
      const parentKey = remote.entity_type === "sample"
        ? queryKeys.samples(payload?.experiment_id as number)
        : queryKeys.exposures(payload?.sample_id as number);
      qc.invalidateQueries({ queryKey: parentKey });
      break;
    }
```

with:

```ts
    case "add_tag":
    case "remove_tag": {
      if (remote.entity_type === "sample") {
        // A sample tag appears in two cached projections: the per-experiment
        // samples list AND the corpus-wide contact-sheet list. The route
        // always emits experiment_id, so the experiment key still invalidates
        // (a harmless no-op if not cached); the corpusSamples invalidation is
        // what refreshes the contact sheet from a foreign tab (#159).
        qc.invalidateQueries({
          queryKey: queryKeys.samples(payload?.experiment_id as number),
        });
        qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
      } else {
        qc.invalidateQueries({
          queryKey: queryKeys.exposures(payload?.sample_id as number),
        });
      }
      break;
    }
```

(`queryKeys` is already imported in this file; `corpusSamples` is a property on it.)

- [ ] **Step 4: Run the test to verify it passes**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- test/queue/sseEventPayload.contract.test.ts`
Expected: PASS — sample tags invalidate both keys, exposure tags invalidate one.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/applyRemoteToCache.ts \
        packages/HimalayaUI/frontend/test/queue/sseEventPayload.contract.test.ts
git commit -m "$(cat <<'EOF'
Invalidate the corpus key on sample-tag SSE frames (#159)

applyRemoteToCache's add_tag/remove_tag branch now invalidates both
queryKeys.samples(experiment_id) and queryKeys.corpusSamples for
entity_type="sample" frames, so a foreign tab's tag refreshes the
corpus contact sheet. The exposure branch is unchanged.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: SSE-wins cross-mutator handoff test

**Files:**
- Test: `packages/HimalayaUI/frontend/test/queue/mutatorOnSseWins.test.ts`

This pins the path where this tab's own corpus `add_tag` is confirmed by an SSE frame before HTTP returns: `synthesizeResponseFromSse` resolves the synth via `resolveMutatorForEvent("add_tag","sample")` → `addSampleTagMutator.synthesizeFromSse`, and the pending mutation's own `onSuccess` (the **corpus** mutator) consumes it. No existing kind splits synth and `onSuccess` across two mutators — this is novel to #159.

- [ ] **Step 1: Write the failing test**

In `packages/HimalayaUI/frontend/test/queue/mutatorOnSseWins.test.ts`, add to the import from `../../src/lib/queue/mutators/trivial` (currently `import { updateSampleMutator } from ...`) so it reads:

```ts
import {
  updateSampleMutator,
  addSampleTagMutator,
  addCorpusSampleTagMutator,
} from "../../src/lib/queue/mutators/trivial";
```

Then add this test inside the `describe("mutator onSuccess on SSE-wins synthetic responses", ...)` block:

```ts
  it("corpus add_tag SSE-wins: addSampleTagMutator synth feeds the corpus onSuccess", () => {
    // Own corpus add_tag confirmed by SSE before HTTP returns.
    // synthesizeResponseFromSse picks the synth via
    // resolveMutatorForEvent("add_tag","sample") -> addSampleTagMutator,
    // but the pending mutation is the CORPUS mutator, so ITS onSuccess
    // consumes that synth. Pin this cross-mutator handoff.
    const sample = {
      id: 10, experiment_id: 1, name: "n", display_name: "D1", notes: null,
      q_units: "nm^-1",
      tags: [{ id: -1, key: "lipid", value: "DOPC", source: "manual" }],
    };
    qc.setQueryData(queryKeys.corpusSamples, [sample]);

    // The SSE frame the server broadcasts for this tag.
    const remote = {
      id: 77, kind: "add_tag", entity_type: "sample", entity_id: 10,
      client_op_id: "op-corpus-1",
      payload: { tag_id: 500, key: "lipid", value: "DOPC", experiment_id: 1 },
    } as any;
    const base = {
      event_id: 77, client_op_id: "op-corpus-1",
      analysis_inputs_hash: undefined,
    };
    // synthesizeResponseFromSse routes through addSampleTagMutator's synth.
    const synth = addSampleTagMutator.synthesizeFromSse!(remote, base as any);
    expect(synth).toMatchObject({
      id: 500, key: "lipid", value: "DOPC", source: "manual",
    });

    // The pending mutation is the corpus mutator — its onSuccess runs.
    addCorpusSampleTagMutator.onSuccess(
      { sampleId: 10, key: "lipid", value: "DOPC",
        username: "u", clientId: "c", clientOpId: "op-corpus-1" } as any,
      synth as any,
      qc,
    );

    const list = qc.getQueryData<any[]>(queryKeys.corpusSamples);
    // The optimistic placeholder (id: -1) is replaced by the canonical row.
    expect(list![0].tags).toEqual([
      { id: 500, key: "lipid", value: "DOPC", source: "manual" },
    ]);
    // q_units preserved — the mutator never reconstructs the CorpusSample row.
    expect(list![0].q_units).toBe("nm^-1");
  });
```

- [ ] **Step 2: Run the test to verify it fails first, then passes**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- test/queue/mutatorOnSseWins.test.ts`

This task has no implementation step — it pins behaviour delivered by Tasks 2 and 6. Expected: **PASS** if Tasks 2 and 6 are complete. If it FAILS, that is a real defect:
- `synth` undefined → `addSampleTagMutator.synthesizeFromSse` is missing or its `tag_id` guard is wrong (Task 2 / pre-existing).
- placeholder not replaced → the corpus `onSuccess`'s `replacePlaceholder` match predicate is wrong (fix in Task 2's `addCorpusSampleTagMutator`).
- `q_units` lost → the corpus `onSuccess` reconstructs the row instead of spreading `...s` (fix in Task 2).

Fix the named source file, not the test.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/test/queue/mutatorOnSseWins.test.ts
git commit -m "$(cat <<'EOF'
Pin the corpus add_tag SSE-wins cross-mutator handoff (#159)

The own-op SSE confirmation builds its response via
addSampleTagMutator.synthesizeFromSse while the corpus mutator's
onSuccess consumes it — a synth/onSuccess split across two mutators
that no other kind has. Pinned as a contract.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: Full verification

**Files:** none — verification only.

- [ ] **Step 1: Run the full frontend queue test suite**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- test/queue/`
Expected: PASS — all queue tests green, including every file this plan touched.

- [ ] **Step 2: Run the frontend build**

Run (from `packages/HimalayaUI/frontend/`): `npm run build`
Expected: `tsc --noEmit` reports no type errors and `vite build` completes.

- [ ] **Step 3: Run the backend route test**

```bash
julia --project=packages/HimalayaUI -e '
using HimalayaUI, Test, HTTP, JSON3, SQLite, DBInterface, Tables
include("packages/HimalayaUI/test/test_http.jl")
include("packages/HimalayaUI/test/test_routes_samples.jl")
include("packages/HimalayaUI/test/test_route_response_shapes.jl")
'
```
Expected: PASS — `samples routes` testsets (incl. the new `source` testset) and the `POST /api/samples/:id/tags → SampleTag exactly` shape test all green.

- [ ] **Step 4: Verify acceptance criteria**

Confirm against `docs/superpowers/specs/2026-05-18-corpus-sample-tag-mutations-design.md` "Acceptance criteria":
- Corpus tag op routes to a sample-tag mutator (Task 6 test).
- Corpus tag add/remove round-trips + patches `corpusSamples` (Tasks 2, 3, 8 tests).
- Corpus tag SSE frame invalidates `corpusSamples` (Task 7 test).
- `resolveMutator` tri-scope; `resolveMutatorForEvent` carries the verified-2-arm comment (Task 6).
- `applyRemoteToCache` invalidates the corpus key for sample frames (Task 7).
- Tag `source` parameterized at the route (Task 1).
- Six-layer contract tests pass (Tasks 1–8).

- [ ] **Step 5: Commit (only if Step 1–3 surfaced a fix)**

If verification was clean, there is nothing to commit. If a fix was needed, commit it with a descriptive message and the `Co-Authored-By` trailer.

---

## Self-review notes

- **Spec coverage:** every spec layer maps to a task — Layer 1 → Task 1; Layer 2 `resolveMutator` → Task 6; Layer 3 corpus mutators → Task 2 (+ hooks Task 5); Layer 4 `applyRemoteToCache` → Task 7; Layers 5–6 → Tasks 2/3/8. The six-layer test matrix → Tasks 1–8. The `resolveMutatorForEvent` 2-arm decision → Task 6 Step 4.
- **Sequencing / contention:** this plan modifies `mutatorRegistry.ts` and `applyRemoteToCache.ts`, which the series-event-kind issues #166/#167/#168 also touch. Per issue-map §3, land this branch and rebase it **last** onto whatever that cluster produced — the modified `add_tag`/`remove_tag` branches may have moved. Not a code change, a merge-process note.
- **Type consistency:** `addCorpusSampleTagMutator` / `removeCorpusSampleTagMutator`, `AddCorpusSampleTagInput` / `RemoveCorpusSampleTagInput`, and `useAddCorpusSampleTag` / `useRemoveCorpusSampleTag` are used with identical names in every task that references them.
- **Task 3 & 4 are test-only** by design — the mutators they exercise are delivered in Task 2. If executed strictly in order they pass immediately at Step 2; their Step 2/3 wording covers the case where a real defect surfaces.
