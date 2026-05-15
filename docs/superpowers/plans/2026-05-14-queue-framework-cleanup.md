# Queue Framework Cleanups Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Resolve issue #129 — co-locate per-kind SSE synth with each mutator, dedupe five placeholder-replace loops behind one helper, strip queue framework metadata from `onSuccess` callbacks via a shared helper.

**Architecture:** Three independent helpers, applied bottom-up. `replacePlaceholder<T>` extracts the optimistic→server row replacement loop. `synthesizeFromSse` becomes an optional `Mutator` interface method; `replayCoordinator` dispatches via `resolveMutatorForEvent(kind, entity_type)` (a sibling to the existing `resolveMutator(op)`). `stripQueueMetadata<T>` returns `{meta, payload}` so `onSuccess` consumers stop hand-destructuring framework fields.

**Tech Stack:** TypeScript, Vitest, TanStack Query v5, React.

**Working dir:** `packages/HimalayaUI/frontend/` (run all `npm` commands from here).

**Run commands:**
- Unit tests: `npm test`
- Single test file: `npm test -- path/to/test.test.ts`
- Type check + build: `npm run build`

**Spec:** [`docs/superpowers/specs/2026-05-14-queue-framework-cleanup-design.md`](../specs/2026-05-14-queue-framework-cleanup-design.md)

---

## File Structure

**Created:**
- `packages/HimalayaUI/frontend/src/lib/queue/replacePlaceholder.ts` — generic placeholder→server-row replace helper
- `packages/HimalayaUI/frontend/src/lib/queue/queueMeta.ts` — `stripQueueMetadata` + `QueueResponseMeta` type
- `packages/HimalayaUI/frontend/test/queue/replacePlaceholder.test.ts`
- `packages/HimalayaUI/frontend/test/queue/queueMeta.test.ts`

**Modified:**
- `packages/HimalayaUI/frontend/src/lib/queue/types.ts` — add `synthesizeFromSse?` and `eventKinds?` to `Mutator<>`
- `packages/HimalayaUI/frontend/src/lib/queue/mutatorRegistry.ts` — add `resolveMutatorForEvent(kind, entity_type)`
- `packages/HimalayaUI/frontend/src/lib/queue/replayCoordinator.ts` — replace switch with mutator dispatch
- `packages/HimalayaUI/frontend/src/lib/queue/mutators/peakAdd.ts` — adopt `replacePlaceholder`, `stripQueueMetadata`, declare `synthesizeFromSse`
- `packages/HimalayaUI/frontend/src/lib/queue/mutators/peakSetExcluded.ts` — adopt `stripQueueMetadata`, declare `synthesizeFromSse`
- `packages/HimalayaUI/frontend/src/lib/queue/mutators/indexGroup.ts` — adopt `stripQueueMetadata` (×2 sites)
- `packages/HimalayaUI/frontend/src/lib/queue/mutators/saveComparison.ts` — declare `synthesizeFromSse`, `eventKinds`
- `packages/HimalayaUI/frontend/src/lib/queue/mutators/deleteComparison.ts` — declare `synthesizeFromSse`
- `packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts` — adopt `replacePlaceholder` (×4 sites); declare `synthesizeFromSse` for add_tag mutators
- `packages/HimalayaUI/frontend/test/queue/sseEventPayload.contract.test.ts` — add coverage assertion (no failing test churn; just an extra it-block)

---

## Phase A — `replacePlaceholder` helper

### Task A1: Write replacePlaceholder unit tests

**Files:**
- Test: `packages/HimalayaUI/frontend/test/queue/replacePlaceholder.test.ts`

- [ ] **Step 1: Write the failing test file**

Create `packages/HimalayaUI/frontend/test/queue/replacePlaceholder.test.ts`:

```ts
import { describe, it, expect } from "vitest";
import { replacePlaceholder } from "../../src/lib/queue/replacePlaceholder";

interface Row { id: number; q: number; label?: string }

describe("replacePlaceholder", () => {
  it("replaces a single negative-id placeholder matching the predicate", () => {
    const list: Row[] = [
      { id: -1, q: 0.1, label: "optimistic" },
      { id:  5, q: 0.2, label: "existing" },
    ];
    const server: Row = { id: 42, q: 0.1, label: "server" };
    const out = replacePlaceholder(list, server, (r) => r.q === 0.1);
    expect(out).toEqual([
      { id: 42, q: 0.1, label: "server" },
      { id:  5, q: 0.2, label: "existing" },
    ]);
  });

  it("appends server item when no placeholder matches", () => {
    const list: Row[] = [{ id: 5, q: 0.2 }];
    const server: Row = { id: 42, q: 0.1 };
    const out = replacePlaceholder(list, server, (r) => r.q === 0.1);
    expect(out).toEqual([{ id: 5, q: 0.2 }, { id: 42, q: 0.1 }]);
  });

  it("dedupes against a concurrent SSE that already inserted the server id", () => {
    const list: Row[] = [
      { id: -1, q: 0.1, label: "optimistic" },
      { id: 42, q: 0.1, label: "from-sse" },
    ];
    const server: Row = { id: 42, q: 0.1, label: "server" };
    const out = replacePlaceholder(list, server, (r) => r.q === 0.1);
    // Placeholder replaced with server; the pre-existing positive-id 42 is dropped.
    expect(out).toEqual([{ id: 42, q: 0.1, label: "server" }]);
  });

  it("replaces only the first placeholder when multiple match", () => {
    const list: Row[] = [
      { id: -1, q: 0.1, label: "a" },
      { id: -2, q: 0.1, label: "b" },
      { id:  5, q: 0.2, label: "keep" },
    ];
    const server: Row = { id: 42, q: 0.1, label: "server" };
    const out = replacePlaceholder(list, server, (r) => r.q === 0.1);
    expect(out).toEqual([
      { id: 42, q: 0.1, label: "server" },
      { id: -2, q: 0.1, label: "b" },
      { id:  5, q: 0.2, label: "keep" },
    ]);
  });

  it("handles an empty list (appends server item)", () => {
    const out = replacePlaceholder<Row>([], { id: 42, q: 0.1 }, () => true);
    expect(out).toEqual([{ id: 42, q: 0.1 }]);
  });

  it("ignores positive-id rows even when matches() returns true", () => {
    const list: Row[] = [
      { id: 7, q: 0.1, label: "real" },
    ];
    const server: Row = { id: 42, q: 0.1, label: "server" };
    const out = replacePlaceholder(list, server, (r) => r.q === 0.1);
    // No placeholder; server appended; pre-existing positive row preserved.
    expect(out).toEqual([
      { id: 7, q: 0.1, label: "real" },
      { id: 42, q: 0.1, label: "server" },
    ]);
  });

  it("preserves placeholder position (does not move the row to the end)", () => {
    const list: Row[] = [
      { id:  1, q: 0.2 },
      { id: -1, q: 0.1, label: "opt" },
      { id:  2, q: 0.3 },
    ];
    const server: Row = { id: 42, q: 0.1, label: "server" };
    const out = replacePlaceholder(list, server, (r) => r.q === 0.1);
    expect(out).toEqual([
      { id:  1, q: 0.2 },
      { id: 42, q: 0.1, label: "server" },
      { id:  2, q: 0.3 },
    ]);
  });

  describe("isDuplicate override", () => {
    interface Tagged { id: number; source: "manual" | "auto" }

    it("preserves rows the override classifies as non-duplicates", () => {
      const list: Tagged[] = [
        { id: 42, source: "auto" },       // same id as server, but different source
        { id: -1, source: "manual" },     // placeholder
      ];
      const server: Tagged = { id: 42, source: "manual" };
      const out = replacePlaceholder(
        list,
        server,
        (r) => r.source === "manual",
        { isDuplicate: (r) => r.source === "manual" && r.id === server.id },
      );
      // Auto row with id=42 survives; placeholder replaced with server row.
      expect(out).toEqual([
        { id: 42, source: "auto" },
        { id: 42, source: "manual" },
      ]);
    });

    it("drops rows the override classifies as duplicates", () => {
      const list: Tagged[] = [
        { id: 42, source: "manual" },     // SSE-race duplicate
        { id: -1, source: "manual" },     // placeholder
      ];
      const server: Tagged = { id: 42, source: "manual" };
      const out = replacePlaceholder(
        list,
        server,
        (r) => r.source === "manual",
        { isDuplicate: (r) => r.source === "manual" && r.id === server.id },
      );
      expect(out).toEqual([{ id: 42, source: "manual" }]);
    });

    it("preserves peakAdd's exact three-row race: placeholder + auto-same-id + manual-different-id", () => {
      // Mirrors the original peakAdd loop's defended scenario:
      // placeholder gets replaced; auto peak sharing id namespace survives;
      // unrelated manual peak survives.
      const list: Tagged[] = [
        { id: -1, source: "manual" },     // placeholder for the q the user clicked
        { id: 42, source: "auto"   },     // auto peak that happens to have id=42
        { id: 15, source: "manual" },     // unrelated existing manual peak
      ];
      const server: Tagged = { id: 42, source: "manual" };
      const out = replacePlaceholder(
        list,
        server,
        (r) => r.source === "manual",
        { isDuplicate: (r) => r.source === "manual" && r.id === server.id },
      );
      expect(out).toEqual([
        { id: 42, source: "manual" },     // placeholder replaced
        { id: 42, source: "auto"   },     // auto survives — id collision is allowed across sources
        { id: 15, source: "manual" },     // unrelated manual survives
      ]);
    });
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npm test -- test/queue/replacePlaceholder.test.ts`
Expected: FAIL with "Failed to load url ... replacePlaceholder" (module does not exist yet).

- [ ] **Step 3: Write minimal implementation**

Create `packages/HimalayaUI/frontend/src/lib/queue/replacePlaceholder.ts`:

```ts
/**
 * Replace a single negative-id optimistic placeholder with the server row,
 * deduping against any concurrent SSE that already inserted the server id.
 *
 * Algorithm:
 *   - Walk the list once.
 *   - The FIRST item with id < 0 whose `matches(item)` predicate returns true
 *     is replaced in place with `serverItem`. Order is preserved.
 *   - Any item where `isDuplicate(item)` returns true is dropped. By default
 *     this matches `item.id === serverItem.id` (the common case: an SSE-wins
 *     race already inserted the server row). Override for cases like
 *     `peakAdd` where the id namespace is shared across kinds and dedup
 *     must be scoped — e.g. only drop the duplicate if it's a manual peak.
 *   - If no placeholder matched, `serverItem` is appended.
 *
 * Invariants:
 *   - At most one row matching `isDuplicate` survives in the output (it
 *     becomes the serverItem if it also matched the placeholder predicate).
 *   - At most one placeholder is consumed per call. Stale placeholders left
 *     over (rare) survive — a later mutator call or invalidate sweeps them.
 *   - Rows that do not match `isDuplicate` are preserved verbatim.
 */
export function replacePlaceholder<T extends { id: number }>(
  list: T[],
  serverItem: T,
  matches: (item: T) => boolean,
  options?: { isDuplicate?: (item: T) => boolean },
): T[] {
  const isDup = options?.isDuplicate ?? ((item: T) => item.id === serverItem.id);
  let replaced = false;
  const out: T[] = [];
  for (const item of list) {
    if (isDup(item)) continue;
    if (!replaced && item.id < 0 && matches(item)) {
      out.push(serverItem);
      replaced = true;
      continue;
    }
    out.push(item);
  }
  if (!replaced) out.push(serverItem);
  return out;
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npm test -- test/queue/replacePlaceholder.test.ts`
Expected: PASS — 7 tests green.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/replacePlaceholder.ts \
        packages/HimalayaUI/frontend/test/queue/replacePlaceholder.test.ts
git commit -m "feat(queue): add replacePlaceholder helper

Generic helper that replaces one negative-id placeholder with the
server row, deduping against concurrent SSE inserts. Will replace
five hand-rolled copies in peakAdd + trivial. (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task A2: Migrate peakAdd to replacePlaceholder

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutators/peakAdd.ts:50-88`

- [ ] **Step 1: Pre-check — run existing peakAdd tests baseline**

Run: `npm test -- test/queue/`
Expected: PASS — record the failing/passing baseline so we can detect regressions.

- [ ] **Step 2: Replace the onSuccess body**

In `packages/HimalayaUI/frontend/src/lib/queue/mutators/peakAdd.ts`, find the `onSuccess` block (currently lines 50–88) and replace **only** the loop body. The destructure + the closing `qc.setQueryData<Exposure>(...)` stay.

Before (current state):

```ts
  onSuccess: (p, response, qc) => {
    // Strip queue-framework metadata before treating the response as a Peak.
    // The route returns a flat `Peak & {event_id, view_row_id, analysis_inputs_hash}`;
    // the cache holds Peak[], so we don't want event_id/etc. polluting it.
    const { event_id: _e, view_row_id: _v, analysis_inputs_hash, ...serverPeak } = response;
    void _e; void _v;
    const peaksKey = queryKeys.peaks(p.exposureId);
    qc.setQueryData<Peak[]>(peaksKey, (old) => {
      const list = old ?? [];
      // Drop the most recent negative-id placeholder for this q (within tol),
      // and dedupe against any concurrent SSE that already inserted the row.
      // Dedup is scoped to MANUAL peaks because auto_peaks.id and
      // peak_curations.id share a namespace on the wire — an auto peak with
      // the same id would otherwise falsely register as "already inserted."
      const next: Peak[] = [];
      let replaced = false;
      const seenManual = new Set<number>();
      for (const pk of list) {
        if (pk.id < 0 && !replaced
            && Math.abs(pk.q - p.q) < peakQTol(p.q)
            && pk.exposure_id === p.exposureId) {
          if (!seenManual.has(serverPeak.id)) {
            next.push(serverPeak);
            seenManual.add(serverPeak.id);
          }
          replaced = true;
          continue;
        }
        if (pk.source === "manual" && seenManual.has(pk.id)) continue;
        next.push(pk);
        if (pk.source === "manual") seenManual.add(pk.id);
      }
      if (!replaced && !seenManual.has(serverPeak.id)) next.push(serverPeak);
      return next;
    });
    qc.setQueryData<Exposure>(queryKeys.exposure(p.exposureId), (old) =>
      old ? { ...old, analysis_inputs_hash } : old);
  },
```

After:

```ts
  onSuccess: (p, response, qc) => {
    // Strip queue-framework metadata before treating the response as a Peak.
    const { event_id: _e, view_row_id: _v, analysis_inputs_hash, ...serverPeak } = response;
    void _e; void _v;
    const peaksKey = queryKeys.peaks(p.exposureId);
    qc.setQueryData<Peak[]>(peaksKey, (old) =>
      // Dedup is scoped to MANUAL peaks: auto_peaks.id and peak_curations.id
      // share a wire namespace, so an auto peak with the same id must survive.
      replacePlaceholder(
        old ?? [],
        serverPeak,
        (pk) => pk.source !== "auto"
             && Math.abs(pk.q - p.q) < peakQTol(p.q)
             && pk.exposure_id === p.exposureId,
        { isDuplicate: (pk) => pk.source === "manual" && pk.id === serverPeak.id },
      ),
    );
    qc.setQueryData<Exposure>(queryKeys.exposure(p.exposureId), (old) =>
      old ? { ...old, analysis_inputs_hash } : old);
  },
```

Add import at the top of the file (alongside the existing `peakQTol` import):

```ts
import { replacePlaceholder } from "../replacePlaceholder";
```

- [ ] **Step 3: Run peakAdd tests + rollbackSymmetry + cache-shape**

Run: `npm test -- test/queue/`
Expected: PASS — all existing tests still green. If `peakAddRollback` or `cache-shape` fail, check the predicate (`source !== "auto"` must keep auto rows in place when their id matches the server id — verify by tracing `seenManual.has(pk.id)` semantics from the old loop).

- [ ] **Step 4: Run typecheck**

Run: `npm run build`
Expected: PASS — no type errors. (The build script runs `tsc --noEmit` before vite.)

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/peakAdd.ts
git commit -m "refactor(queue): migrate peakAdd to replacePlaceholder

The hand-rolled loop is exactly replacePlaceholder's contract.
The 'manual-only dedup' wrinkle moves into the predicate. (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task A3: Migrate postSampleMessage to replacePlaceholder

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts:286-309`

- [ ] **Step 1: Replace the onSuccess body**

In `trivial.ts`, find `postSampleMessageMutator.onSuccess` (lines ~286–309) and replace:

Before:

```ts
  onSuccess: (p, response, qc) => {
    const key = queryKeys.messages(p.sampleId);
    const list = qc.getQueryData<SampleMessage[]>(key) ?? [];
    // Replace the most recent negative-id placeholder for this body, and
    // dedupe against any concurrent SSE that already inserted the real msg.
    const seen = new Set<number>();
    const next: SampleMessage[] = [];
    let replaced = false;
    for (const m of list) {
      if (m.id < 0 && !replaced && m.body === response.body
          && m.sample_id === response.sample_id) {
        if (!seen.has(response.id)) { next.push(response); seen.add(response.id); }
        replaced = true;
        continue;
      }
      if (seen.has(m.id)) continue;
      next.push(m);
      seen.add(m.id);
    }
    if (!replaced && !seen.has(response.id)) next.push(response);
    qc.setQueryData<SampleMessage[]>(key, next);
  },
```

After:

```ts
  onSuccess: (p, response, qc) => {
    const key = queryKeys.messages(p.sampleId);
    qc.setQueryData<SampleMessage[]>(key, (list) =>
      replacePlaceholder(
        list ?? [],
        response,
        (m) => m.body === response.body && m.sample_id === response.sample_id,
      ),
    );
  },
```

Add import at the top of `trivial.ts`:

```ts
import { replacePlaceholder } from "../replacePlaceholder";
```

(Add the import once; it will be reused in A4–A6.)

- [ ] **Step 2: Run tests**

Run: `npm test -- test/queue/`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts
git commit -m "refactor(queue): migrate postSampleMessage to replacePlaceholder (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task A4: Migrate postComparisonMessage to replacePlaceholder

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts:346-367`

- [ ] **Step 1: Replace the onSuccess body**

Before:

```ts
  onSuccess: (p, response, qc) => {
    const key = queryKeys.comparisonMessages(p.comparisonId);
    const list = qc.getQueryData<ComparisonMessage[]>(key) ?? [];
    const seen = new Set<number>();
    const next: ComparisonMessage[] = [];
    let replaced = false;
    for (const m of list) {
      if (m.id < 0 && !replaced && m.body === response.body
          && m.comparison_id === response.comparison_id) {
        if (!seen.has(response.id)) { next.push(response); seen.add(response.id); }
        replaced = true;
        continue;
      }
      if (seen.has(m.id)) continue;
      next.push(m);
      seen.add(m.id);
    }
    if (!replaced && !seen.has(response.id)) next.push(response);
    qc.setQueryData<ComparisonMessage[]>(key, next);
  },
```

After:

```ts
  onSuccess: (p, response, qc) => {
    const key = queryKeys.comparisonMessages(p.comparisonId);
    qc.setQueryData<ComparisonMessage[]>(key, (list) =>
      replacePlaceholder(
        list ?? [],
        response,
        (m) => m.body === response.body && m.comparison_id === response.comparison_id,
      ),
    );
  },
```

- [ ] **Step 2: Run tests**

Run: `npm test -- test/queue/`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts
git commit -m "refactor(queue): migrate postComparisonMessage to replacePlaceholder (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task A5: Migrate addSampleTag to replacePlaceholder (nested list)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts:133-155`

Note: tag lists live nested inside the parent `Sample.tags`. The helper applies to the inner list; the outer `samples.map` survives. This changes behavior slightly — placeholders now keep their position instead of being appended at the end. This is an intentional improvement.

- [ ] **Step 1: Replace the onSuccess body**

Before:

```ts
  onSuccess: (p, response, qc) => {
    // The route emits `{id, sample_id, key, value, source}` (routes_samples.jl)
    // but the SampleTag type omits `sample_id`. Strip it so the cache row
    // matches the type — otherwise tag entries pollute with sample_id.
    const tag: SampleTag = {
      id: response.id, key: response.key, value: response.value,
      source: response.source,
    };
    const samplesKey = queryKeys.samples(p.experimentId);
    qc.setQueryData<Sample[]>(samplesKey, (list) => {
      if (!list) return list;
      return list.map((s) => {
        if (s.id !== p.sampleId) return s;
        const filtered = s.tags.filter((t) =>
          !(t.id < 0 && t.key === p.key && t.value === p.value)
          && t.id !== tag.id,
        );
        return { ...s, tags: [...filtered, tag] };
      });
    });
  },
```

After:

```ts
  onSuccess: (p, response, qc) => {
    // The route emits `{id, sample_id, key, value, source}` (routes_samples.jl)
    // but SampleTag omits `sample_id`. Strip it so the cache row matches type.
    const tag: SampleTag = {
      id: response.id, key: response.key, value: response.value,
      source: response.source,
    };
    const samplesKey = queryKeys.samples(p.experimentId);
    qc.setQueryData<Sample[]>(samplesKey, (list) => {
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
```

- [ ] **Step 2: Run tests**

Run: `npm test -- test/queue/`
Expected: PASS. If any test asserts position-at-end for tags (search `s.tags[s.tags.length - 1]`), update that assertion to find-by-id since the placeholder now keeps its insertion order.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts
git commit -m "refactor(queue): migrate addSampleTag to replacePlaceholder

Side benefit: placeholder tags now keep their insertion position
rather than jumping to the end of the list. (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task A6: Migrate addExposureTag to replacePlaceholder

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts:211-230`

- [ ] **Step 1: Replace the onSuccess body**

Before:

```ts
  onSuccess: (p, response, qc) => {
    const tag: ExposureTag = {
      id: response.id, key: response.key, value: response.value,
      source: response.source,
    };
    rewriteExposureLists(qc, p.sampleId, (list) =>
      list.map((e) => {
        if (e.id !== p.exposureId) return e;
        const filtered = e.tags.filter((t) =>
          !(t.id < 0 && t.key === p.key && t.value === p.value)
          && t.id !== tag.id,
        );
        return { ...e, tags: [...filtered, tag] };
      }),
    );
  },
```

After:

```ts
  onSuccess: (p, response, qc) => {
    const tag: ExposureTag = {
      id: response.id, key: response.key, value: response.value,
      source: response.source,
    };
    rewriteExposureLists(qc, p.sampleId, (list) =>
      list.map((e) =>
        e.id !== p.exposureId
          ? e
          : {
              ...e,
              tags: replacePlaceholder(
                e.tags,
                tag,
                (t) => t.key === p.key && t.value === p.value,
              ),
            },
      ),
    );
  },
```

- [ ] **Step 2: Run tests**

Run: `npm test -- test/queue/`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts
git commit -m "refactor(queue): migrate addExposureTag to replacePlaceholder (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Phase B — `synthesizeFromSse` on the Mutator interface

### Task B1: Extend the Mutator interface

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/types.ts`

- [ ] **Step 1: Add the new fields + type**

Insert after the existing `SseEvent` interface declaration and before `Mutator<>`:

```ts
/**
 * The framework-owned fields that the replay coordinator threads into every
 * synthesized response. Mutators MAY destructure these off via
 * `stripQueueMetadata`. The values come from the SSE frame, not the payload.
 */
export interface QueueResponseMeta {
  event_id: number;
  client_op_id: string | null | undefined;
  analysis_inputs_hash: string | undefined;
  view_row_id?: number;
}
```

Add the two new optional fields to the `Mutator<TInput, TScope, TResponse>` interface (currently lines 128–147):

```ts
export interface Mutator<TInput, TScope, TResponse> {
  kind: OpKind;
  /**
   * The set of SSE event kinds (strings on the wire) this mutator can answer
   * for. Defaults to `[kind]` when omitted. Override when the OpKind name
   * differs from the event kind, e.g. `comparison_save` emits
   * `comparison_created` and `comparison_submitted`.
   */
  eventKinds?: string[];
  onMutate: (payload: FlatPayload<TInput, TScope>, qc: QueryClient) => RollbackContext;
  request: (payload: FlatPayload<TInput, TScope>, signal: AbortSignal) => Promise<TResponse>;
  onSuccess: (payload: FlatPayload<TInput, TScope>, response: TResponse, qc: QueryClient) => void;
  /**
   * Build a synthetic response for the SSE-wins path. When the SSE for our
   * own pending op lands before the HTTP response, the deferred resolves
   * with this value so `onSuccess` can apply the cache effect without
   * waiting for the (now-aborted) HTTP call. Return `undefined` to fall
   * back to the generic `{ ...base, ...payload }` shape.
   *
   * The shape contract for what this returns IS the same shape `onSuccess`
   * expects in its second argument — co-locating the two prevents drift.
   */
  synthesizeFromSse?: (remote: SseEvent, base: QueueResponseMeta) => TResponse | undefined;
  affectsExposurePeaks?: (payload: FlatPayload<TInput, TScope>, exposureId: number) => boolean;
  /**
   * When set, an HTTP 404 response is treated as a no-op success: the
   * framework skips both the toast and the rollback, leaving the optimistic
   * cache state in place. Use on idempotent remove/delete operations where
   * "the server doesn't have it" is the desired end state.
   *
   * Motivation: under 5xx-then-retry, a successful first attempt deletes the
   * row; the client retries on the missed response and gets a 404 from a
   * server that's already done the work. Without this flag, the rollback
   * re-inserts the row that the server (and SSE) have already removed,
   * producing a phantom row visible until the next refetch.
   */
  treats404AsSuccess?: boolean;
}
```

- [ ] **Step 2: Run typecheck**

Run: `npm run build`
Expected: PASS — `synthesizeFromSse?` and `eventKinds?` are optional, so no existing mutator must change.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/types.ts
git commit -m "feat(queue): extend Mutator with synthesizeFromSse + eventKinds

Optional fields — no existing mutators need to change. Sets up the
shape contract for moving per-kind synth from replayCoordinator into
the mutator. (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task B2: Add resolveMutatorForEvent dispatcher

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutatorRegistry.ts`
- Test: `packages/HimalayaUI/frontend/test/queue/mutatorRegistry.test.ts` (extend)

The dispatcher must handle:
- 1:1 OpKind ↔ event kind for most mutators (`peak_added`, etc.)
- `add_tag` event splits by `entity_type`: `"sample"` → `addSampleTagMutator`, `"exposure"` → `addExposureTagMutator`. Same for `remove_tag` → remove variants.
- `post_message` event splits by `entity_type`: `"sample"` → `postSampleMessageMutator`, `"comparison"` → `postComparisonMessageMutator`.
- `comparison_created` and `comparison_submitted` both → `saveComparisonMutator`.
- `analyze_run` event → `reanalyzeExposureMutator`.
- `speculative_deleted` event → `deleteIndexMutator` (the user gesture is `delete_index`, the wire event is `speculative_deleted`).

- [ ] **Step 1: Update test file imports**

The test file currently imports trivial mutators (addSampleTag, removeSampleTag, addExposureTag, removeExposureTag, updateSample, postSampleMessage, setExposureStatus, selectExposure), peakAdd, peakRemove, peakExclude, peakUnexclude, the three indexGroup variants, createSpeculative, and reanalyzeExposure. **Add three new imports:**

```ts
// Inside the existing trivial-mutators import block: add postComparisonMessageMutator.
import {
  addSampleTagMutator,
  removeSampleTagMutator,
  addExposureTagMutator,
  removeExposureTagMutator,
  updateSampleMutator,
  postSampleMessageMutator,
  postComparisonMessageMutator,    // ← new
  setExposureStatusMutator,
  selectExposureMutator,
} from "../../src/lib/queue/mutators/trivial";

// New separate imports near the existing mutator imports:
import { saveComparisonMutator }   from "../../src/lib/queue/mutators/saveComparison";
import { deleteComparisonMutator } from "../../src/lib/queue/mutators/deleteComparison";
```

Also add the named export of the new function to the existing top-of-file import:

```ts
import { resolveMutator, resolveMutatorForEvent } from "../../src/lib/queue/mutatorRegistry";
```

- [ ] **Step 2: Add the unit test**

Append to `packages/HimalayaUI/frontend/test/queue/mutatorRegistry.test.ts`:

```ts
describe("resolveMutatorForEvent", () => {
  // Imports must match the existing test file's mutator import block.
  it("routes peak_added to peakAddMutator", () => {
    expect(resolveMutatorForEvent("peak_added", "exposure")).toBe(peakAddMutator);
  });

  it("routes peak_excluded and peak_unexcluded to their dedicated mutators", () => {
    expect(resolveMutatorForEvent("peak_excluded",   "exposure")).toBe(peakExcludeMutator);
    expect(resolveMutatorForEvent("peak_unexcluded", "exposure")).toBe(peakUnexcludeMutator);
  });

  it("splits add_tag by entity_type", () => {
    expect(resolveMutatorForEvent("add_tag", "sample"))  .toBe(addSampleTagMutator);
    expect(resolveMutatorForEvent("add_tag", "exposure")).toBe(addExposureTagMutator);
  });

  it("splits remove_tag by entity_type", () => {
    expect(resolveMutatorForEvent("remove_tag", "sample"))  .toBe(removeSampleTagMutator);
    expect(resolveMutatorForEvent("remove_tag", "exposure")).toBe(removeExposureTagMutator);
  });

  it("splits post_message by entity_type", () => {
    expect(resolveMutatorForEvent("post_message", "sample"))    .toBe(postSampleMessageMutator);
    expect(resolveMutatorForEvent("post_message", "comparison")).toBe(postComparisonMessageMutator);
  });

  it("routes both comparison_created and comparison_submitted to saveComparison", () => {
    expect(resolveMutatorForEvent("comparison_created",   "comparison")).toBe(saveComparisonMutator);
    expect(resolveMutatorForEvent("comparison_submitted", "comparison")).toBe(saveComparisonMutator);
  });

  it("routes comparison_deleted to deleteComparison", () => {
    expect(resolveMutatorForEvent("comparison_deleted", "comparison")).toBe(deleteComparisonMutator);
  });

  it("routes analyze_run to reanalyzeExposure", () => {
    expect(resolveMutatorForEvent("analyze_run", "exposure")).toBe(reanalyzeExposureMutator);
  });

  it("returns undefined for unknown event kinds", () => {
    expect(resolveMutatorForEvent("not_a_kind", "exposure")).toBeUndefined();
  });
});

describe("resolveMutator ↔ resolveMutatorForEvent consistency", () => {
  // Cross-check: for each mutator, the event kind it canonically emits
  // (from `eventKinds[0]` or `kind`) must resolve back to the same mutator.
  // Catches drift if someone adds a mutator + updates resolveMutator but
  // forgets resolveMutatorForEvent (or vice versa).
  const cases = [
    { mutator: peakAddMutator,          eventKind: "peak_added",          entityType: "exposure"   },
    { mutator: peakRemoveMutator,       eventKind: "peak_removed",        entityType: "exposure"   },
    { mutator: peakExcludeMutator,      eventKind: "peak_excluded",       entityType: "exposure"   },
    { mutator: peakUnexcludeMutator,    eventKind: "peak_unexcluded",     entityType: "exposure"   },
    { mutator: addIndexToGroupMutator,  eventKind: "index_confirmed",     entityType: "exposure"   },
    { mutator: removeIndexFromGroupMutator, eventKind: "index_unconfirmed", entityType: "exposure" },
    { mutator: deleteIndexMutator,      eventKind: "speculative_deleted", entityType: "exposure"   },
    { mutator: createSpeculativeMutator, eventKind: "speculative_created", entityType: "exposure"  },
    { mutator: reanalyzeExposureMutator, eventKind: "analyze_run",        entityType: "exposure"   },
    { mutator: saveComparisonMutator,   eventKind: "comparison_created",  entityType: "comparison" },
    { mutator: saveComparisonMutator,   eventKind: "comparison_submitted", entityType: "comparison" },
    { mutator: deleteComparisonMutator, eventKind: "comparison_deleted",  entityType: "comparison" },
    { mutator: addSampleTagMutator,     eventKind: "add_tag",             entityType: "sample"     },
    { mutator: addExposureTagMutator,   eventKind: "add_tag",             entityType: "exposure"   },
    { mutator: removeSampleTagMutator,  eventKind: "remove_tag",          entityType: "sample"     },
    { mutator: removeExposureTagMutator, eventKind: "remove_tag",         entityType: "exposure"   },
    { mutator: postSampleMessageMutator, eventKind: "post_message",       entityType: "sample"     },
    { mutator: postComparisonMessageMutator, eventKind: "post_message",   entityType: "comparison" },
    { mutator: updateSampleMutator,     eventKind: "update_sample",       entityType: "sample"     },
    { mutator: setExposureStatusMutator, eventKind: "set_exposure_status", entityType: "exposure"  },
    { mutator: selectExposureMutator,   eventKind: "select_exposure",     entityType: "exposure"   },
  ];

  it.each(cases)(
    "$eventKind / $entityType resolves to $mutator.kind",
    ({ mutator, eventKind, entityType }) => {
      expect(resolveMutatorForEvent(eventKind, entityType)).toBe(mutator);
    },
  );
});
```

- [ ] **Step 3: Run test to verify it fails**

Run: `npm test -- test/queue/mutatorRegistry.test.ts`
Expected: FAIL — `resolveMutatorForEvent` is not exported.

- [ ] **Step 4: Implement resolveMutatorForEvent**

Append to `packages/HimalayaUI/frontend/src/lib/queue/mutatorRegistry.ts`:

```ts
/**
 * Resolve the mutator that owns a given SSE event kind. Used by the replay
 * coordinator to dispatch `synthesizeFromSse` per event kind.
 *
 * `entity_type` is required because three event kinds (`add_tag`,
 * `remove_tag`, `post_message`) are shared across mutators that differ only
 * by the entity scope they operate on.
 *
 * Returns undefined for unknown event kinds — replayCoordinator will fall
 * back to the generic `{ ...base, ...payload }` shape.
 */
export function resolveMutatorForEvent(
  eventKind: string,
  entityType: string,
): Mutator<any, any, any> | undefined {
  switch (eventKind) {
    case "peak_added":         return peakAddMutator;
    case "peak_removed":       return peakRemoveMutator;
    case "peak_excluded":      return peakExcludeMutator;
    case "peak_unexcluded":    return peakUnexcludeMutator;
    case "index_confirmed":    return addIndexToGroupMutator;
    case "index_unconfirmed":  return removeIndexFromGroupMutator;
    case "speculative_created":return createSpeculativeMutator;
    case "speculative_deleted":return deleteIndexMutator;
    case "analyze_run":        return reanalyzeExposureMutator;
    case "comparison_created":
    case "comparison_submitted":
      return saveComparisonMutator;
    case "comparison_deleted": return deleteComparisonMutator;
    case "update_sample":      return updateSampleMutator;
    case "set_exposure_status":return setExposureStatusMutator;
    case "select_exposure":    return selectExposureMutator;
    case "add_tag":
      return entityType === "sample" ? addSampleTagMutator : addExposureTagMutator;
    case "remove_tag":
      return entityType === "sample" ? removeSampleTagMutator : removeExposureTagMutator;
    case "post_message":
      return entityType === "comparison" ? postComparisonMessageMutator : postSampleMessageMutator;
    default:
      return undefined;
  }
}
```

The imports at the top of `mutatorRegistry.ts` already include every mutator name used above (the existing `resolveMutator` switch references them). No new imports needed unless verifying with the existing file shows a missing one.

- [ ] **Step 5: Run test to verify it passes**

Run: `npm test -- test/queue/mutatorRegistry.test.ts`
Expected: PASS — all new describe block assertions green, plus the consistency cross-check.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutatorRegistry.ts \
        packages/HimalayaUI/frontend/test/queue/mutatorRegistry.test.ts
git commit -m "feat(queue): add resolveMutatorForEvent dispatcher

Maps SSE event kind + entity_type to the owning mutator, the inverse
index of resolveMutator. Used next to dispatch synthesizeFromSse from
replayCoordinator. (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task B3: Wire replayCoordinator to dispatch via the mutator

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/replayCoordinator.ts:135-212`

This is the riskiest commit. Keep the old switch as a fallback so the existing behavior is preserved exactly — the next commits will move the per-kind branches into individual mutators and reduce the fallback to a default-only path.

- [ ] **Step 1: Add mutator dispatch BEFORE the existing switch**

Replace the `synthesizeResponseFromSse` body (lines 135–212):

```ts
function synthesizeResponseFromSse(remote: SseEvent): unknown {
  const payload = (remote.payload as Record<string, unknown> | undefined) ?? {};
  const base: QueueResponseMeta = {
    event_id: remote.id,
    client_op_id: remote.client_op_id,
    analysis_inputs_hash: remote.post_state?.analysis_inputs_hash,
  };

  // Phase B preferred path: ask the owning mutator to build the synth shape.
  const mutator = resolveMutatorForEvent(remote.kind, remote.entity_type);
  const synth = mutator?.synthesizeFromSse?.(remote, base);
  if (synth !== undefined) return synth;

  // Legacy fallback: per-kind branches still here until each mutator owns
  // its own synthesizeFromSse (tasks B4–B8). Once those land, this block
  // collapses to just `return { ...base, ...payload };`.
  if (remote.kind === "peak_added") {
    const peakId = payload.peak_curation_id as number | undefined;
    return {
      ...base,
      id: peakId,
      exposure_id: remote.entity_id,
      q: payload.q as number,
      intensity: null,
      prominence: null,
      sharpness: null,
      source: "manual",
      excluded: false,
      view_row_id: peakId,
    };
  }
  if (remote.kind === "add_tag") {
    return {
      ...base,
      id: payload.tag_id,
      key: payload.key,
      value: payload.value,
      source: "manual",
    };
  }
  if (remote.kind === "comparison_created" || remote.kind === "comparison_submitted") {
    return {
      ...base,
      ...payload,
      id: remote.entity_id,
    };
  }
  if (remote.kind === "comparison_deleted") {
    return { ...base, id: remote.entity_id };
  }
  if (remote.kind === "peak_excluded" || remote.kind === "peak_unexcluded") {
    return {
      ...base,
      id: payload.auto_peak_id,
      q: payload.q,
      source: "auto",
      excluded: remote.kind === "peak_excluded",
    };
  }
  return { ...base, ...payload };
}
```

Add imports at the top of `replayCoordinator.ts`:

```ts
import { resolveMutatorForEvent } from "./mutatorRegistry";
import type { QueueResponseMeta } from "./types";
```

- [ ] **Step 2: Run all queue tests**

Run: `npm test -- test/queue/`
Expected: PASS — behavior unchanged (legacy fallback handles every kind). `sseEventPayload.contract.test.ts`, `replayCoordinator.test.ts`, `mutatorOnSseWins.test.tsx`, `sseResolvedThenLate404.test.tsx` should all stay green.

- [ ] **Step 3: Run typecheck**

Run: `npm run build`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/replayCoordinator.ts
git commit -m "feat(queue): dispatch synthesizeFromSse via mutator first

replayCoordinator now asks the owning mutator to synthesize the SSE
response, falling back to the legacy per-kind switch when the mutator
hasn't declared synthesizeFromSse yet. Behavior identical. (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task B4: peakAdd owns its synth

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutators/peakAdd.ts`

- [ ] **Step 1: Add `synthesizeFromSse` to peakAddMutator**

Insert into the mutator object literal, alongside `kind`, `onMutate`, `request`, `onSuccess`:

```ts
  synthesizeFromSse: (remote, base) => {
    const payload = (remote.payload as Record<string, unknown> | undefined) ?? {};
    const peakId = payload.peak_curation_id as number | undefined;
    if (peakId === undefined) return undefined;
    return {
      ...base,
      id: peakId,
      exposure_id: remote.entity_id,
      q: payload.q as number,
      intensity: null,
      prominence: null,
      sharpness: null,
      source: "manual",
      excluded: false,
      view_row_id: peakId,
    } as PeakAddResponse;
  },
```

(`PeakAddResponse` is the existing imported TResponse type at the top of the file — `import type { Peak, PeakAddResponse, Exposure, AuthOpts } from "../../../api"`. No new imports needed for the cast.)

- [ ] **Step 2: Remove the legacy peak_added branch from replayCoordinator**

In `replayCoordinator.ts` `synthesizeResponseFromSse`, delete:

```ts
  if (remote.kind === "peak_added") {
    const peakId = payload.peak_curation_id as number | undefined;
    return {
      ...base,
      id: peakId,
      exposure_id: remote.entity_id,
      q: payload.q as number,
      intensity: null,
      prominence: null,
      sharpness: null,
      source: "manual",
      excluded: false,
      view_row_id: peakId,
    };
  }
```

- [ ] **Step 3: Run tests**

Run: `npm test -- test/queue/`
Expected: PASS. Particular focus on `mutatorOnSseWins.test.tsx`, `sseEventPayload.contract.test.ts`, and `replayCoordinator.test.ts`.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/peakAdd.ts \
        packages/HimalayaUI/frontend/src/lib/queue/replayCoordinator.ts
git commit -m "refactor(queue): peakAdd owns its SSE synth (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task B5: Tag mutators own their synth (add_tag for sample + exposure)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts`

Both `addSampleTagMutator` and `addExposureTagMutator` produce the same synth shape — only the `entity_type` discriminates which one runs (handled by `resolveMutatorForEvent`). The bodies are identical; declare both.

- [ ] **Step 1: Add `synthesizeFromSse` to both add-tag mutators**

For `addSampleTagMutator`, add:

```ts
  synthesizeFromSse: (remote, base) => {
    const payload = (remote.payload as Record<string, unknown> | undefined) ?? {};
    if (payload.tag_id === undefined) return undefined;
    return {
      ...base,
      id: payload.tag_id as number,
      key: payload.key as string,
      value: (payload.value as string) ?? null,
      source: "manual",
    } as SampleTag;
  },
```

For `addExposureTagMutator`, add:

```ts
  synthesizeFromSse: (remote, base) => {
    const payload = (remote.payload as Record<string, unknown> | undefined) ?? {};
    if (payload.tag_id === undefined) return undefined;
    return {
      ...base,
      id: payload.tag_id as number,
      key: payload.key as string,
      value: (payload.value as string) ?? null,
      source: "manual",
    } as ExposureTag;
  },
```

`SampleTag` and `ExposureTag` are the actual `TResponse` types declared on each mutator (`Mutator<…, …, SampleTag>` and `Mutator<…, …, ExposureTag>` respectively, with the types imported from `../../../api`). The fact that the synth object also carries `event_id` / `analysis_inputs_hash` / `client_op_id` is irrelevant to the cast — those framework fields are stripped in `onSuccess` (the route response is `SampleTag & sample_id & framework-meta`, and `onSuccess` already filters the inner shape before writing the cache).

- [ ] **Step 2: Remove the legacy `add_tag` branch from replayCoordinator**

Delete:

```ts
  if (remote.kind === "add_tag") {
    return {
      ...base,
      id: payload.tag_id,
      key: payload.key,
      value: payload.value,
      source: "manual",
    };
  }
```

- [ ] **Step 3: Run tests**

Run: `npm test -- test/queue/`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/trivial.ts \
        packages/HimalayaUI/frontend/src/lib/queue/replayCoordinator.ts
git commit -m "refactor(queue): add-tag mutators own their SSE synth (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task B6: saveComparison owns its synth (two event kinds)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutators/saveComparison.ts`

- [ ] **Step 1: Declare eventKinds + synthesizeFromSse**

Add to `saveComparisonMutator`:

```ts
  eventKinds: ["comparison_created", "comparison_submitted"],
  synthesizeFromSse: (remote, base) => {
    const payload = (remote.payload as Record<string, unknown> | undefined) ?? {};
    // Partial Comparison shape — onSuccess's looksFull detector trips the
    // invalidate fallback because `members` is absent. id is required so
    // the cache-key targeting works.
    return {
      ...base,
      ...payload,
      id: remote.entity_id,
    } as Comparison;
  },
```

(`Comparison` is the response type — confirm in the existing file.)

- [ ] **Step 2: Remove the legacy comparison_created/comparison_submitted branch**

Delete from `replayCoordinator.ts`:

```ts
  if (remote.kind === "comparison_created" || remote.kind === "comparison_submitted") {
    return {
      ...base,
      ...payload,
      id: remote.entity_id,
    };
  }
```

- [ ] **Step 3: Run tests**

Run: `npm test -- test/queue/`
Expected: PASS. `saveComparison.test.tsx` covers the looksFull-vs-invalidate branch.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/saveComparison.ts \
        packages/HimalayaUI/frontend/src/lib/queue/replayCoordinator.ts
git commit -m "refactor(queue): saveComparison owns its SSE synth (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task B7: deleteComparison owns its synth

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutators/deleteComparison.ts`

- [ ] **Step 1: Declare synthesizeFromSse**

Add to `deleteComparisonMutator`:

```ts
  synthesizeFromSse: (remote, base) => ({
    ...base,
    id: remote.entity_id,
    deleted: true,
  } as DeleteResponse),
```

(`DeleteResponse` is the file-local type declared in `deleteComparison.ts:29` — `{ id: number; deleted: boolean; event_id: number }`. The synth must include `deleted: true` to satisfy the type; the server-side semantics of the event already imply deletion succeeded.)

- [ ] **Step 2: Remove the legacy comparison_deleted branch**

Delete from `replayCoordinator.ts`:

```ts
  if (remote.kind === "comparison_deleted") {
    return { ...base, id: remote.entity_id };
  }
```

- [ ] **Step 3: Run tests**

Run: `npm test -- test/queue/`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/deleteComparison.ts \
        packages/HimalayaUI/frontend/src/lib/queue/replayCoordinator.ts
git commit -m "refactor(queue): deleteComparison owns its SSE synth (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task B8: peakSetExcluded owns its synth (two event kinds)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutators/peakSetExcluded.ts`

Both `peakExcludeMutator` and `peakUnexcludeMutator` are exported from this file. The synth shape differs only by the boolean field. The mapping by event kind is already in `resolveMutatorForEvent` (B2). Each mutator declares only its own kind.

- [ ] **Step 1: Add synthesizeFromSse to both**

For `peakExcludeMutator`, add:

```ts
  synthesizeFromSse: (remote, base) => {
    const payload = (remote.payload as Record<string, unknown> | undefined) ?? {};
    if (payload.auto_peak_id === undefined) return undefined;
    return {
      ...base,
      id: payload.auto_peak_id as number,
      q: payload.q as number,
      source: "auto",
      excluded: true,
    } as PeakUpdatedResponse;
  },
```

For `peakUnexcludeMutator`, add the same body but `excluded: false`.

(`PeakUpdatedResponse` is the existing TResponse type — confirm in the file.)

- [ ] **Step 2: Remove the legacy peak_excluded/peak_unexcluded branch**

Delete from `replayCoordinator.ts`:

```ts
  if (remote.kind === "peak_excluded" || remote.kind === "peak_unexcluded") {
    return {
      ...base,
      id: payload.auto_peak_id,
      q: payload.q,
      source: "auto",
      excluded: remote.kind === "peak_excluded",
    };
  }
```

- [ ] **Step 3: Run tests**

Run: `npm test -- test/queue/`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/peakSetExcluded.ts \
        packages/HimalayaUI/frontend/src/lib/queue/replayCoordinator.ts
git commit -m "refactor(queue): peakSetExcluded owns its SSE synth (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task B9: Replace synthesizeResponseFromSse body with mutator-only dispatch

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/replayCoordinator.ts`

After B4–B8 land, every per-kind branch has been migrated. The function reduces to ~10 lines.

- [ ] **Step 1: Replace the function body**

Final state of `synthesizeResponseFromSse`:

```ts
/**
 * Build a synthetic response object mirroring what the HTTP route would
 * have returned. The deferred is awaited inside useQueueMutation's
 * mutationFn; the resolution lets the mutation transition to "success"
 * without the HTTP call having returned. event_id is lifted from the SSE
 * frame; analysis_inputs_hash from post_state if present; remaining
 * fields come from the owning mutator's synthesizeFromSse method.
 *
 * Dispatch: lookup by (event_kind, entity_type) via resolveMutatorForEvent,
 * then ask the mutator to build its shape. Falls back to a generic
 * `{...base, ...payload}` for forward-scaffolded kinds that have no client
 * mutator yet (post_message, set_exposure_status, update_sample,
 * select_exposure, remove_tag). These emit SSE today (applyRemoteToCache
 * merges them) but no UI gesture queues them, so the generic shape suffices.
 * When a future plan wires a UI mutation that emits one of these kinds, add
 * `synthesizeFromSse` to the corresponding mutator and the generic fallback
 * stops being hit for that kind automatically.
 *
 * Ordering note: `applyPostStateOnly(remote)` runs BEFORE the deferred is
 * resolved with this synth (see handleRemoteEvent). That ordering keeps the
 * post-mutation indices cache fresh before the mutator's onSuccess fires.
 * This refactor does NOT change that ordering — it only changes how the
 * synth shape is produced.
 */
function synthesizeResponseFromSse(remote: SseEvent): unknown {
  const base: QueueResponseMeta = {
    event_id: remote.id,
    client_op_id: remote.client_op_id,
    analysis_inputs_hash: remote.post_state?.analysis_inputs_hash,
  };
  const mutator = resolveMutatorForEvent(remote.kind, remote.entity_type);
  const synth = mutator?.synthesizeFromSse?.(remote, base);
  if (synth !== undefined) return synth;
  const payload = (remote.payload as Record<string, unknown> | undefined) ?? {};
  return { ...base, ...payload };
}
```

Replace the function's original doc-comment (lines 120–134 in the legacy file) with the doc-comment above so the new dispatch model is documented.

- [ ] **Step 2: Run all queue tests**

Run: `npm test -- test/queue/`
Expected: PASS.

- [ ] **Step 3: Run typecheck**

Run: `npm run build`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/replayCoordinator.ts
git commit -m "refactor(queue): replayCoordinator drops legacy synth switch

Every kind now has either a mutator-owned synthesizeFromSse or hits
the generic fallback. The 80-line per-kind switch is gone. (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task B10: Add coverage assertion for synth contract

**Files:**
- Modify: `packages/HimalayaUI/frontend/test/queue/sseEventPayload.contract.test.ts`

Pin the invariant: every SSE event kind that the backend can emit is reachable through `resolveMutatorForEvent`, and the kinds that historically had bespoke synth still produce a non-undefined result.

- [ ] **Step 1: Add the assertion**

Append a new `describe` block:

```ts
describe("synthesizeFromSse coverage (resolveMutatorForEvent contract)", () => {
  // The exact set of kinds that had a bespoke branch in the legacy switch.
  // Each MUST resolve to a mutator that returns non-undefined for a
  // minimally-shaped SseEvent.
  const cases: Array<{ kind: string; entity_type: string; payload: object }> = [
    { kind: "peak_added",         entity_type: "exposure",   payload: { peak_curation_id: 99, q: 0.123 } },
    { kind: "add_tag",            entity_type: "sample",     payload: { tag_id: 7, key: "tag", value: "v" } },
    { kind: "add_tag",            entity_type: "exposure",   payload: { tag_id: 7, key: "tag", value: "v" } },
    { kind: "comparison_created", entity_type: "comparison", payload: { title: "T" } },
    { kind: "comparison_submitted", entity_type: "comparison", payload: { title: "T" } },
    { kind: "comparison_deleted", entity_type: "comparison", payload: {} },
    { kind: "peak_excluded",      entity_type: "exposure",   payload: { auto_peak_id: 1, q: 0.1 } },
    { kind: "peak_unexcluded",    entity_type: "exposure",   payload: { auto_peak_id: 1, q: 0.1 } },
  ];
  it.each(cases)("$kind/$entity_type resolves and synthesizes", ({ kind, entity_type, payload }) => {
    const mutator = resolveMutatorForEvent(kind, entity_type);
    expect(mutator).toBeDefined();
    const synth = mutator!.synthesizeFromSse?.(
      { id: 1, kind, entity_type, entity_id: 42, payload } as SseEvent,
      { event_id: 1, client_op_id: "x", analysis_inputs_hash: undefined },
    );
    expect(synth).toBeDefined();
  });

  // Forward-scaffolded kinds: their mutators MUST NOT declare synthesizeFromSse.
  // If a future contributor accidentally adds a half-baked synth to one of
  // these mutators before the rest of the pipeline is ready, this assertion
  // catches it. These kinds emit SSE today but no client mutator queues them,
  // so synthesis falls through to the generic `{...base, ...payload}` shape.
  const forwardScaffolded: Array<{ kind: string; entity_type: string }> = [
    { kind: "post_message",       entity_type: "sample"     },
    { kind: "post_message",       entity_type: "comparison" },
    { kind: "set_exposure_status", entity_type: "exposure"   },
    { kind: "update_sample",      entity_type: "sample"     },
    { kind: "select_exposure",    entity_type: "exposure"   },
    { kind: "remove_tag",         entity_type: "sample"     },
    { kind: "remove_tag",         entity_type: "exposure"   },
  ];
  it.each(forwardScaffolded)(
    "$kind/$entity_type stays on the generic fallback (no mutator.synthesizeFromSse)",
    ({ kind, entity_type }) => {
      const mutator = resolveMutatorForEvent(kind, entity_type);
      expect(mutator).toBeDefined();
      expect(mutator!.synthesizeFromSse).toBeUndefined();
    },
  );
});
```

Add the missing import at the top of the file:

```ts
import { resolveMutatorForEvent } from "../../src/lib/queue/mutatorRegistry";
```

- [ ] **Step 2: Run the test**

Run: `npm test -- test/queue/sseEventPayload.contract.test.ts`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/test/queue/sseEventPayload.contract.test.ts
git commit -m "test(queue): pin synthesizeFromSse contract (#129)

Asserts every event kind with bespoke synth resolves to a mutator
that returns a non-undefined synth.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Phase C — `stripQueueMetadata` helper

### Task C1: Write stripQueueMetadata helper + tests

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/queue/queueMeta.ts`
- Test: `packages/HimalayaUI/frontend/test/queue/queueMeta.test.ts`

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/frontend/test/queue/queueMeta.test.ts`:

```ts
import { describe, it, expect } from "vitest";
import { stripQueueMetadata } from "../../src/lib/queue/queueMeta";

describe("stripQueueMetadata", () => {
  it("splits response into meta + payload", () => {
    const response = {
      event_id: 42,
      view_row_id: 99,
      analysis_inputs_hash: "abc",
      client_op_id: "op-1",
      id: 7,
      q: 0.12,
      source: "manual" as const,
    };
    const { meta, payload } = stripQueueMetadata(response);
    expect(meta).toEqual({
      event_id: 42,
      view_row_id: 99,
      analysis_inputs_hash: "abc",
      client_op_id: "op-1",
    });
    expect(payload).toEqual({ id: 7, q: 0.12, source: "manual" });
  });

  it("leaves undefined meta fields as undefined", () => {
    const response = { id: 7 };
    const { meta, payload } = stripQueueMetadata(response);
    expect(meta).toEqual({
      event_id: undefined,
      view_row_id: undefined,
      analysis_inputs_hash: undefined,
      client_op_id: undefined,
    });
    expect(payload).toEqual({ id: 7 });
  });

  it("does not mutate the input", () => {
    const response = { event_id: 1, id: 2 };
    stripQueueMetadata(response);
    expect(response).toEqual({ event_id: 1, id: 2 });
  });
});
```

- [ ] **Step 2: Run to verify it fails**

Run: `npm test -- test/queue/queueMeta.test.ts`
Expected: FAIL — module does not exist.

- [ ] **Step 3: Write the helper**

Create `packages/HimalayaUI/frontend/src/lib/queue/queueMeta.ts`:

```ts
import type { QueueResponseMeta } from "./types";

/**
 * Split a queue mutator response into its framework-owned metadata and
 * the mutator-specific payload. Use in `onSuccess` so the per-mutator
 * cache write doesn't have to spell out every framework field by hand —
 * and a new field added to `QueueResponseMeta` later flows through
 * automatically.
 *
 *   const { meta, payload } = stripQueueMetadata(response);
 *   qc.setQueryData(...payload row...)
 *   writeExposureHash(qc, exposureId, meta.analysis_inputs_hash);
 */
export function stripQueueMetadata<T extends Partial<QueueResponseMeta>>(
  response: T,
): { meta: QueueResponseMeta; payload: Omit<T, keyof QueueResponseMeta> } {
  const {
    event_id,
    view_row_id,
    analysis_inputs_hash,
    client_op_id,
    ...payload
  } = response as T & QueueResponseMeta;
  return {
    meta: { event_id, view_row_id, analysis_inputs_hash, client_op_id },
    payload: payload as Omit<T, keyof QueueResponseMeta>,
  };
}
```

- [ ] **Step 4: Run to verify it passes**

Run: `npm test -- test/queue/queueMeta.test.ts`
Expected: PASS — 3 tests green.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/queueMeta.ts \
        packages/HimalayaUI/frontend/test/queue/queueMeta.test.ts
git commit -m "feat(queue): add stripQueueMetadata helper (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task C2: Migrate peakAdd to stripQueueMetadata

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutators/peakAdd.ts`

- [ ] **Step 1: Replace the destructure**

In `onSuccess`, change:

```ts
    const { event_id: _e, view_row_id: _v, analysis_inputs_hash, ...serverPeak } = response;
    void _e; void _v;
```

To:

```ts
    const { meta, payload: serverPeak } = stripQueueMetadata(response);
```

And update the line that writes the exposure hash:

```ts
    qc.setQueryData<Exposure>(queryKeys.exposure(p.exposureId), (old) =>
      old ? { ...old, analysis_inputs_hash: meta.analysis_inputs_hash } : old);
```

Add import:

```ts
import { stripQueueMetadata } from "../queueMeta";
```

- [ ] **Step 2: Run tests**

Run: `npm test -- test/queue/`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/peakAdd.ts
git commit -m "refactor(queue): peakAdd uses stripQueueMetadata (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task C3: Migrate peakSetExcluded to stripQueueMetadata

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutators/peakSetExcluded.ts:50-70`

- [ ] **Step 1: Replace the destructure**

Change:

```ts
      const {
        // eslint-disable-next-line @typescript-eslint/no-unused-vars
        event_id, view_row_id, analysis_inputs_hash, client_op_id,
        ...peakFields
      } = response as PeakUpdatedResponse & { client_op_id?: string };
```

To:

```ts
      const { meta, payload: peakFields } = stripQueueMetadata(response);
```

And update the exposure-hash write:

```ts
      qc.setQueryData<Exposure>(queryKeys.exposure(p.exposureId), (old) =>
        old ? { ...old, analysis_inputs_hash: meta.analysis_inputs_hash } : old);
```

Add import:

```ts
import { stripQueueMetadata } from "../queueMeta";
```

- [ ] **Step 2: Run tests**

Run: `npm test -- test/queue/`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/peakSetExcluded.ts
git commit -m "refactor(queue): peakSetExcluded uses stripQueueMetadata (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task C4: Migrate both indexGroup onSuccess sites to stripQueueMetadata

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutators/indexGroup.ts:69-91, 118-135`

- [ ] **Step 1: Replace both destructures**

In both `addIndexToGroupMutator.onSuccess` and `removeIndexFromGroupMutator.onSuccess`, change:

```ts
    const { event_id: _e, view_row_id: _v, ...row } = response;
    void _e; void _v;
```

To (in both places):

```ts
    const { payload: row } = stripQueueMetadata(response);
```

(The `meta` object is unused in indexGroup — its onSuccess paths don't write `analysis_inputs_hash`. Just discard the meta key.)

Add import:

```ts
import { stripQueueMetadata } from "../queueMeta";
```

- [ ] **Step 2: Run tests**

Run: `npm test -- test/queue/`
Expected: PASS. Notable test: `indexGroupHasMatchInvariant` or similar — ensures the looksFull/invalidate fallback path still triggers when `row.id` is undefined.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/lib/queue/mutators/indexGroup.ts
git commit -m "refactor(queue): indexGroup uses stripQueueMetadata (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task C5: Pin client_op_id absence in cache-shape contract

**Files:**
- Modify: `packages/HimalayaUI/frontend/test/queue/cache-shape.test.ts`

`stripQueueMetadata` strips `client_op_id` from the payload, which is a silent improvement for `peakAdd` and `indexGroup` (both previously left it in `serverPeak` / `row`). Pin the property so a future regression that bypasses `stripQueueMetadata` and writes the raw response into the cache fails this test.

- [ ] **Step 1: Add assertions to the existing cache-shape tests**

In `cache-shape.test.ts`, the file has per-mutator `describe` blocks. Inside the `peakAdd` and `indexGroup` (both confirm + unconfirm) describes, add a test asserting the cache row does NOT carry queue metadata. Pattern:

```ts
it("does not write queue metadata fields into the cache row (peakAdd)", () => {
  const qc = new QueryClient();
  const exposureId = 1;
  const response: PeakAddResponse = {
    id: 7,
    exposure_id: exposureId,
    q: 0.12,
    intensity: null,
    prominence: null,
    sharpness: null,
    source: "manual",
    excluded: false,
    view_row_id: 7,
    event_id: 42,
    analysis_inputs_hash: "abc",
    // @ts-expect-error — client_op_id is plumbing, NOT part of the response type
    client_op_id: "op-1",
  };
  peakAddMutator.onSuccess(
    { exposureId, q: 0.12, kind: "peak_added", payload: { q: 0.12 }, clientOpId: "op-1" } as any,
    response,
    qc,
  );
  const row = qc.getQueryData<Peak[]>(queryKeys.peaks(exposureId))![0];
  expect(row).not.toHaveProperty("event_id");
  expect(row).not.toHaveProperty("view_row_id");
  expect(row).not.toHaveProperty("analysis_inputs_hash");
  expect(row).not.toHaveProperty("client_op_id");
});
```

Add an analogous test for `addIndexToGroupMutator.onSuccess` and `removeIndexFromGroupMutator.onSuccess`. The exact response shape comes from the existing test fixtures in the file — match those.

- [ ] **Step 2: Run tests**

Run: `npm test -- test/queue/cache-shape.test.ts`
Expected: PASS — fields are stripped, assertions hold.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/test/queue/cache-shape.test.ts
git commit -m "test(queue): pin client_op_id absence post-stripQueueMetadata (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Phase D — Verification

### Task D1: Full test sweep + build

- [ ] **Step 1: Run full frontend test suite**

Run: `npm test`
Expected: PASS — entire Vitest suite green.

- [ ] **Step 2: Run typecheck + build**

Run: `npm run build`
Expected: PASS — `tsc --noEmit` clean, Vite build succeeds.

- [ ] **Step 3: Run mocked Playwright E2E (smoke)**

Run: `npm run e2e`
Expected: PASS — all spec files green.

- [ ] **Step 4: Run live Playwright spec for peak_add (one canonical SSE-wins path)**

Prereq: a running backend on the live-test port (see `packages/HimalayaUI/frontend/e2e/live/README.md` for setup). With backend + Vite up:

Run: `npm run e2e:live -- peak-add-no-stale-banner`
Expected: PASS — the SSE-wins path that this refactor touches most heavily is exercised end-to-end (peakAdd's synth + replacePlaceholder + stripQueueMetadata all participate).

If the live-test infrastructure isn't available locally, document the skip in the PR description — the mocked E2E in Step 3 plus the unit-level contract tests cover the same paths at the Layer 1–5 level.

- [ ] **Step 5: Grep for remaining hand-rolled placeholder loops**

Run: `grep -n "id < 0 && !replaced" packages/HimalayaUI/frontend/src/lib/queue/`
Expected: zero matches. (If `peakRemove.ts` or any other file still has a hand-rolled loop, address before closing.)

Run: `grep -nE "event_id: _e, view_row_id: _v" packages/HimalayaUI/frontend/src/lib/queue/`
Expected: zero matches.

- [ ] **Step 6: Grep for the legacy synth switch**

Run: `grep -nE 'remote\.kind === "(peak_added|add_tag|comparison_created|comparison_deleted|peak_excluded|peak_unexcluded)"' packages/HimalayaUI/frontend/src/lib/queue/replayCoordinator.ts`
Expected: zero matches.

- [ ] **Step 7: Final commit if any cleanup landed**

Only if cleanup landed:

```bash
git add -p   # review changes
git commit -m "chore(queue): final cleanup pass (#129)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Acceptance criteria (issue #129)

- [x] `replayCoordinator.synthesizeResponseFromSse` no longer dispatches on a per-kind switch — Task B9 reduces it to ~10 lines using `resolveMutatorForEvent` + fallback.
- [x] Mutators with non-trivial synth own their synth as a mutator method — Tasks B4–B8 migrate the five bespoke kinds.
- [x] Five placeholder-replace loops collapse to one helper call each — Tasks A2–A6 migrate the five sites.
- [x] All queue contract tests green; rollbackSymmetry / cache-shape / sseEventPayload / authHeaders unchanged — each migration task explicitly runs `npm test -- test/queue/` before committing.
