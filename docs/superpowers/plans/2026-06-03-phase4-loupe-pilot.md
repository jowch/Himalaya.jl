# Phase-4 Loupe Pilot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the legacy Loupe page with a greenfield `src/print/pages/LoupePage.tsx` assembled from existing print composites, establishing `src/print/pages/` and the repeatable per-route cutover recipe.

**Architecture:** The greenfield page mounts body-only inside the carried legacy `CorpusShell` `<Outlet>`. It composes three **already-built** print composites — `BigFrame`, `ThumbnailGallery`, `LoupeSidePanel` — plus `PlateHeader`, owns the `activeId` selection state, and wires the CARRIED query hooks (`useCorpusSamples`/`useExposures`/`useSetExposureStatus`/`useSelectExposure`/`useAddCorpusSampleTag`/`useRemoveCorpusSampleTag`). All `Exposure → view-model` mapping lives in a pure, unit-tested adapters module. The legacy `usePeaks`/signal-meter is **dropped** (the greenfield `LoupeSidePanel` has no signal block; the mockup confirms none).

**Tech Stack:** React 18, TypeScript (`exactOptionalPropertyTypes: true`), TanStack Query, react-router-dom, Vitest + RTL, Playwright (mocked e2e), boneyard-js skeletons, Tailwind with the closed-look design-guard (`npm run lint:design`).

**Provenance discipline (from the strategy spec):** every file is NEW (`src/print/**`), OLD (`src/pages/**`, `src/components/**` — deleted in this plan), or CARRIED (`queries.ts`, hooks — imported as-is). `src/print/pages/LoupePage.tsx` must import **only NEW + CARRIED, never OLD**. Verify with a grep before the final gate.

**Standing constraints (do not violate):**
- Commit ONLY specifically-named files. NEVER `git add -A`/`git add .`.
- NEVER stage `src/bones/registry.ts` or `src/bones/contact-sheet.bones.json`.
- This plan file lives under `docs/superpowers/plans/` — NEVER stage it.
- Every commit's exact last line: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
- Typecheck: `npx tsc --noEmit -p tsconfig.build.json`. Work from `packages/HimalayaUI/frontend`.
- Design guard: `src/print/pages/**` is NOT exempt — placement/layout/named-token classes only; no arbitrary `text-[…]`/`rounded-[…]`/raw-colour/side-stripe literals.

---

## File map

- **Create** `src/print/pages/loupeAdapters.ts` — pure `Exposure`/`SampleTag` → view-model helpers (no React).
- **Create** `test/print-pages/loupeAdapters.test.ts` — unit tests for the adapters (new dir `test/print-pages/`).
- **Create** `src/print/pages/LoupePage.tsx` — the greenfield page.
- **Create** `test/print-pages/LoupePage.test.tsx` — RTL tests (render, keyboard, drop/rep, tags, not-found, no-exposures).
- **Modify** `src/components/AppRoutes.tsx` — repoint the `/samples/loupe/:sampleId` import to the greenfield page.
- **Delete** `src/pages/LoupePage.tsx`, `src/components/LoupeFrame.tsx`, `src/components/LoupeSidebar.tsx` (+ any `LoupeTagsEditor` and the legacy tests that target them) — after grep-confirming no other importers.
- **Capture/commit** `src/bones/loupe.bones.json` (commit only this file; never `registry.ts`).
- **Check/update** any mocked loupe spec under `e2e/`.

---

## Reference: verified APIs (do not re-derive)

```ts
// src/api.ts
export interface Exposure {
  id: number; sample_id: number; filename: string | null;
  kind: "file" | "averaged" | "background_subtracted";
  selected: boolean; status: "accepted" | "rejected" | null;
  image_path: string | null; image_version: string;
  tags: ExposureTag[]; sources: unknown[];
  trace_hash: string | null; analysis_inputs_hash: string | null;
}
export interface Sample { id: number; experiment_id: number; name: string | null; display_name: string | null; notes: string | null; tags: SampleTag[]; }
export interface CorpusSample extends Sample { q_units: string; screened?: boolean; phase?: string | null; }
export interface SampleTag { id: number; key: string; value: string; source: string; }

// src/print/ui  (MetaList.tsx, tag.ts)
export interface MetaEntry { key: string; value: ReactNode; }
export interface Tag { key: string; value?: string; }

// src/print/components/ThumbnailGallery.tsx
export interface GalleryExposure { id: number; src: string | null; frameNo?: string | number; representative?: boolean; rejected?: boolean; }
// <ThumbnailGallery exposures selectedId onSelect size="lg" align="center" className />

// src/print/components/BigFrame.tsx
// <BigFrame src={string|null} caption={ReactNode} rejected={boolean} className />

// src/print/components/LoupeSidePanel.tsx
// <LoupeSidePanel meta dropped onToggleDrop isRepresentative onSetRepresentative tags onAddTag onRemoveTag className />

// src/print/components/PlateHeader.tsx
// <PlateHeader title subtitle as="h2" className />   // omit kicker + children for loupe

// src/queries.ts — CARRIED hooks (signatures confirmed against legacy LoupePage)
// useCorpusSamples()                 -> { data?: CorpusSample[], isLoading }
// useExposures(sampleId | undefined) -> { data?: Exposure[], isLoading }
// useSetExposureStatus(sampleId)     -> .mutate({ exposureId, status: "rejected" | null })
// useSelectExposure(sampleId)        -> .mutate(exposureId)
// useAddCorpusSampleTag(sampleId)    -> .mutate({ key, value })
// useRemoveCorpusSampleTag(sampleId) -> .mutate(tagId)

// Image URL (replicate exactly — cache-coherence):
//   full : /api/exposures/{id}/image?v={image_version}
//   thumb: /api/exposures/{id}/image?thumb=1&v={image_version}
//   if image_version is "" → omit ?v ; if image_path is null → src is null (placeholder)
```

---

## Task 1: Pure adapters module

**Files:**
- Create: `src/print/pages/loupeAdapters.ts`
- Test: `test/print-pages/loupeAdapters.test.ts`

- [ ] **Step 1: Write the failing tests**

Create `test/print-pages/loupeAdapters.test.ts`:

```ts
import { describe, it, expect } from "vitest";
import type { Exposure, SampleTag } from "../../src/api";
import {
  defaultExposureId,
  buildExposureImageUrl,
  toGalleryExposures,
  toMetaEntries,
  toLoupeTags,
  findSampleTagId,
} from "../../src/print/pages/loupeAdapters";

function exp(over: Partial<Exposure>): Exposure {
  return {
    id: 1, sample_id: 1, filename: "JC000-001.dat", kind: "file",
    selected: false, status: "accepted", image_path: "/x.tif", image_version: "v9",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null, ...over,
  };
}

describe("defaultExposureId", () => {
  it("prefers the representative (selected)", () => {
    expect(defaultExposureId([exp({ id: 1 }), exp({ id: 2, selected: true })])).toBe(2);
  });
  it("falls back to first accepted", () => {
    expect(defaultExposureId([exp({ id: 1, status: "rejected" }), exp({ id: 2, status: "accepted" })])).toBe(2);
  });
  it("falls back to first exposure", () => {
    expect(defaultExposureId([exp({ id: 7, status: null }), exp({ id: 8, status: null })])).toBe(7);
  });
  it("returns undefined for an empty list", () => {
    expect(defaultExposureId([])).toBeUndefined();
  });
});

describe("buildExposureImageUrl", () => {
  it("builds a full image url with the version param", () => {
    expect(buildExposureImageUrl(exp({ id: 5, image_version: "v9" }))).toBe("/api/exposures/5/image?v=v9");
  });
  it("builds a thumb url", () => {
    expect(buildExposureImageUrl(exp({ id: 5, image_version: "v9" }), { thumb: true })).toBe("/api/exposures/5/image?thumb=1&v=v9");
  });
  it("omits ?v when image_version is empty", () => {
    expect(buildExposureImageUrl(exp({ id: 5, image_version: "" }))).toBe("/api/exposures/5/image");
  });
  it("returns null when there is no image", () => {
    expect(buildExposureImageUrl(exp({ id: 5, image_path: null }))).toBeNull();
  });
});

describe("toGalleryExposures", () => {
  it("maps id/src(thumb)/frameNo/rejected/representative", () => {
    const out = toGalleryExposures([
      exp({ id: 10, selected: true, status: "accepted" }),
      exp({ id: 11, status: "rejected" }),
    ]);
    expect(out[0]).toEqual({ id: 10, src: "/api/exposures/10/image?thumb=1&v=v9", frameNo: 1, representative: true, rejected: false });
    expect(out[1]).toEqual({ id: 11, src: "/api/exposures/11/image?thumb=1&v=v9", frameNo: 2, representative: false, rejected: true });
  });
});

describe("toMetaEntries", () => {
  it("builds frame position + the deferred-field placeholders, no signal row", () => {
    const exposures = [exp({ id: 1 }), exp({ id: 2 }), exp({ id: 3 })];
    expect(toMetaEntries(exposures[1]!, exposures)).toEqual([
      { key: "frame", value: "2 of 3" },
      { key: "integration", value: "—" },
      { key: "collected", value: "—" },
    ]);
  });
});

describe("toLoupeTags / findSampleTagId", () => {
  const tags: SampleTag[] = [
    { id: 100, key: "LL37", value: "", source: "user" },
    { id: 101, key: "temp", value: "37C", source: "user" },
  ];
  it("maps SampleTag -> Tag, omitting empty value (exactOptionalPropertyTypes)", () => {
    expect(toLoupeTags(tags)).toEqual([{ key: "LL37" }, { key: "temp", value: "37C" }]);
  });
  it("finds the tag id by key+value for removal", () => {
    expect(findSampleTagId(tags, { key: "temp", value: "37C" })).toBe(101);
    expect(findSampleTagId(tags, { key: "LL37" })).toBe(100);
    expect(findSampleTagId(tags, { key: "missing" })).toBeUndefined();
  });
});
```

- [ ] **Step 2: Run to verify it fails**

Run: `npm test -- print-pages/loupeAdapters`
Expected: FAIL — cannot resolve `../../src/print/pages/loupeAdapters`.

- [ ] **Step 3: Implement the adapters**

Create `src/print/pages/loupeAdapters.ts`:

```ts
import type { Exposure, SampleTag } from "../../api";
import type { GalleryExposure } from "../components/ThumbnailGallery";
import type { MetaEntry, Tag } from "../ui";

/** Default exposure when the loupe opens: representative → first accepted → first. */
export function defaultExposureId(exposures: Exposure[]): number | undefined {
  const representative = exposures.find((e) => e.selected);
  if (representative) return representative.id;
  const firstAccepted = exposures.find((e) => e.status === "accepted");
  if (firstAccepted) return firstAccepted.id;
  return exposures[0]?.id;
}

/**
 * Detector image URL. Mirrors the legacy DetectorImage builder exactly so the
 * browser cache key matches (cache-coherence): `?thumb=1` for the strip,
 * `?v=<image_version>` only when present. null image_path → null (placeholder).
 */
export function buildExposureImageUrl(
  exposure: Exposure,
  opts?: { thumb?: boolean },
): string | null {
  if (exposure.image_path === null) return null;
  const params = new URLSearchParams();
  if (opts?.thumb) params.set("thumb", "1");
  if (exposure.image_version) params.set("v", exposure.image_version);
  const qs = params.toString();
  return `/api/exposures/${exposure.id}/image${qs ? `?${qs}` : ""}`;
}

/** Map the per-sample exposure list to the filmstrip view-model. */
export function toGalleryExposures(exposures: Exposure[]): GalleryExposure[] {
  return exposures.map((e, i) => ({
    id: e.id,
    src: buildExposureImageUrl(e, { thumb: true }),
    frameNo: i + 1,
    representative: e.selected,
    rejected: e.status === "rejected",
  }));
}

/**
 * "This exposure" metadata rows. `integration`/`collected` are placeholders
 * until the backend lands those fields (#256). The legacy signal-meter row is
 * intentionally dropped — the greenfield LoupeSidePanel has no signal block.
 */
export function toMetaEntries(active: Exposure, exposures: Exposure[]): MetaEntry[] {
  const idx = exposures.findIndex((e) => e.id === active.id);
  const position = idx >= 0 ? `${idx + 1} of ${exposures.length}` : "—";
  return [
    { key: "frame", value: position },
    { key: "integration", value: "—" },
    { key: "collected", value: "—" },
  ];
}

/** SampleTag[] → greenfield Tag[]; omit empty value (exactOptionalPropertyTypes). */
export function toLoupeTags(tags: SampleTag[]): Tag[] {
  return tags.map((t) => (t.value ? { key: t.key, value: t.value } : { key: t.key }));
}

/** Resolve a greenfield Tag back to its SampleTag id for removal (key+value match). */
export function findSampleTagId(tags: SampleTag[], tag: Tag): number | undefined {
  const want = tag.value ?? "";
  return tags.find((t) => t.key === tag.key && (t.value ?? "") === want)?.id;
}
```

- [ ] **Step 4: Run to verify it passes**

Run: `npm test -- print-pages/loupeAdapters`
Expected: PASS (all cases).

- [ ] **Step 5: Typecheck + design guard**

Run: `npx tsc --noEmit -p tsconfig.build.json` → no errors.
Run: `npm run lint:design` → exit 0 (the adapters file has no JSX/classes, so it is design-guard-neutral).

- [ ] **Step 6: Commit**

```bash
git add src/print/pages/loupeAdapters.ts test/print-pages/loupeAdapters.test.ts
git commit -m "$(printf 'Phase-4 Loupe: pure Exposure->view-model adapters\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 2: The greenfield LoupePage component

**Files:**
- Create: `src/print/pages/LoupePage.tsx`
- Test: `test/print-pages/LoupePage.test.tsx`

- [ ] **Step 1: Write the failing tests**

Create `test/print-pages/LoupePage.test.tsx`. Mock the CARRIED query hooks so the test is deterministic (the queue/network never runs):

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import type { Exposure, CorpusSample } from "../../src/api";

const setStatusMutate = vi.fn();
const selectMutate = vi.fn();
const addTagMutate = vi.fn();
const removeTagMutate = vi.fn();

const state = {
  samples: [] as CorpusSample[],
  exposures: [] as Exposure[],
  loading: false,
};

vi.mock("../../src/queries", () => ({
  useCorpusSamples: () => ({ data: state.samples, isLoading: state.loading }),
  useExposures: () => ({ data: state.exposures, isLoading: state.loading }),
  useSetExposureStatus: () => ({ mutate: setStatusMutate }),
  useSelectExposure: () => ({ mutate: selectMutate }),
  useAddCorpusSampleTag: () => ({ mutate: addTagMutate }),
  useRemoveCorpusSampleTag: () => ({ mutate: removeTagMutate }),
}));

// boneyard Skeleton: render children when not loading (avoid capture machinery in JSDOM).
vi.mock("boneyard-js/react", () => ({
  Skeleton: ({ children }: { children: React.ReactNode }) => <>{children}</>,
}));

import { LoupePage } from "../../src/print/pages/LoupePage";

function exp(over: Partial<Exposure>): Exposure {
  return {
    id: 1, sample_id: 1, filename: "JC000-001.dat", kind: "file",
    selected: false, status: "accepted", image_path: "/x.tif", image_version: "v9",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null, ...over,
  };
}

function renderAt(sampleId: number) {
  return render(
    <MemoryRouter initialEntries={[`/samples/loupe/${sampleId}`]}>
      <Routes>
        <Route path="/samples/loupe/:sampleId" element={<LoupePage />} />
        <Route path="/samples" element={<div data-testid="sheet">sheet</div>} />
      </Routes>
    </MemoryRouter>,
  );
}

beforeEach(() => {
  vi.clearAllMocks();
  state.samples = [{
    id: 42, experiment_id: 1, name: "JC042", display_name: "JC042 — LL37",
    notes: null, q_units: "A-1",
    tags: [{ id: 100, key: "LL37", value: "", source: "user" }],
  }];
  state.exposures = [
    exp({ id: 1, selected: true }),
    exp({ id: 2, status: "rejected" }),
  ];
  state.loading = false;
});

describe("LoupePage", () => {
  it("renders the headline, frame, side panel and filmstrip", () => {
    renderAt(42);
    expect(screen.getByTestId("loupe-page")).toBeInTheDocument();
    expect(screen.getByRole("heading", { name: /JC042 — LL37/ })).toBeInTheDocument();
    expect(screen.getByTestId("big-frame")).toBeInTheDocument();
    expect(screen.getByTestId("loupe-side-panel")).toBeInTheDocument();
  });

  it("opens on the representative exposure (id 1)", () => {
    renderAt(42);
    // big-frame is not rejected because exposure 1 is accepted+selected
    expect(screen.getByTestId("big-frame")).not.toHaveAttribute("data-rejected");
  });

  it("drop toggle mutates status to rejected for the active exposure", () => {
    renderAt(42);
    fireEvent.keyDown(window, { key: "x" });
    expect(setStatusMutate).toHaveBeenCalledWith({ exposureId: 1, status: "rejected" });
  });

  it("R sets the representative", () => {
    renderAt(42);
    fireEvent.keyDown(window, { key: "r" });
    expect(selectMutate).toHaveBeenCalledWith(1);
  });

  it("ArrowRight flips to the next exposure", () => {
    renderAt(42);
    fireEvent.keyDown(window, { key: "ArrowRight" });
    // exposure 2 is rejected → big-frame shows the dropped state
    expect(screen.getByTestId("big-frame")).toHaveAttribute("data-rejected", "true");
  });

  it("Escape navigates back to the sheet", () => {
    renderAt(42);
    fireEvent.keyDown(window, { key: "Escape" });
    expect(screen.getByTestId("sheet")).toBeInTheDocument();
  });

  it("Back button navigates to the sheet", () => {
    renderAt(42);
    fireEvent.click(screen.getByTestId("loupe-back"));
    expect(screen.getByTestId("sheet")).toBeInTheDocument();
  });

  it("removing a sample tag resolves its id via key+value", () => {
    renderAt(42);
    const panel = screen.getByTestId("loupe-side-panel");
    // TagList renders a remove control per tag; find the LL37 tag's remove button.
    const removeBtn = within(panel).getByRole("button", { name: /remove .*LL37/i });
    fireEvent.click(removeBtn);
    expect(removeTagMutate).toHaveBeenCalledWith(100);
  });

  it("shows not-found when the sample is missing", () => {
    state.samples = [];
    renderAt(999);
    expect(screen.getByTestId("loupe-not-found")).toBeInTheDocument();
  });

  it("shows the no-exposures state", () => {
    state.exposures = [];
    renderAt(42);
    expect(screen.getByText(/no exposures/i)).toBeInTheDocument();
  });
});
```

> Note for the implementer: the exact accessible name of TagList's remove control (`/remove .*LL37/i`) must match what `src/print/ui` `TagList` actually renders — open `TagList.tsx`, confirm the `aria-label`, and adjust the matcher. If TagList exposes a different remove affordance, query it the way the component actually labels it (do not assert on class names).

- [ ] **Step 2: Run to verify it fails**

Run: `npm test -- print-pages/LoupePage`
Expected: FAIL — cannot resolve `../../src/print/pages/LoupePage`.

- [ ] **Step 3: Implement the page**

Create `src/print/pages/LoupePage.tsx`:

```tsx
import { useCallback, useEffect, useMemo, useState } from "react";
import { useNavigate, useParams, useSearchParams } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import {
  useCorpusSamples,
  useExposures,
  useSetExposureStatus,
  useSelectExposure,
  useAddCorpusSampleTag,
  useRemoveCorpusSampleTag,
} from "../../queries";
import type { Tag } from "../ui";
import { BigFrame } from "../components/BigFrame";
import { ThumbnailGallery } from "../components/ThumbnailGallery";
import { LoupeSidePanel } from "../components/LoupeSidePanel";
import { PlateHeader } from "../components/PlateHeader";
import {
  defaultExposureId,
  buildExposureImageUrl,
  toGalleryExposures,
  toMetaEntries,
  toLoupeTags,
  findSampleTagId,
} from "./loupeAdapters";

// Boneyard fixture — a real render with mock props so the headless capture CLI
// measures the greenfield loupe body. image_path:null → DetectorImage takes its
// placeholder branch (a plain rectangle at the frame's fixed aspect).
const FIXTURE_EXPOSURE = {
  id: 0, sample_id: 0, filename: "JC000-001.dat", kind: "file" as const,
  selected: false, status: "accepted" as const, image_path: null,
  image_version: "", tags: [], sources: [],
  trace_hash: null, analysis_inputs_hash: null,
};
const LOUPE_FIXTURE = (
  <div className="grid grid-cols-[1fr_286px] gap-7">
    <div>
      <BigFrame src={null} caption="frame 1 of 1 · kept" />
      <ThumbnailGallery
        exposures={[{ id: 0, src: null, frameNo: 1 }]}
        selectedId={0}
        size="lg"
        align="center"
        className="mt-3"
      />
    </div>
    <LoupeSidePanel
      meta={toMetaEntries(FIXTURE_EXPOSURE, [FIXTURE_EXPOSURE])}
      dropped={false}
      isRepresentative={false}
      tags={[]}
    />
  </div>
);

/**
 * LoupePage (greenfield) — the sample loupe at /samples/loupe/:sampleId.
 * URL-owned: the sample id is the route param, never Zustand `activeSampleId`.
 * Mounts body-only inside the carried CorpusShell <Outlet>.
 */
export function LoupePage(): JSX.Element {
  const { sampleId: sampleIdParam } = useParams<{ sampleId: string }>();
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();
  const sampleId = Number(sampleIdParam);
  const hasValidId = Number.isFinite(sampleId);

  const corpusQ = useCorpusSamples();
  const exposuresQ = useExposures(hasValidId ? sampleId : undefined);

  const sample = corpusQ.data?.find((s) => s.id === sampleId);
  const exposures = useMemo(() => exposuresQ.data ?? [], [exposuresQ.data]);
  const isLoading = corpusQ.isLoading || exposuresQ.isLoading;

  // Active exposure — local state, defaulted by defaultExposureId, reset on
  // sample change so the next sample picks its own default.
  const [activeId, setActiveId] = useState<number | undefined>(undefined);
  useEffect(() => { setActiveId(undefined); }, [sampleId]);
  const computedDefault = defaultExposureId(exposures);
  useEffect(() => {
    if (activeId === undefined && computedDefault !== undefined) setActiveId(computedDefault);
  }, [activeId, computedDefault]);
  const activeExposure = exposures.find((e) => e.id === activeId);

  const frameIndex = exposures.findIndex((e) => e.id === activeId);
  const exposurePosition =
    frameIndex >= 0 ? `exposure ${frameIndex + 1} of ${exposures.length}` : "—";

  const setStatus = useSetExposureStatus(hasValidId ? sampleId : 0);
  const setRepresentative = useSelectExposure(hasValidId ? sampleId : 0);
  const addTag = useAddCorpusSampleTag(hasValidId ? sampleId : 0);
  const removeTag = useRemoveCorpusSampleTag(hasValidId ? sampleId : 0);

  const handleDropToggle = useCallback(() => {
    if (!activeExposure) return;
    setStatus.mutate({
      exposureId: activeExposure.id,
      status: activeExposure.status === "rejected" ? null : "rejected",
    });
  }, [activeExposure, setStatus]);

  const handleSetRepresentative = useCallback(() => {
    if (!activeExposure) return;
    setRepresentative.mutate(activeExposure.id);
  }, [activeExposure, setRepresentative]);

  const handleAddTag = useCallback((t: Tag) => {
    addTag.mutate({ key: t.key, value: t.value ?? "" });
  }, [addTag]);

  const handleRemoveTag = useCallback((t: Tag) => {
    if (!sample) return;
    const id = findSampleTagId(sample.tags, t);
    if (id !== undefined) removeTag.mutate(id); // optimistic-add w/o id → no-op (ledger risk)
  }, [removeTag, sample]);

  const flip = useCallback((delta: number) => {
    if (activeId === undefined || exposures.length === 0) return;
    const idx = exposures.findIndex((e) => e.id === activeId);
    if (idx < 0) return;
    const next = Math.min(Math.max(idx + delta, 0), exposures.length - 1);
    setActiveId(exposures[next]!.id);
  }, [activeId, exposures]);

  const goBack = useCallback(() => {
    const beamtime = searchParams.get("beamtime");
    navigate(beamtime ? `/samples?beamtime=${beamtime}` : "/samples");
  }, [navigate, searchParams]);

  useEffect(() => {
    function onKeyDown(e: KeyboardEvent): void {
      const tag = (e.target as HTMLElement | null)?.tagName;
      if (tag === "INPUT" || tag === "TEXTAREA") return;
      if (e.key === "ArrowLeft") flip(-1);
      else if (e.key === "ArrowRight") flip(1);
      else if (e.key === "x" || e.key === "X") handleDropToggle();
      else if (e.key === "r" || e.key === "R") handleSetRepresentative();
      else if (e.key === "Escape") goBack();
    }
    window.addEventListener("keydown", onKeyDown);
    return () => window.removeEventListener("keydown", onKeyDown);
  }, [flip, handleDropToggle, handleSetRepresentative, goBack]);

  if (!corpusQ.isLoading && !sample) {
    return (
      <div data-testid="loupe-page" className="mx-auto max-w-[1100px] px-8 py-7">
        <div data-testid="loupe-not-found" className="rounded border border-hair-strong p-8 text-sm text-ink-faint">
          Sample not found.{" "}
          <button onClick={goBack} className="font-semibold text-print-accent hover:underline">
            Back to the sheet
          </button>
        </div>
      </div>
    );
  }

  const isDropped = activeExposure?.status === "rejected";

  return (
    <div data-testid="loupe-page" className="mx-auto max-w-[1100px] px-8 py-7">
      <button data-testid="loupe-back" onClick={goBack} className="mb-3.5 text-sm font-semibold text-print-accent hover:underline">
        ← Back to the sheet
      </button>
      <PlateHeader
        title={sample?.display_name ?? sample?.name ?? "—"}
        subtitle={`${sample?.name ?? "—"} · ${exposurePosition}`}
        className="mb-5"
      />
      <Skeleton name="loupe" className="block" loading={isLoading} stagger={50} transition={200}
        fixture={LOUPE_FIXTURE}
        fallback={<div data-testid="loupe-skeleton" className="p-8 text-sm italic text-ink-faint">Loading sample…</div>}>
        <div className="grid grid-cols-[1fr_286px] gap-7">
          {sample && activeExposure ? (
            <>
              <div>
                <BigFrame
                  src={buildExposureImageUrl(activeExposure)}
                  caption={`frame ${frameIndex + 1} of ${exposures.length} · ${isDropped ? "dropped" : "kept"}`}
                  rejected={isDropped}
                />
                <ThumbnailGallery
                  exposures={toGalleryExposures(exposures)}
                  selectedId={activeId}
                  onSelect={setActiveId}
                  size="lg"
                  align="center"
                  className="mt-3"
                />
              </div>
              <LoupeSidePanel
                meta={toMetaEntries(activeExposure, exposures)}
                dropped={!!isDropped}
                isRepresentative={activeExposure.selected}
                tags={toLoupeTags(sample.tags)}
                onToggleDrop={handleDropToggle}
                onSetRepresentative={handleSetRepresentative}
                onAddTag={handleAddTag}
                onRemoveTag={handleRemoveTag}
              />
            </>
          ) : (
            <div className="col-span-2 p-8 text-sm text-ink-faint">This sample has no exposures.</div>
          )}
        </div>
      </Skeleton>
    </div>
  );
}
```

- [ ] **Step 4: Run to verify it passes**

Run: `npm test -- print-pages/LoupePage`
Expected: PASS. If the TagList remove-name matcher fails, open `src/print/ui/TagList.tsx`, read the actual remove `aria-label`, and fix the matcher (test bug, not code bug).

- [ ] **Step 5: Typecheck + design guard**

Run: `npx tsc --noEmit -p tsconfig.build.json` → no errors (watch `exactOptionalPropertyTypes`: the conditional handlers are always passed here, which is fine — `LoupeSidePanel` already guards its optional callbacks).
Run: `npm run lint:design` → exit 0. If it flags `text-print-accent` on the back button (named token, should pass), extract a `Button` `link`/`text` variant in `src/print/ui/` (design-guard-exempt) as a refactor-on-contact and use it instead — keep the page placement-only.

- [ ] **Step 6: Commit**

```bash
git add src/print/pages/LoupePage.tsx test/print-pages/LoupePage.test.tsx
git commit -m "$(printf 'Phase-4 Loupe: greenfield LoupePage assembled from print composites\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 3: Repoint the route

**Files:**
- Modify: `src/components/AppRoutes.tsx`

- [ ] **Step 1: Confirm the current import + route**

Run: `grep -n "LoupePage" src/components/AppRoutes.tsx`
Expected: an import `import { LoupePage } from "../pages/LoupePage";` and `<Route path="/samples/loupe/:sampleId" element={<LoupePage />} />`.

- [ ] **Step 2: Repoint the import (atomic swap — the route element name stays the same)**

Edit `src/components/AppRoutes.tsx`: change the import path only.

```diff
-import { LoupePage } from "../pages/LoupePage";
+import { LoupePage } from "../print/pages/LoupePage";
```

Leave the `<Route ...>` line unchanged.

- [ ] **Step 3: Typecheck + the mocked e2e**

Run: `npx tsc --noEmit -p tsconfig.build.json` → no errors.
Run: `grep -rln "loupe" e2e --include=*.ts` to find any mocked loupe spec. If one exists, run it: `npm run e2e -- <loupe-spec>`. Update selectors only if the test asserts legacy testids that the greenfield page renamed (the page keeps `loupe-page` + `loupe-back`; the headline is now a `PlateHeader` heading — query by `getByRole("heading")` instead of a legacy `loupe-headline` testid). Do NOT weaken assertions to pass; fix them to the new structure.

- [ ] **Step 4: Commit**

```bash
git add src/components/AppRoutes.tsx
# add any updated e2e spec by exact path if you changed one
git commit -m "$(printf 'Phase-4 Loupe: route /samples/loupe/:sampleId -> greenfield page\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 4: Delete the OLD

**Files:**
- Delete: `src/pages/LoupePage.tsx`, `src/components/LoupeFrame.tsx`, `src/components/LoupeSidebar.tsx`, any `src/components/LoupeTagsEditor.tsx`, and the legacy tests that target them.

- [ ] **Step 1: Confirm nothing else imports them**

Run:
```bash
grep -rn "pages/LoupePage\|components/LoupeFrame\|components/LoupeSidebar\|LoupeTagsEditor" src test e2e
```
Expected: the ONLY remaining references are the files themselves + their own tests. If anything else imports them, STOP and report — a non-loupe surface depends on them (unexpected; investigate before deleting).

- [ ] **Step 2: Identify the legacy tests**

Run: `ls test | grep -i loupe ; grep -rln "LoupeFrame\|LoupeSidebar\|LoupeTagsEditor" test`
Note every legacy loupe test file (component-level tests for the deleted components). Page-level legacy tests that asserted the OLD page body are superseded by `test/print-pages/LoupePage.test.tsx`.

- [ ] **Step 3: Delete the files**

```bash
git rm src/pages/LoupePage.tsx src/components/LoupeFrame.tsx src/components/LoupeSidebar.tsx
# git rm src/components/LoupeTagsEditor.tsx   # only if it exists
# git rm <each legacy loupe test file identified in Step 2 by exact path>
```

- [ ] **Step 4: Verify the build is whole**

Run: `npx tsc --noEmit -p tsconfig.build.json` → no errors (no dangling imports).
Run: `npm test -- loupe` → only the new `print-pages` loupe tests run and pass; no references to deleted modules remain.

- [ ] **Step 5: Commit**

```bash
# stage exactly the deletions listed above (no -A)
git commit -m "$(printf 'Phase-4 Loupe: delete legacy LoupePage + LoupeFrame + LoupeSidebar\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 5: Boneyard skeleton recapture

**Files:**
- Capture/commit: `src/bones/loupe.bones.json` (ONLY this file).

> This step is semi-manual (the boneyard Vite plugin captures during `npm run dev`). It also regenerates `src/bones/registry.ts`, which the standing constraints forbid committing. The implementer should do the capture, then stage ONLY `loupe.bones.json`. If wiring the new bone truly requires a `registry.ts` change, STOP and surface it to the human (do not commit `registry.ts`).

- [ ] **Step 1: Capture**

Run `npm run dev`, open `/samples/loupe/<a real sample id>` against a backend (or the mocked dev data), and trigger the cold-loading state so the boneyard plugin measures the greenfield body and writes `src/bones/loupe.bones.json`. Confirm the file appeared: `ls -l src/bones/loupe.bones.json`.

- [ ] **Step 2: Verify geometry**

Reload with the skeleton visible; confirm the captured bones match the new 2-col `1fr / 286px` layout (big frame + strip on the left, side panel on the right) rather than the legacy layout.

- [ ] **Step 3: Commit ONLY the bone**

```bash
git add src/bones/loupe.bones.json
git status --short src/bones/   # registry.ts may show as modified — do NOT stage it
git commit -m "$(printf 'Phase-4 Loupe: capture greenfield loupe skeleton bones\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 6: Final gate + visual fidelity

- [ ] **Step 1: Provenance grep (NEW imports only)**

Run:
```bash
grep -rn "from \"\.\./\.\./components\|from \"\.\./\.\./pages\|/src/components/\|/src/pages/" src/print/pages/
```
Expected: **no matches** — the greenfield page imports only NEW (`../components`, `../ui`) + CARRIED (`../../queries`, `../../api`). Any hit is an OLD import that violates the boundary rule — fix it.

- [ ] **Step 2: Full local gate**

Run each; all must pass:
```bash
npm run lint:design                                   # exit 0
npx tsc --noEmit -p tsconfig.build.json               # no errors
npm test -- print-pages                               # adapters + page green
npm run e2e -- <loupe-spec if present>                # mocked e2e green
npm run build                                         # exit 0
```

- [ ] **Step 3: Manual visual fidelity check**

Run `npm run dev`, open the loupe, and compare side-by-side against the mockup `docs/redesign-mockups/sample-table.html` (the `.loupe-shell` view — open it in a browser and toggle the loupe view). Confirm: back link, serif headline + mono subtitle, 2-col body, big frame with frame caption + (for a rejected exposure) the "Dropped" pill + grease-pencil ✕, the centered filmstrip with rep/rejected markers, and the side panel's five blocks (This exposure / Verdict / Representative / Sample tags / Keys). Note any visible divergence; fix placement-only issues in the page, appearance issues in the relevant `ui/` primitive.

- [ ] **Step 4: Done**

No extra commit unless Step 3 required a fix. The pilot is complete: `/samples/loupe/:sampleId` renders the greenfield body under the carried shell, the OLD is deleted, and `src/print/pages/` + the per-route recipe are established for Focus next.

---

## Self-review checklist (run before declaring done)

- **Spec coverage:** assemble (T2) · wire CARRIED (T2) · adapter for URL/meta/tags (T1) · route swap (T3) · delete OLD (T4) · boneyard recapture (T5) · gate incl. visual fidelity (T6) · provenance grep (T6.1). The `usePeaks`/signal-meter drop is realized by simply not wiring it (T2). ✅
- **No placeholders:** every code step has complete code; commands have expected output. ✅
- **Type consistency:** `buildExposureImageUrl`, `toGalleryExposures`, `toMetaEntries`, `toLoupeTags`, `findSampleTagId`, `defaultExposureId` names match between T1 implementation, its tests, and T2's imports. Hook names + `.mutate` shapes match the verified `queries.ts` signatures. ✅
- **Provenance:** the page imports only NEW + CARRIED; T6.1 enforces it by grep. ✅
```
