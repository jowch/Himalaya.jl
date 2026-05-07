/**
 * Compare draft state + sessionStorage persistence (Plan §Phase 4, Task 4.3).
 *
 * Asserts:
 * - Round-trip persistence: write actions mutate state and the
 *   sessionStorage mirror.
 * - Schema-version mismatch on reload → ignored cleanly (draft stays null).
 * - Each named action behaves as documented.
 * - `loadDraftFromComparison` recompute-snapshot recovery path: snapshots
 *   are derived from the current TanStack cache, NOT the saved snapshot.
 */
import { describe, it, expect, beforeEach } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import {
  loadDraftFromSession,
  COMPARE_DRAFT_KEY,
  COMPARE_DRAFT_VERSION,
  type ActiveDraft,
  type DraftMember,
} from "../src/lib/comparison/draft";
import { useAppState } from "../src/state";
import { queryKeys } from "../src/queries";
import type {
  Peak, IndexEntry, GroupEntry, Exposure, Comparison, ComparisonMember,
} from "../src/api";

function makeMember(
  overrides: Partial<DraftMember> = {},
): DraftMember {
  return {
    id: undefined,
    exposure_id: 100,
    display_order: 0,
    band_height: 1.0,
    y_offset: 0.0,
    normalization: "none",
    color_override: undefined,
    label_override: undefined,
    q_window_min: undefined,
    q_window_max: undefined,
    peak_display: undefined,
    snapshot: undefined,
    ...overrides,
  };
}

function resetStore(): void {
  // Restore the persisted slice to defaults; clear ephemeral draft.
  useAppState.setState({
    activeDraft: null,
  });
}

describe("compare draft state + sessionStorage", () => {
  beforeEach(() => {
    sessionStorage.clear();
    resetStore();
  });

  it("startNewDraft initializes an empty draft with id and baseHash undefined", () => {
    useAppState.getState().startNewDraft();
    const d = useAppState.getState().activeDraft;
    expect(d).not.toBeNull();
    expect(d!.id).toBeUndefined();
    expect(d!.baseHash).toBeUndefined();
    expect(d!.title).toBe("");
    expect(d!.description).toBe("");
    expect(d!.members).toEqual([]);
  });

  it("setDraftTitle and setDraftDescription mutate fields", () => {
    useAppState.getState().startNewDraft();
    useAppState.getState().setDraftTitle("My comparison");
    useAppState.getState().setDraftDescription("hello");
    const d = useAppState.getState().activeDraft!;
    expect(d.title).toBe("My comparison");
    expect(d.description).toBe("hello");
  });

  it("addMember appends with the next display_order; removeMember removes by index", () => {
    const qc = new QueryClient();
    useAppState.getState().startNewDraft();
    useAppState.getState().addMember(101, qc);
    useAppState.getState().addMember(102, qc);
    let members = useAppState.getState().activeDraft!.members;
    expect(members).toHaveLength(2);
    expect(members[0]!.exposure_id).toBe(101);
    expect(members[0]!.display_order).toBe(0);
    expect(members[1]!.exposure_id).toBe(102);
    expect(members[1]!.display_order).toBe(1);

    useAppState.getState().removeMember(0);
    members = useAppState.getState().activeDraft!.members;
    expect(members).toHaveLength(1);
    // After removal, display_order is renumbered.
    expect(members[0]!.exposure_id).toBe(102);
    expect(members[0]!.display_order).toBe(0);
  });

  it("updateMember partially mutates the indexed member", () => {
    const qc = new QueryClient();
    useAppState.getState().startNewDraft();
    useAppState.getState().addMember(101, qc);
    useAppState.getState().updateMember(0, { y_offset: 0.5, normalization: "max" });
    const m = useAppState.getState().activeDraft!.members[0]!;
    expect(m.y_offset).toBe(0.5);
    expect(m.normalization).toBe("max");
    expect(m.exposure_id).toBe(101);
  });

  it("reorderMembers reorders by index list and renumbers display_order", () => {
    const qc = new QueryClient();
    useAppState.getState().startNewDraft();
    useAppState.getState().addMember(101, qc);
    useAppState.getState().addMember(102, qc);
    useAppState.getState().addMember(103, qc);
    useAppState.getState().reorderMembers([2, 0, 1]);
    const ms = useAppState.getState().activeDraft!.members;
    expect(ms.map((m) => m.exposure_id)).toEqual([103, 101, 102]);
    expect(ms.map((m) => m.display_order)).toEqual([0, 1, 2]);
  });

  it("resizeBands shifts a band's height and clamps the neighbor", () => {
    const qc = new QueryClient();
    useAppState.getState().startNewDraft();
    useAppState.getState().addMember(101, qc);
    useAppState.getState().addMember(102, qc);
    // memberIdx=0 grows by 20px out of 200px total → +0.10 ratio; neighbor shrinks.
    useAppState.getState().resizeBands(0, 20, 200);
    const ms = useAppState.getState().activeDraft!.members;
    expect(ms[0]!.band_height).toBeCloseTo(1.10, 5);
    expect(ms[1]!.band_height).toBeCloseTo(0.90, 5);
  });

  it("discardDraft clears the slot and the sessionStorage mirror", () => {
    const qc = new QueryClient();
    useAppState.getState().startNewDraft();
    useAppState.getState().addMember(101, qc);
    expect(sessionStorage.getItem(COMPARE_DRAFT_KEY)).not.toBeNull();
    useAppState.getState().discardDraft();
    expect(useAppState.getState().activeDraft).toBeNull();
    expect(sessionStorage.getItem(COMPARE_DRAFT_KEY)).toBeNull();
  });

  it("round-trips through sessionStorage with the current schema version", () => {
    useAppState.getState().startNewDraft();
    useAppState.getState().setDraftTitle("title");
    const raw = sessionStorage.getItem(COMPARE_DRAFT_KEY);
    expect(raw).not.toBeNull();
    const parsed = JSON.parse(raw!);
    expect(parsed.version).toBe(COMPARE_DRAFT_VERSION);
    expect(parsed.draft.title).toBe("title");

    const loaded = loadDraftFromSession();
    expect(loaded).not.toBeNull();
    expect(loaded!.title).toBe("title");
  });

  it("loadDraftFromSession ignores an unknown schema version", () => {
    sessionStorage.setItem(
      COMPARE_DRAFT_KEY,
      JSON.stringify({
        version: 99999, // future incompatible version
        draft: { id: undefined, baseHash: undefined, title: "stale", description: "", members: [] },
      }),
    );
    const loaded = loadDraftFromSession();
    expect(loaded).toBeNull();
  });

  it("loadDraftFromSession ignores malformed JSON", () => {
    sessionStorage.setItem(COMPARE_DRAFT_KEY, "{not json");
    expect(loadDraftFromSession()).toBeNull();
  });

  it("startNewDraft is a no-op when the active draft is already a fresh one (id undefined)", () => {
    // Push the URL→draft hydration guards into the action body so the
    // ComparePageEdit hydration effect can drop `draft` from its read set
    // (and from its deps) without an eslint-disable. See the effect at
    // ComparePageEdit.tsx hydration step.
    useAppState.getState().startNewDraft();
    const before = useAppState.getState().activeDraft;
    expect(before).not.toBeNull();
    // Mutate the draft to give us a way to detect a re-create (re-create
    // would clobber the title back to "").
    useAppState.getState().setDraftTitle("in-progress");
    expect(useAppState.getState().activeDraft!.title).toBe("in-progress");

    // Re-call startNewDraft — should NOT clobber the in-progress draft.
    useAppState.getState().startNewDraft();
    expect(useAppState.getState().activeDraft!.title).toBe("in-progress");
  });

  it("startNewDraft DOES replace a non-new draft (one with a comparison id)", () => {
    // Seed a draft tied to a comparison id (mimic post-loadDraftFromComparison).
    useAppState.setState({
      activeDraft: {
        id: 42,
        baseHash: "h",
        title: "loaded",
        description: "",
        members: [],
        forkedFromId: undefined,
        forkedAtHash: undefined,
      },
    });
    useAppState.getState().startNewDraft();
    const after = useAppState.getState().activeDraft!;
    expect(after.id).toBeUndefined();
    expect(after.title).toBe("");
  });

  it("loadDraftFromComparison is a no-op when the active draft already matches the comparison id", () => {
    const qc = new QueryClient();
    // Seed a draft ALREADY tied to comparison 7 + mutate title in-progress.
    useAppState.setState({
      activeDraft: {
        id: 7,
        baseHash: "h-server",
        title: "in-progress edit",
        description: "",
        members: [],
        forkedFromId: undefined,
        forkedAtHash: undefined,
      },
    });

    const comparison: Comparison = {
      id: 7,
      title: "Stored",
      description: "from server",
      content_hash: "h-server",
      created_by: 1,
      created_at: "2026-04-01T00:00:00Z",
      updated_at: "2026-04-15T00:00:00Z",
      forked_from_id: null,
      forked_at_hash: null,
      forked_from_title: null,
      members: [],
    };

    useAppState.getState().loadDraftFromComparison(comparison, qc);
    // No clobber: in-progress title preserved.
    expect(useAppState.getState().activeDraft!.title).toBe("in-progress edit");
  });

  it("loadDraftFromComparison DOES replace a draft tied to a different comparison id", () => {
    const qc = new QueryClient();
    useAppState.setState({
      activeDraft: {
        id: 99,
        baseHash: "h-other",
        title: "different",
        description: "",
        members: [],
        forkedFromId: undefined,
        forkedAtHash: undefined,
      },
    });

    const comparison: Comparison = {
      id: 7,
      title: "Stored",
      description: "from server",
      content_hash: "h-server",
      created_by: 1,
      created_at: "2026-04-01T00:00:00Z",
      updated_at: "2026-04-15T00:00:00Z",
      forked_from_id: null,
      forked_at_hash: null,
      forked_from_title: null,
      members: [],
    };

    useAppState.getState().loadDraftFromComparison(comparison, qc);
    expect(useAppState.getState().activeDraft!.id).toBe(7);
    expect(useAppState.getState().activeDraft!.title).toBe("Stored");
  });

  it("loadDraftFromComparison DOES replace a null draft", () => {
    const qc = new QueryClient();
    expect(useAppState.getState().activeDraft).toBeNull();
    const comparison: Comparison = {
      id: 7,
      title: "Stored",
      description: "from server",
      content_hash: "h-server",
      created_by: 1,
      created_at: "2026-04-01T00:00:00Z",
      updated_at: "2026-04-15T00:00:00Z",
      forked_from_id: null,
      forked_at_hash: null,
      forked_from_title: null,
      members: [],
    };
    useAppState.getState().loadDraftFromComparison(comparison, qc);
    expect(useAppState.getState().activeDraft!.id).toBe(7);
  });

  it("loadDraftFromComparison computes fresh snapshots against the cache", () => {
    const qc = new QueryClient();
    const exposureId = 200;
    // Seed cache with peaks/indices/groups/exposure that determine the snapshot.
    const peaks: Peak[] = [
      { id: 1, exposure_id: exposureId, q: 0.10, intensity: 1.0, sharpness: 0.5, source: "auto", excluded: false },
      { id: 2, exposure_id: exposureId, q: 0.20, intensity: 0.8, sharpness: 0.3, source: "auto", excluded: false },
    ];
    const indices: IndexEntry[] = [];
    const groups: GroupEntry[] = [];
    const exposure: Exposure = {
      id: exposureId,
      sample_id: 1,
      filename: "x.dat",
      kind: "file",
      selected: false,
      status: "accepted",
      image_path: null,
      image_version: "",
      tags: [],
      sources: [],
      trace_hash: "tr",
      analysis_inputs_hash: "abcd",
    };
    qc.setQueryData(queryKeys.peaks(exposureId), peaks);
    qc.setQueryData(queryKeys.indices(exposureId), indices);
    qc.setQueryData(queryKeys.groups(exposureId), groups);
    qc.setQueryData(queryKeys.exposure(exposureId), exposure);

    const savedMember: ComparisonMember = {
      id: 999,
      comparison_id: 7,
      exposure_id: exposureId,
      display_order: 0,
      band_height: 1.5,
      y_offset: 0.2,
      normalization: "max",
      color_override: "#ff0000",
      label_override: "Lipid A",
      q_window_min: 0.05,
      q_window_max: 0.30,
      peak_display: { hidden: [], labeled: [] },
      // The saved snapshot is intentionally STALE so we can verify the
      // recovery path picks the cache shape rather than echoing this.
      snapshot: {
        effective_peaks: [
          { id: 999, q: 9.99, intensity: 9.99, sharpness: 9.99, source: "auto" },
        ],
        confirmed_index: null,
        analysis_inputs_hash: "STALE",
      },
      is_stale: true,
      created_by: 1,
      created_at: "2026-04-01T00:00:00Z",
    };

    const comparison: Comparison = {
      id: 7,
      title: "Stored",
      description: "from server",
      content_hash: "h-server",
      created_by: 1,
      created_at: "2026-04-01T00:00:00Z",
      updated_at: "2026-04-15T00:00:00Z",
      forked_from_id: null,
      forked_at_hash: null,
      forked_from_title: null,
      members: [savedMember],
    };

    useAppState.getState().loadDraftFromComparison(comparison, qc);
    const d = useAppState.getState().activeDraft!;
    expect(d.id).toBe(7);
    expect(d.baseHash).toBe("h-server");
    expect(d.title).toBe("Stored");
    expect(d.description).toBe("from server");
    expect(d.members).toHaveLength(1);
    const m = d.members[0]!;
    // Render parameters preserved from the server state
    expect(m.id).toBe(999);
    expect(m.exposure_id).toBe(exposureId);
    expect(m.band_height).toBe(1.5);
    expect(m.y_offset).toBe(0.2);
    expect(m.normalization).toBe("max");
    expect(m.color_override).toBe("#ff0000");
    expect(m.label_override).toBe("Lipid A");
    expect(m.q_window_min).toBe(0.05);
    expect(m.q_window_max).toBe(0.30);
    // Snapshot recovery: cache-derived, NOT the stale stored snapshot
    expect(m.snapshot).not.toBeUndefined();
    expect(m.snapshot!.analysis_inputs_hash).toBe("abcd");
    expect(m.snapshot!.effective_peaks).toHaveLength(2);
    expect(m.snapshot!.effective_peaks[0]!.id).toBe(1);
  });
});

// Type smoke
const _: ActiveDraft | null = null;
void _;
