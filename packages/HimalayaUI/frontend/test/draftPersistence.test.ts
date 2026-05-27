/**
 * Compare-era draft state + sessionStorage persistence (Plan §Phase 4, Task 4.3).
 *
 * I5.3 (#184): the Compare-only create/fork/membership actions (startNewDraft,
 * setDraftTitle/Description, addMember, removeMember, discardDraft,
 * loadDraftFromComparison) were removed with the retired Compare page. This file
 * now covers only the KEPT slice — the member-edit actions consumed by the
 * shared series-builder render core (updateMember / reorderMembers / resizeBands)
 * and the draft.ts sessionStorage round-trip + schema-version guards. Drafts are
 * seeded directly via setState; the kept actions auto-mirror to sessionStorage.
 */
import { describe, it, expect, beforeEach } from "vitest";
import {
  loadDraftFromSession,
  persistDraftToSession,
  emptyDraft,
  COMPARE_DRAFT_KEY,
  COMPARE_DRAFT_VERSION,
  type ActiveDraft,
  type DraftMember,
} from "../src/lib/comparison/draft";
import { useAppState } from "../src/state";

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

/** Seed an active draft with the given members (display_order = index). */
function seedDraft(members: DraftMember[]): void {
  const draft: ActiveDraft = {
    ...emptyDraft(),
    members: members.map((m, i) => ({ ...m, display_order: i })),
  };
  useAppState.setState({ activeDraft: draft });
}

function resetStore(): void {
  useAppState.setState({ activeDraft: null });
}

describe("compare-era draft member edits + sessionStorage", () => {
  beforeEach(() => {
    sessionStorage.clear();
    resetStore();
  });

  it("updateMember partially mutates the indexed member", () => {
    seedDraft([makeMember({ exposure_id: 101 })]);
    useAppState.getState().updateMember(0, { y_offset: 0.5, normalization: "max" });
    const m = useAppState.getState().activeDraft!.members[0]!;
    expect(m.y_offset).toBe(0.5);
    expect(m.normalization).toBe("max");
    expect(m.exposure_id).toBe(101);
  });

  it("reorderMembers reorders by index list and renumbers display_order", () => {
    seedDraft([
      makeMember({ exposure_id: 101 }),
      makeMember({ exposure_id: 102 }),
      makeMember({ exposure_id: 103 }),
    ]);
    useAppState.getState().reorderMembers([2, 0, 1]);
    const ms = useAppState.getState().activeDraft!.members;
    expect(ms.map((m) => m.exposure_id)).toEqual([103, 101, 102]);
    expect(ms.map((m) => m.display_order)).toEqual([0, 1, 2]);
  });

  it("resizeBands shifts a band's height and clamps the neighbor", () => {
    seedDraft([
      makeMember({ exposure_id: 101 }),
      makeMember({ exposure_id: 102 }),
    ]);
    // memberIdx=0 grows by 20px out of 200px total → +0.10 ratio; neighbor shrinks.
    useAppState.getState().resizeBands(0, 20, 200);
    const ms = useAppState.getState().activeDraft!.members;
    expect(ms[0]!.band_height).toBeCloseTo(1.10, 5);
    expect(ms[1]!.band_height).toBeCloseTo(0.90, 5);
  });

  it("a kept member edit mirrors the draft to sessionStorage", () => {
    seedDraft([makeMember({ exposure_id: 101 })]);
    useAppState.getState().updateMember(0, { y_offset: 0.25 });
    const raw = sessionStorage.getItem(COMPARE_DRAFT_KEY);
    expect(raw).not.toBeNull();
    const parsed = JSON.parse(raw!);
    expect(parsed.version).toBe(COMPARE_DRAFT_VERSION);
    expect(parsed.draft.members[0].y_offset).toBe(0.25);
  });

  it("round-trips a persisted draft through loadDraftFromSession", () => {
    const draft: ActiveDraft = { ...emptyDraft(), title: "title" };
    persistDraftToSession(draft);
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
    expect(loadDraftFromSession()).toBeNull();
  });

  it("loadDraftFromSession ignores malformed JSON", () => {
    sessionStorage.setItem(COMPARE_DRAFT_KEY, "{not json");
    expect(loadDraftFromSession()).toBeNull();
  });
});

// Type smoke
const _: ActiveDraft | null = null;
void _;
