/**
 * Compare baseHash capture-moment — Compare UX C-15 Step 4c.
 *
 * `draft.baseHash` is the comparison `content_hash` the draft was forked off.
 * It rides the eventual save as `expected_content_hash` so the backend can
 * detect a conflicting concurrent update.
 *
 * This file pins the capture-moment contract on the EXISTING behavior:
 *
 *   - baseHash is captured by the draft-creation path (`loadDraftFromComparison`
 *     → `fromComparison`, which copies `c.content_hash`). It is NOT recaptured
 *     by subsequent edit actions (`setDraftTitle` etc. only spread `...cur`).
 *   - A foreign update BEFORE the draft is created ⇒ the draft is built from
 *     the latest cached comparison (h2).
 *   - A foreign update AFTER the draft is created ⇒ baseHash stays pinned to
 *     the hash at creation time (h1); `loadDraftFromComparison` no-ops once a
 *     draft for that id already exists, so the drift cannot rebase baseHash.
 *
 * If a future refactor changes baseHash semantics, these tests are the
 * tripwire — do not silently relax them.
 */
import { describe, it, expect, beforeEach } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { useAppState } from "../src/state";
import type { Comparison } from "../src/api";

function makeQc(): QueryClient {
  return new QueryClient({
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
      mutations: { retry: false },
    },
  });
}

function makeComparison(over: Partial<Comparison> = {}): Comparison {
  return {
    id: 42,
    title: "Cubic vs Hex",
    description: null,
    content_hash: "h1",
    created_by: 7,
    created_at: "2026-01-01",
    updated_at: "2026-01-01",
    forked_from_id: null,
    forked_at_hash: null,
    forked_from_title: null,
    view_grouping_mode: null,
    view_show_peak_ticks: null,
    view_show_peak_labels: null,
    last_event_at: null,
    members: [],
    ...over,
  };
}

describe("Compare baseHash capture-moment — Compare UX C-15", () => {
  beforeEach(() => {
    sessionStorage.clear();
    localStorage.clear();
    useAppState.setState({ activeDraft: null });
  });

  it("captures baseHash at draft-creation time, not on later edits", () => {
    const qc = makeQc();
    const comparison = makeComparison({ content_hash: "h1" });
    // First edit — the draft-creation path captures content_hash → baseHash.
    useAppState.getState().loadDraftFromComparison(comparison, qc);
    // A later edit action only spreads `...cur` — it must not touch baseHash.
    useAppState.getState().setDraftTitle("renamed");
    const draft = useAppState.getState().activeDraft;
    expect(draft?.baseHash).toBe("h1");
    expect(draft?.title).toBe("renamed");
  });

  it("a foreign update BEFORE first edit ⇒ draft built from the latest (h2)", () => {
    const qc = makeQc();
    // The cache already reflects the foreign update by the time the user
    // makes the first edit — the draft forks off h2.
    const latest = makeComparison({ content_hash: "h2" });
    useAppState.getState().loadDraftFromComparison(latest, qc);
    expect(useAppState.getState().activeDraft?.baseHash).toBe("h2");
  });

  it("a foreign update AFTER first edit ⇒ baseHash stays pinned (h1)", () => {
    const qc = makeQc();
    // First edit forks the draft off h1.
    useAppState.getState().loadDraftFromComparison(makeComparison({ content_hash: "h1" }), qc);
    expect(useAppState.getState().activeDraft?.baseHash).toBe("h1");

    // A foreign update lands. The hydration effect re-runs
    // `loadDraftFromComparison` with the drifted comparison — but it no-ops
    // because a draft for id 42 is already active. baseHash stays pinned.
    useAppState.getState().loadDraftFromComparison(makeComparison({ content_hash: "h2" }), qc);
    // A subsequent edit also leaves baseHash untouched.
    useAppState.getState().setDraftTitle("edited after drift");
    const draft = useAppState.getState().activeDraft;
    expect(draft?.baseHash).toBe("h1");
    expect(draft?.title).toBe("edited after drift");
  });
});
