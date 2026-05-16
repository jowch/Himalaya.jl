/**
 * Compare-page e2e smoke (Plan §Phase 13, Task 13.3).
 *
 * Single end-to-end happy-path covering:
 *   1. Navigate to /experiments/1/compare
 *   2. Click "+ New" → enter edit mode
 *   3. Add 2 members via the picker
 *   4. Save → verify navigation to review URL
 *   5. Reload → verify saved state survives
 *   6. Fork (testing under a *different* current user) → land in edit mode
 *   7. Save the fork → verify navigation to fork's review URL
 *
 * All `/api/*` endpoints are mocked via `page.route`. The conflict-modal
 * flow is intentionally NOT covered here — the unit-tested ConflictModal
 * exercises that contract more thoroughly than a Playwright run can.
 *
 * Selectors are stable `data-testid`s per the Phase 0 selector table.
 */
import { test, expect, type Page, type Route } from "@playwright/test";

// ─── Fixtures ──────────────────────────────────────────────────────────────

const EXPERIMENT = {
  id: 1, name: "SSRL Test", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-05-01",
};
const SAMPLES = [
  { id: 10, experiment_id: 1, display_name: "D1", name: "cubic_run03", notes: null, tags: [] },
  { id: 11, experiment_id: 1, display_name: "D2", name: "hex_run01",   notes: null, tags: [] },
];
const EXPOSURES = [
  {
    id: 100, sample_id: 10, filename: "cubic.dat", kind: "file",
    selected: true, status: "accepted", image_path: null, image_version: "",
    tags: [], sources: [],
    trace_hash: "tr-100", analysis_inputs_hash: "h-100",
  },
  {
    id: 101, sample_id: 11, filename: "hex.dat", kind: "file",
    selected: true, status: "accepted", image_path: null, image_version: "",
    tags: [], sources: [],
    trace_hash: "tr-101", analysis_inputs_hash: "h-101",
  },
];
const TRACE = { q: [0.1, 0.2, 0.3], I: [10, 5, 2], sigma: [1, 1, 1] };

const PICKER_SAMPLES = SAMPLES.map((s) => {
  const all_exposures = EXPOSURES.filter((e) => e.sample_id === s.id).map((e) => ({
    id: e.id,
    sample_id: e.sample_id,
    filename: e.filename,
    selected: e.selected,
  }));
  const indexing = all_exposures.find((e) => e.selected) ?? all_exposures[all_exposures.length - 1];
  return {
    sample: s,
    indexing_exposure_id: indexing ? indexing.id : null,
    all_exposures,
  };
});

// Initial users — only "alice" exists. Fork test creates "bob" via setting
// localStorage username before reload (see test 2 below).
const USERS = [{ id: 1, username: "alice", first_name: null, last_name: null }];

const EMPTY_COMPARISONS_LIST: ReadonlyArray<unknown> = [];

// Listing-projection fields (#136/#137). Every comparison object that lands
// in `comparisonsByExp` is also served as a sidebar listing row, and the
// redesigned ComparisonSidebar (Compare UX Phase F) consumes these fields,
// so each mock comparison must carry the projection shape.
const PROJECTION_DEFAULTS = {
  view_grouping_mode: null,
  view_show_peak_ticks: null,
  view_show_peak_labels: null,
  last_event_at: "2026-05-06T00:00:00Z",
  author_username: "alice",
  member_count: 2,
  member_phases: [] as string[],
  member_phase_count: 0,
  has_stale_members: false,
};

interface ServerState {
  users: typeof USERS;
  comparisonsByExp: Map<number, Array<{ id: number; title: string; description: string | null;
    content_hash: string; created_by: number | null; created_at: string;
    updated_at: string; forked_from_id: number | null; forked_at_hash: string | null;
  }>>;
  comparisonsById: Map<number, unknown>;
  pinsByUser: Map<string, number[]>;
  nextComparisonId: number;
}

function makeState(): ServerState {
  return {
    users: [...USERS],
    comparisonsByExp: new Map([[1, []]]),
    comparisonsById: new Map(),
    pinsByUser: new Map(),
    nextComparisonId: 1,
  };
}

// ─── Mock helpers ──────────────────────────────────────────────────────────

async function jsonOK(route: Route, body: unknown, status = 200): Promise<void> {
  await route.fulfill({
    status,
    contentType: "application/json",
    body: JSON.stringify(body),
  });
}

async function mockApi(page: Page, state: ServerState): Promise<void> {
  // Users.
  await page.route("**/api/users", (r) => jsonOK(r, state.users));

  // Experiments + samples + exposures (read-only for this smoke).
  await page.route("**/api/experiments", (r) => jsonOK(r, [EXPERIMENT]));
  await page.route("**/api/experiments/1", (r) => jsonOK(r, EXPERIMENT));
  await page.route("**/api/experiments/1/samples", (r) => jsonOK(r, SAMPLES));
  await page.route("**/api/experiments/1/exposures*", (r) =>
    jsonOK(r, EXPOSURES));
  await page.route("**/api/experiments/1/sample-tags", (r) => jsonOK(r, []));
  await page.route("**/api/experiments/1/picker-samples", (r) =>
    jsonOK(r, PICKER_SAMPLES));
  await page.route("**/api/samples/10/exposures*", (r) =>
    jsonOK(r, EXPOSURES.filter((e) => e.sample_id === 10)));
  await page.route("**/api/samples/11/exposures*", (r) =>
    jsonOK(r, EXPOSURES.filter((e) => e.sample_id === 11)));
  await page.route("**/api/samples/10/messages", (r) => jsonOK(r, []));
  await page.route("**/api/samples/11/messages", (r) => jsonOK(r, []));

  // Per-exposure trace + peaks + indices + groups + entity for the cache
  // warm-up needed by `useMemberTraces` / `computeMemberSnapshot`.
  for (const e of EXPOSURES) {
    await page.route(`**/api/exposures/${e.id}`, (r) => jsonOK(r, e));
    await page.route(`**/api/exposures/${e.id}/trace`, (r) => jsonOK(r, TRACE));
    await page.route(`**/api/exposures/${e.id}/peaks`, (r) => jsonOK(r, []));
    await page.route(`**/api/exposures/${e.id}/indices`, (r) => jsonOK(r, []));
    await page.route(`**/api/exposures/${e.id}/groups`, (r) => jsonOK(r, []));
  }

  // Picker support routes.
  await page.route("**/api/users/*/recently-picked-exposures*", (r) => jsonOK(r, []));

  // Comparison list endpoints.
  await page.route("**/api/experiments/1/comparisons", (r) =>
    jsonOK(r, state.comparisonsByExp.get(1) ?? []));
  await page.route("**/api/comparisons", async (r) => {
    const req = r.request();
    if (req.method() === "POST") {
      const body = req.postDataJSON() as {
        title: string; description?: string | null;
        members: Array<{ exposure_id: number | null; display_order: number }>;
        forked_from_id?: number | null; forked_at_hash?: string | null;
      };
      const id = state.nextComparisonId++;
      const comp = {
        id,
        ...PROJECTION_DEFAULTS,
        member_count: body.members.length,
        title: body.title,
        description: body.description ?? null,
        content_hash: `h-${id}`,
        created_by: 1, // alice's id
        created_at: "2026-05-06T00:00:00Z",
        updated_at: "2026-05-06T00:00:00Z",
        forked_from_id: body.forked_from_id ?? null,
        forked_at_hash: body.forked_at_hash ?? null,
        forked_from_title: null,
        members: body.members.map((m, i) => ({
          id: i + 1,
          comparison_id: id,
          exposure_id: m.exposure_id,
          display_order: m.display_order,
          band_height: 1.0,
          y_offset: 0.0,
          normalization: "none",
          color_override: null,
          label_override: null,
          q_window_min: null,
          q_window_max: null,
          peak_display: null,
          snapshot: { effective_peaks: [], confirmed_index: null,
                      analysis_inputs_hash: `h-${m.exposure_id}` },
          is_stale: false,
          created_by: 1,
          created_at: "2026-05-06T00:00:00Z",
        })),
      };
      state.comparisonsById.set(id, comp);
      const list = state.comparisonsByExp.get(1) ?? [];
      list.push(comp);
      state.comparisonsByExp.set(1, list);
      await jsonOK(r, comp);
      return;
    }
    // GET: global list
    await jsonOK(r, Array.from(state.comparisonsById.values()));
  });

  await page.route("**/api/comparisons/*", async (r) => {
    const url = new URL(r.request().url());
    const m = url.pathname.match(/\/api\/comparisons\/(\d+)$/);
    if (m) {
      const id = Number(m[1]);
      const comp = state.comparisonsById.get(id);
      if (!comp) return jsonOK(r, { error: "not found" }, 404);
      return jsonOK(r, comp);
    }
    return r.continue();
  });

  await page.route("**/api/comparisons/*/forks", (r) => jsonOK(r, []));
  await page.route("**/api/comparisons/*/messages", async (r) => {
    if (r.request().method() === "POST") {
      const body = r.request().postDataJSON() as { body: string };
      return jsonOK(r, {
        id: 1, comparison_id: 1, author_id: 1, author: "alice",
        body: body.body, created_at: "2026-05-06T00:00:00Z",
      }, 201);
    }
    return jsonOK(r, []);
  });

  await page.route("**/api/comparisons/*/submit", async (r) => {
    if (r.request().method() === "POST") {
      const url = new URL(r.request().url());
      const m = url.pathname.match(/\/api\/comparisons\/(\d+)\/submit$/);
      const id = Number(m![1]);
      const body = r.request().postDataJSON() as {
        title: string; description?: string | null;
        members: Array<{ exposure_id: number | null; display_order: number }>;
        expected_content_hash?: string;
      };
      const prior = state.comparisonsById.get(id) as { content_hash?: string } | undefined;
      // Surface a 409 if the client's expected hash doesn't match the
      // server's current. This smoke does NOT exercise the conflict modal,
      // but the path is plumbed for completeness.
      if (body.expected_content_hash !== undefined
          && prior?.content_hash !== body.expected_content_hash) {
        return jsonOK(r, { error: "conflict",
          current_hash: prior?.content_hash, current_state: prior }, 409);
      }
      const next = { ...(prior ?? {}), id, title: body.title,
        description: body.description ?? null,
        content_hash: `h-${id}-v2`, updated_at: "2026-05-06T01:00:00Z",
        created_by: 1, created_at: "2026-05-06T00:00:00Z",
        forked_from_id: null, forked_at_hash: null, forked_from_title: null,
        members: body.members.map((m2, i) => ({
          id: i + 1, comparison_id: id,
          exposure_id: m2.exposure_id, display_order: m2.display_order,
          band_height: 1.0, y_offset: 0.0, normalization: "none",
          color_override: null, label_override: null,
          q_window_min: null, q_window_max: null, peak_display: null,
          snapshot: { effective_peaks: [], confirmed_index: null,
                      analysis_inputs_hash: `h-${m2.exposure_id}` },
          is_stale: false, created_by: 1,
          created_at: "2026-05-06T00:00:00Z",
        })),
      };
      state.comparisonsById.set(id, next);
      return jsonOK(r, next);
    }
    return r.continue();
  });

  await page.route("**/api/comparisons/*/pin", async (r) => {
    const m = r.request().url().match(/\/api\/comparisons\/(\d+)\/pin$/);
    const id = Number(m![1]);
    const username = r.request().headers()["x-username"] ?? "alice";
    const list = state.pinsByUser.get(username) ?? [];
    if (r.request().method() === "POST") {
      state.pinsByUser.set(username, [id, ...list.filter((x) => x !== id)]);
      return jsonOK(r, { comparison_id: id, pinned: true });
    }
    state.pinsByUser.set(username, list.filter((x) => x !== id));
    return jsonOK(r, { comparison_id: id, pinned: false });
  });
  await page.route("**/api/users/me/comparison-pins", (r) => {
    const username = r.request().headers()["x-username"] ?? "alice";
    return jsonOK(r, state.pinsByUser.get(username) ?? []);
  });

  // SSE — drain immediately so the EventSource doesn't hang the page.
  await page.route("**/api/events", (r) =>
    r.fulfill({
      status: 200,
      contentType: "text/event-stream",
      body: "",
    }));
}

async function seedState(page: Page, extra: Record<string, unknown> = {}): Promise<void> {
  await page.addInitScript((state) => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({
        state: {
          username: "alice",
          activePage: "compare",
          tutorialSeen: true,
          theme: "dark",
          activeExperimentId: 1,
          ...state,
        },
        version: 3,
      }),
    );
  }, extra);
}

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() => { localStorage.clear(); sessionStorage.clear(); });
});

// ─── Smoke: full create flow ────────────────────────────────────────────────

test("compare smoke: create → submit → review", async ({ page }) => {
  const state = makeState();
  await mockApi(page, state);
  await seedState(page);
  await page.goto("/experiments/1/compare");

  // 1. Sidebar visible, empty state.
  await expect(page.getByTestId("comparison-sidebar")).toBeVisible();
  await expect(page.getByTestId("comparison-sidebar-empty")).toBeVisible();

  // 2. Click "+ New" → land in edit mode.
  await page.getByTestId("comparison-new").click();
  await expect(page).toHaveURL(/\/experiments\/1\/compare\/new$/);
  await expect(page.getByTestId("compare-page-edit")).toBeVisible();

  // 3. Fill the title (Compare UX C-13 — title is now an InlineEditableText
  // inside CompareTitleStrip; click the rest span to enter edit mode).
  await page.getByTestId("compare-title").click();
  await page.getByTestId("compare-title").fill("My first comparison");

  // 4. Inline picker panel is in the right slot (PR2). Tick each sample
  // row to immediate-commit add — no modal, no batch "Add N selected" footer.
  await expect(page.getByTestId("comparison-picker-panel")).toBeVisible();
  // Both samples 10 and 11 surface as picker rows in experiment-1's panel.
  const rows = page.getByTestId("picker-row");
  await expect(rows).toHaveCount(2);
  for (let i = 0; i < 2; i++) {
    await rows.nth(i).locator('input[type="checkbox"]').check();
  }
  // Both members are now in the draft (immediate-commit).

  // 5. Save — verify the request lands and we navigate to review.
  await page.getByTestId("save-pill").click();
  // Server state must contain the new comparison.
  await expect(page).toHaveURL(/\/experiments\/1\/compare\/\d+$/);
  await expect(page.getByTestId("compare-review-plot")).toBeVisible();
});

// ─── Smoke: reload preserves saved state ────────────────────────────────────

test("compare smoke: review survives full page reload", async ({ page }) => {
  const state = makeState();
  // Pre-seed a saved comparison.
  const comp = {
    id: 1,
    ...PROJECTION_DEFAULTS,
    member_count: 1,
    title: "Saved comparison",
    description: null,
    content_hash: "h-1",
    created_by: 1,
    created_at: "2026-05-06T00:00:00Z",
    updated_at: "2026-05-06T00:00:00Z",
    forked_from_id: null,
    forked_at_hash: null,
    forked_from_title: null,
    members: [{
      id: 1, comparison_id: 1, exposure_id: 100, display_order: 0,
      band_height: 1.0, y_offset: 0.0, normalization: "none",
      color_override: null, label_override: null,
      q_window_min: null, q_window_max: null, peak_display: null,
      snapshot: { effective_peaks: [], confirmed_index: null,
                  analysis_inputs_hash: "h-100" },
      is_stale: false, created_by: 1, created_at: "2026-05-06T00:00:00Z",
    }],
  };
  state.comparisonsById.set(1, comp);
  state.comparisonsByExp.set(1, [comp]);
  state.nextComparisonId = 2;

  await mockApi(page, state);
  await seedState(page);
  await page.goto("/experiments/1/compare/1");
  await expect(page.getByTestId("compare-review-plot")).toBeVisible();

  // Reload — review page should re-render against the same mocked server.
  await page.reload();
  await expect(page.getByTestId("compare-review-plot")).toBeVisible();
  // Edit affordance is present (alice authored the comparison).
  await expect(page.getByTestId("comparison-edit")).toBeVisible();
});

// ─── Smoke: fork as a different user ────────────────────────────────────────

test("compare smoke: non-author sees Fork → can save fork", async ({ page }) => {
  const state = makeState();
  // alice authored the comparison; we'll log in as "bob" who can only fork.
  const comp = {
    id: 1,
    ...PROJECTION_DEFAULTS,
    member_count: 1,
    title: "Alice's comparison",
    description: null,
    content_hash: "h-1",
    created_by: 1, // alice
    created_at: "2026-05-06T00:00:00Z",
    updated_at: "2026-05-06T00:00:00Z",
    forked_from_id: null,
    forked_at_hash: null,
    forked_from_title: null,
    members: [{
      id: 1, comparison_id: 1, exposure_id: 100, display_order: 0,
      band_height: 1.0, y_offset: 0.0, normalization: "none",
      color_override: null, label_override: null,
      q_window_min: null, q_window_max: null, peak_display: null,
      snapshot: { effective_peaks: [], confirmed_index: null,
                  analysis_inputs_hash: "h-100" },
      is_stale: false, created_by: 1, created_at: "2026-05-06T00:00:00Z",
    }],
  };
  state.comparisonsById.set(1, comp);
  state.comparisonsByExp.set(1, [comp]);
  state.nextComparisonId = 2;
  state.users = [
    ...USERS,
    { id: 2, username: "bob", first_name: null, last_name: null },
  ];

  await mockApi(page, state);
  await seedState(page, { username: "bob" });
  await page.goto("/experiments/1/compare/1");

  // Non-author → Fork (mutually exclusive with Edit).
  await expect(page.getByTestId("comparison-fork")).toBeVisible();
  await expect(page.getByTestId("comparison-edit")).toHaveCount(0);

  await page.getByTestId("comparison-fork").click();
  // Fork lands the user in /new (the create flow with lineage pre-populated).
  await expect(page).toHaveURL(/\/experiments\/1\/compare\/new$/);
  await expect(page.getByTestId("compare-page-edit")).toBeVisible();

  // Edit the title and save the fork (Compare UX C-13 — InlineEditableText
  // title + SavePill replace the legacy input + Save button).
  await page.getByTestId("compare-title").click();
  await page.getByTestId("compare-title").fill("Bob's fork");
  await page.getByTestId("save-pill").click();
  // Lands on the new fork's review page (id=2 since nextComparisonId starts at 2).
  await expect(page).toHaveURL(/\/experiments\/1\/compare\/\d+$/);
  await expect(page.getByTestId("compare-review-plot")).toBeVisible();
});
