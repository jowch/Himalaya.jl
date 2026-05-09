/**
 * Bundle B paper-cuts smoke (Task 6).
 *
 * Headless Playwright coverage for three user-visible Bundle B fixes:
 *   - #77 cold-load: persisted activePage='compare' at URL '/' must mount
 *     the Compare page (not just the rocker).
 *   - #78 edit-mode hint: the right-slot hint at /compare/new must show
 *     for an empty draft and disappear once members are added.
 *   - #75 ForksPopover dismissal: the popover must close on Esc and on
 *     outside-click; aria-expanded reflects state.
 *
 * `mockApi` / `seedState` / `makeState` are inlined (verbatim from
 * `compare.spec.ts`) per the existing project pattern — no shared helper
 * module yet.
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

const USERS = [{ id: 1, username: "alice", first_name: null, last_name: null }];

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
  await page.route("**/api/samples/10/exposures*", (r) =>
    jsonOK(r, EXPOSURES.filter((e) => e.sample_id === 10)));
  await page.route("**/api/samples/11/exposures*", (r) =>
    jsonOK(r, EXPOSURES.filter((e) => e.sample_id === 11)));
  await page.route("**/api/samples/10/messages", (r) => jsonOK(r, []));
  await page.route("**/api/samples/11/messages", (r) => jsonOK(r, []));

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
        title: body.title,
        description: body.description ?? null,
        content_hash: `h-${id}`,
        created_by: 1,
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

  // SSE — drain immediately.
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

// ─── #77: cold-load with persisted activePage='compare' ────────────────────

test("#77 cold-load with persisted activePage='compare' + URL '/' mounts Compare", async ({ page }) => {
  const state = makeState();
  await mockApi(page, state);
  // seedState defaults: activePage='compare', activeExperimentId=1.
  await seedState(page);
  await page.goto("/");

  // Zustand → URL sync effect must redirect to /experiments/1/compare.
  await expect(page).toHaveURL(/\/experiments\/1\/compare(\/|$)/);
  await expect(page.getByTestId("comparison-sidebar")).toBeVisible();
});

test("#77 cold-load with no active experiment redirects to /compare/all", async ({ page }) => {
  const state = makeState();
  await mockApi(page, state);
  // Inline the localStorage write rather than going through `seedState`:
  // Playwright serializes `addInitScript` args via structured clone, which
  // strips `undefined` values — so passing `{ activeExperimentId: undefined }`
  // to seedState's `extra` would NOT override the default `activeExperimentId: 1`.
  // Build the persisted state explicitly to omit the key.
  await page.addInitScript(() => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({
        state: {
          username: "alice",
          activePage: "compare",
          tutorialSeen: true,
          theme: "dark",
          // activeExperimentId intentionally absent → Zustand default (undefined).
        },
        version: 3,
      }),
    );
  });
  await page.goto("/");

  await expect(page).toHaveURL(/\/compare\/all$/);
  await expect(page.getByTestId("comparison-sidebar")).toBeVisible();
});

// ─── #78: edit-mode right-slot hint conditional ────────────────────────────

test("#78 edit-mode right-slot hint hides when draft has members", async ({ page }) => {
  const state = makeState();
  await mockApi(page, state);
  await seedState(page);
  await page.goto("/experiments/1/compare/new");

  await expect(page.getByTestId("compare-page-edit")).toBeVisible();
  // Empty draft → hint visible.
  await expect(page.getByTestId("compare-edit-right-hint")).toBeVisible();

  // Add a member via the picker.
  await page.getByTestId("compare-edit-add-traces").click();
  await expect(page.getByTestId("comparison-picker")).toBeVisible();
  const rows = page.getByTestId("picker-row");
  await expect(rows.first()).toBeVisible();
  await rows.first().locator('input[type="checkbox"]').check();
  await page.getByTestId("comparison-picker-add").click();
  await expect(page.getByTestId("comparison-picker")).not.toBeVisible();

  // Hint must now be gone.
  await expect(page.getByTestId("compare-edit-right-hint")).toHaveCount(0);
});

// ─── #75: ForksPopover dismissal affordances ───────────────────────────────

test("#75 ForksPopover Esc + outside-click close + aria-expanded", async ({ page }) => {
  const state = makeState();

  // Pre-seed parent + fork pair on the mocked server.
  const parent = {
    id: 1,
    title: "Parent comparison",
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
  const fork = {
    id: 2,
    title: "My fork",
    description: null,
    content_hash: "h-2",
    created_by: 1,
    created_at: "2026-05-06T00:00:00Z",
    updated_at: "2026-05-06T00:00:00Z",
    forked_from_id: 1,
    forked_at_hash: "h-1",
    forked_from_title: "Parent comparison",
    members: [{
      id: 2, comparison_id: 2, exposure_id: 101, display_order: 0,
      band_height: 1.0, y_offset: 0.0, normalization: "none",
      color_override: null, label_override: null,
      q_window_min: null, q_window_max: null, peak_display: null,
      snapshot: { effective_peaks: [], confirmed_index: null,
                  analysis_inputs_hash: "h-101" },
      is_stale: false, created_by: 1, created_at: "2026-05-06T00:00:00Z",
    }],
  };
  state.comparisonsById.set(1, parent);
  state.comparisonsById.set(2, fork);
  state.comparisonsByExp.set(1, [parent, fork]);
  state.nextComparisonId = 3;

  await mockApi(page, state);
  // Override the parent's /forks endpoint to return [fork]. Playwright's
  // page.route handlers are last-registered-first, so this wins over the
  // generic `**/api/comparisons/*/forks` registered inside mockApi.
  await page.route("**/api/comparisons/1/forks", (r) =>
    r.fulfill({
      status: 200,
      contentType: "application/json",
      body: JSON.stringify([fork]),
    }));

  await seedState(page);
  await page.goto("/experiments/1/compare/1");
  await expect(page.getByTestId("compare-review-plot")).toBeVisible();

  const trigger = page.getByTestId("comparison-forks-trigger");
  await expect(trigger).toBeVisible();
  await expect(trigger).toHaveAttribute("aria-expanded", "false");

  // Open via click.
  await trigger.click();
  await expect(page.getByTestId("comparison-forks-popover")).toBeVisible();
  await expect(trigger).toHaveAttribute("aria-expanded", "true");
  // Sanity: the fork row landed (so the override is wired correctly and
  // the popover content reflects state, not just visibility).
  await expect(page.getByTestId("comparison-forks-row")).toHaveCount(1);

  // Esc closes.
  await page.keyboard.press("Escape");
  await expect(page.getByTestId("comparison-forks-popover")).not.toBeVisible();
  await expect(trigger).toHaveAttribute("aria-expanded", "false");

  // Reopen + outside-click closes. Click in the sidebar's empty area at
  // a corner away from any interactive element.
  await trigger.click();
  await expect(page.getByTestId("comparison-forks-popover")).toBeVisible();
  await expect(trigger).toHaveAttribute("aria-expanded", "true");

  await page.getByTestId("comparison-sidebar").click({ position: { x: 5, y: 5 } });
  await expect(page.getByTestId("comparison-forks-popover")).not.toBeVisible();
  await expect(trigger).toHaveAttribute("aria-expanded", "false");
});
