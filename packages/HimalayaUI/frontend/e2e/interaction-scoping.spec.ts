/**
 * E2E: SeriesScopingPage interaction migration — mocked Playwright spec (Task 5.3).
 *
 * Exercises the shell-owned InteractionDock (member cursor stepper + up-link),
 * the roving-tabindex cursor over the worksheet rows (useListCursor in headless /
 * scope-focus mode), Alt+↑/↓ keyboard reorder, drag-reorder, undo/redo, and the
 * SR live-region announcement after SeriesScopingPage was migrated to the
 * shell-level interaction registry.
 *
 * Dock testid contract (InteractionDock / DockStepper / DockUpLink):
 *   dock-prev-member / dock-next-member / dock-member-count  — cursor stepper (in-page)
 *   dock-up-link                                             — DockUpLink "‹ All series"
 *   NO dock-primary (scoping has no Enter-target primary action)
 *
 * Scope / row contract:
 *   [data-interaction-scope]  — scope div (tabIndex=-1); auto-focuses after isLoading→false
 *   [data-reorder-row]        — each member row wrapper; spreads cursor.rowProps (tabIndex 0/-1)
 *   [data-testid="scope-sample-row"]  — the ScopeSampleRow child inside each row wrapper
 *   [data-testid="reorder-announcement"]  — polite SR live region for move feedback
 *
 * Covered cases:
 *   IS-1  Exactly one [data-reorder-row][tabindex="0"] at a time; ArrowDown and
 *          dock-next-member both move it.
 *   IS-2  Cold-load Alt+ArrowDown (no prior click): the cursored member's order
 *          changes AND the cursor follows the item to its new slot. Verifies that
 *          the keyboard layer fires without scope interaction AND that cursor is
 *          seeded at page load (useListCursor initializes to ids[0]).
 *   IS-3  Alt+ArrowUp on the second row announces the correct move message in the
 *          SR live region ([data-testid="reorder-announcement"]).
 *   IS-4  HTML5 drag-and-drop reorder changes the display order, proving the
 *          roving-tabindex wiring did not break drag handlers.
 *   IS-5  ⌘Z undoes a keyboard reorder; ⌘⇧Z redoes it.
 */

import { test, expect, type Page } from "@playwright/test";

// ---------------------------------------------------------------------------
// Fixtures — mirrored from series-scoping.spec.ts
// ---------------------------------------------------------------------------

const EXPERIMENT = {
  id: 1, name: "SSRL Test", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-05-01",
};

/** A picker row whose sample has a known value for the ordering key ("ratio"). */
function memberRow(id: number, name: string, ratio: string) {
  return {
    sample: {
      id, experiment_id: 1, name, display_name: null, notes: null,
      tags: [{ id, key: "ratio", value: ratio, source: "manual" }],
    },
    indexing_exposure_id: null,
    all_exposures: [],
  };
}

/** A picker row whose sample lacks the ordering key — rendered as a loose candidate. */
function candidateRow(id: number, name: string) {
  return {
    sample: {
      id, experiment_id: 1, name, display_name: null, notes: null,
      tags: [],
    },
    indexing_exposure_id: null,
    all_exposures: [],
  };
}

/** Seed localStorage so onboarding/tutorial prompts do not gate navigation. */
async function seedState(page: Page): Promise<void> {
  await page.addInitScript(() => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({
        state: { username: "alice", tutorialSeen: true, theme: "dark" },
        version: 5,
      }),
    );
  });
}

/**
 * Register mocked API routes for the SeriesScopingPage.
 * Route set mirrors series-scoping.spec.ts (no batch capture needed here).
 * Two member rows (A_1to1, B_2to1) + one loose candidate (C_unlabelled).
 */
async function mockScoping(page: Page): Promise<void> {
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([EXPERIMENT]) }));

  // "ratio" appears twice → proposed as the ordering variable.
  await page.route("**/api/sample-tags", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([
        { key: "ratio", value: "1:1" },
        { key: "ratio", value: "2:1" },
      ]) }));

  // Two member samples (with "ratio") + one candidate (without).
  await page.route("**/api/picker-samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([
        memberRow(10, "A_1to1", "1:1"),
        memberRow(11, "B_2to1", "2:1"),
        candidateRow(12, "C_unlabelled"),
      ]) }));

  // Batch write endpoint (not called in interaction tests, but stub to avoid
  // unmatched-route errors if confirm is accidentally triggered).
  await page.route("**/api/samples/tags/batch", (r) =>
    r.fulfill({ status: 201, contentType: "application/json", body: "[]" }));

  // Folio listing — empty is fine (not exercised here).
  await page.route("**/api/series", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));

  // SSE: drain immediately so the App-level EventSource doesn't hang the page.
  await page.route("**/api/events*", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
}

// ---------------------------------------------------------------------------
// Setup
// ---------------------------------------------------------------------------

test.beforeEach(async ({ page }) => {
  await seedState(page);
});

// ---------------------------------------------------------------------------
// IS-1: Roving tabindex invariant — exactly one [data-reorder-row][tabindex="0"]
// ---------------------------------------------------------------------------

test("IS-1: exactly one [data-reorder-row][tabindex='0'] at a time; ArrowDown and dock-next-member both advance it", async ({ page }) => {
  await mockScoping(page);
  await page.goto("/series/new");

  // Wait for data to load and the cursor to initialize (ids=[10,11], cursor=10).
  await expect(page.getByTestId("dock-member-count")).toHaveText(/1 \/ 2/);

  // Initial roving state: exactly one row at tabindex=0 (the first, A_1to1).
  await expect(page.locator('[data-reorder-row][tabindex="0"]')).toHaveCount(1);

  // dock-next-member (cursor stepper) advances the cursor to the second row.
  await page.getByTestId("dock-next-member").click();
  await expect(page.getByTestId("dock-member-count")).toHaveText(/2 \/ 2/);
  // Roving invariant holds after stepper click.
  await expect(page.locator('[data-reorder-row][tabindex="0"]')).toHaveCount(1);

  // ArrowDown in the scope advances the cursor back (clock-wrap → or to next).
  // First click moved us to index 1 (last); pressing stepper again wraps or stays.
  // Reset cursor to first: move back via dock-prev-member.
  await page.getByTestId("dock-prev-member").click();
  await expect(page.getByTestId("dock-member-count")).toHaveText(/1 \/ 2/);
  // The scope div auto-focused after load. Press ArrowDown via the scope's onKeyDown.
  await page.keyboard.press("ArrowDown");
  await expect(page.getByTestId("dock-member-count")).toHaveText(/2 \/ 2/);
  // Still exactly one tabindex=0 row.
  await expect(page.locator('[data-reorder-row][tabindex="0"]')).toHaveCount(1);
});

// ---------------------------------------------------------------------------
// IS-2: Cold-load Alt+ArrowDown — reorder fires WITHOUT any prior click
// ---------------------------------------------------------------------------

test("IS-2: cold-load Alt+ArrowDown reorders the cursored member and cursor follows the item", async ({ page }) => {
  // INTENT: Alt+ArrowDown is a chorded key (not bare), so it is NOT scope-gated
  // by inPageScope(). It fires via the keyboard layer regardless of focus position.
  // The cold-load regression here is that the cursor must be seeded at page load
  // (useListCursor initializes to ids[0] = 10 = A_1to1) so the first keypress
  // actually has a cursored member to reorder. No prior click — not even a focus
  // click — must precede this assertion.
  await mockScoping(page);
  await page.goto("/series/new");

  // Wait for data to load: proposal resolves, order=[10,11], cursor=10 (A_1to1).
  // dock-member-count "1 / 2" confirms the cursor is seeded and rows are rendered.
  await expect(page.getByTestId("dock-member-count")).toHaveText(/1 \/ 2/);

  // Initial visual order: A_1to1 first (position 0), B_2to1 second (position 1).
  // Do NOT click anything before this key press — this is the cold-load gate.
  await expect(page.getByTestId("scope-sample-row").first()).toContainText("A_1to1");

  // Press Alt+ArrowDown — fires the `reorderDown` page action.
  // cursor.cursorId = 10 (A_1to1), at order index 0.
  // reorderDown: i=0, i < order.length-1 → moveRow(0, +1) → applyReorder(0,1)
  //   order changes: [10,11] → [11,10] (B_2to1 moves to display position 0)
  //   cursor.setCursor(10) → cursor follows the item to its new slot (index 1).
  await page.keyboard.press("Alt+ArrowDown");

  // Order changed: B_2to1 is now the first row.
  await expect(page.getByTestId("scope-sample-row").first()).toContainText("B_2to1");

  // Cursor followed the item: A_1to1 (id=10) is now at display index 1 but still
  // holds tabindex=0 (cursor-follows-item, the id-based invariant).
  await expect(page.locator('[data-reorder-row][tabindex="0"]')).toHaveCount(1);
  await expect(
    page.locator('[data-reorder-row][tabindex="0"]').first().getByTestId("scope-sample-row"),
  ).toContainText("A_1to1");
});

// ---------------------------------------------------------------------------
// IS-3: Reorder announces in the SR live region
// ---------------------------------------------------------------------------

test("IS-3: Alt+ArrowUp announces the move in [data-testid='reorder-announcement']", async ({ page }) => {
  await mockScoping(page);
  await page.goto("/series/new");
  await expect(page.getByTestId("dock-member-count")).toHaveText(/1 \/ 2/);

  // Move the cursor to B_2to1 (the second member, index 1).
  await page.getByTestId("dock-next-member").click();
  await expect(page.getByTestId("dock-member-count")).toHaveText(/2 \/ 2/);

  // Alt+ArrowUp with cursor on B_2to1 (i=1, delta=-1):
  //   moveRow(1, -1) → applyReorder(1, 0) → order [10,11] → [11,10]
  //   announceReorder("Moved B_2to1 to position 1 of 2.")
  await page.keyboard.press("Alt+ArrowUp");

  // The SR live region (aria-live="polite", sr-only) must contain the position message.
  await expect(page.getByTestId("reorder-announcement")).toContainText(
    "Moved B_2to1 to position 1 of 2.",
  );

  // The order also actually changed (B_2to1 is now first).
  await expect(page.getByTestId("scope-sample-row").first()).toContainText("B_2to1");
});

// ---------------------------------------------------------------------------
// IS-4: Drag-and-drop reorder changes the display order
// ---------------------------------------------------------------------------

test("IS-4: HTML5 drag from last row to first row changes the display order", async ({ page }) => {
  // INTENT: verify the roving-tabindex wiring did NOT break the HTML5 DnD
  // event handlers (useDragReorder: onDragStart / onDragOver / onDrop).
  // Initial order: A_1to1 (index 0), B_2to1 (index 1).
  // Drag from index 0 to index 1 → onReorder(0, 1) → applyReorder(0, 1)
  //   order: [10,11] → [11,10] → B_2to1 first.
  await mockScoping(page);
  await page.goto("/series/new");
  await expect(page.getByTestId("dock-member-count")).toHaveText(/1 \/ 2/);

  const rows = page.locator("[data-reorder-row]");
  // Initial: A_1to1 is first.
  await expect(page.getByTestId("scope-sample-row").first()).toContainText("A_1to1");

  // Drag from the first row to the second row.
  await rows.nth(0).dragTo(rows.nth(1));

  // After drag: B_2to1 moved to the first position.
  await expect(page.getByTestId("scope-sample-row").first()).toContainText("B_2to1");
});

// ---------------------------------------------------------------------------
// IS-5: ⌘Z undoes a keyboard reorder; ⌘⇧Z redoes it
// ---------------------------------------------------------------------------

test("IS-5: ⌘Z undoes a reorder and ⌘⇧Z redoes it", async ({ page }) => {
  await mockScoping(page);
  await page.goto("/series/new");
  await expect(page.getByTestId("dock-member-count")).toHaveText(/1 \/ 2/);

  const rows = page.getByTestId("scope-sample-row");
  // Initial: A_1to1 first.
  await expect(rows.first()).toContainText("A_1to1");

  // Reorder via the grip button (ScopeSampleRow's keyboard path): focus B_2to1's
  // grip button, then press ArrowUp. This pushes a "reorder" entry to the undo stack.
  const grip = page.getByRole("button", { name: /^Reorder B_2to1$/i });
  await grip.focus();
  await page.keyboard.press("ArrowUp");
  await expect(rows.first()).toContainText("B_2to1");

  // ⌘Z restores the original order (prevOrder=[10,11]).
  await page.keyboard.press("Meta+z");
  await expect(rows.first()).toContainText("A_1to1");

  // The Redo affordance appears after an undo.
  await expect(page.getByRole("button", { name: /redo/i })).toBeVisible();

  // ⌘⇧Z replays the reorder (nextOrder=[11,10]).
  await page.keyboard.press("Meta+Shift+z");
  await expect(rows.first()).toContainText("B_2to1");
});
