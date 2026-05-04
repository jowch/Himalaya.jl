/**
 * Live-server Playwright spec for the four #37 bug fixes bundled into PR #36.
 *
 * Lives under `e2e/live/` (excluded from the default `npm run e2e`'s mocked
 * suite via testIgnore in `playwright.config.ts`). Run via `npm run e2e:live`,
 * which uses the dedicated `playwright.live.config.ts`.
 *
 * Runs against a real backend + dev DB — NO route mocks. Discovers an exposure
 * dynamically so the spec is portable across DBs.
 *
 * What this spec covers
 * ─────────────────────
 *   Bug 1a (single-tab, HTTP-wins, new custom group) — implicit in Test 1
 *     via "second confirmation works after the first." Pre-fix: cache stays
 *     on the auto group; second click would either silently no-op or write
 *     to a stale id. Post-fix: invalidate forces fresh fetch, second op
 *     targets the actual custom group.
 *   Bug 1c (foreign-tab, new custom group) — Test 2: two browser contexts
 *     on the same exposure; confirm in tab A; tab B sees the index in
 *     Active set within ~5s. Pre-fix: tab B silent until manual refresh.
 *
 * What this spec does NOT cover (covered by unit tests only)
 * ──────────────────────────────────────────────────────────
 *   Bug 1b (own-op SSE-wins) — Verifying the synth path requires forcing
 *     SSE to win the race against HTTP, which is timing-dependent and
 *     not reliably reproducible end-to-end. The unit tests in
 *     `test/queue/mutatorOnSseWins.test.ts` exercise the synth shape and
 *     onSuccess invalidation deterministically.
 *   Bug 2 (peak_excluded merge-not-replace) — User-visible impact is
 *     "limited" per the issue: optimistic state already has the right
 *     `excluded` flag, so the cache-row merge fix is conceptually correct
 *     but the UI looks the same pre-/post-fix. Verified by
 *     `cache-shape.test.ts` + `mutatorOnSseWins.test.ts`.
 *   Bug 3 (re-run loop throw safety) — Requires injecting a thrown
 *     onMutate, not reachable via real UI flows. Verified by the unit
 *     test "Re-run loop is throw-safe" in `replayCoordinator.test.ts`.
 *   Bug 4 (OP_LOCKS leak on body throw) — Server-side memory leak, not
 *     user-visible. Verified by the unit test in `test_idempotency.jl`.
 *
 * Setup
 * ─────
 *   1. Bring up the backend on a free port:
 *        cd <worktree>
 *        bin/himalaya serve <experiment-dir> --port 8080  # or whatever port
 *      Or without sysimage:
 *        julia --project=packages/HimalayaUI -e 'using HimalayaUI; main(ARGS)' \
 *          -- serve <experiment-dir> --port 8080
 *   2. Bring up Vite pointing at the backend:
 *        cd packages/HimalayaUI/frontend
 *        VITE_API_PORT=8080 npm run dev -- --host 127.0.0.1
 *   3. Run this spec:
 *        cd packages/HimalayaUI/frontend
 *        npm run e2e:live
 *
 * Pre-conditions on the dev DB
 * ────────────────────────────
 *   - At least one experiment with at least one analyzed exposure.
 *   - The chosen exposure should have ≥2 candidate indices (pre-fix Test 1
 *     needs index A and a different index B for the second confirmation).
 *   - The chosen exposure should NOT already have a custom group (i.e.
 *     a previous test run created one). To reset: open a fresh DB OR delete
 *     the custom group:
 *       DELETE FROM index_groups WHERE exposure_id = <id> AND kind = 'custom';
 *       UPDATE index_groups SET active = 1 WHERE exposure_id = <id> AND kind = 'auto';
 *
 * Pass / fail interpretation
 * ──────────────────────────
 *   - Test 1 ("first confirmation + second confirmation"): primary signal.
 *     If both indices land in the Active set AND the backend's groups
 *     query confirms a custom group with both, the cache invalidation in
 *     `addIndexToGroupMutator.onSuccess` is wired correctly.
 *   - Test 2 ("two-tab convergence"): primary signal for Bug 1c. If tab B
 *     sees the confirmed index move to Active set within the timeout, the
 *     `applyRemoteToCache` invalidation for index_confirmed is wired
 *     correctly.
 *   - The console-error check at the bottom catches any unexpected
 *     exceptions during the run (e.g. Bug 3 regressions: an onMutate
 *     throw aborting the re-run loop would log here).
 */
import { test, expect, type Page } from "@playwright/test";

const BACKEND_BASE = process.env["BACKEND_BASE"] ?? "http://127.0.0.1:8080";

type Exposure = { id: number; sample_id: number };
type Sample = { id: number; experiment_id: number };
type IndexEntry = { id: number; kind: string };
type GroupEntry = { id: number; exposure_id: number; kind: string;
                    active: boolean; members: number[] };

async function findTestExposure(): Promise<{ exposureId: number; indices: number[] }> {
  // Walk the API to find an analyzed exposure with ≥2 auto indices and
  // no existing custom group (so we exercise the first-confirmation path).
  const exps = await fetch(`${BACKEND_BASE}/api/experiments`).then(r => r.json());
  for (const exp of exps) {
    const samples: Sample[] = await fetch(
      `${BACKEND_BASE}/api/experiments/${exp.id}/samples`).then(r => r.json());
    for (const sm of samples) {
      const exposures: Exposure[] = await fetch(
        `${BACKEND_BASE}/api/samples/${sm.id}/exposures`).then(r => r.json());
      for (const ex of exposures) {
        const indices: IndexEntry[] = await fetch(
          `${BACKEND_BASE}/api/exposures/${ex.id}/indices`).then(r => r.json());
        const autoIxs = indices.filter(ix => ix.kind === "auto");
        if (autoIxs.length < 2) continue;
        const groups: GroupEntry[] = await fetch(
          `${BACKEND_BASE}/api/exposures/${ex.id}/groups`).then(r => r.json());
        const hasCustom = groups.some(g => g.kind === "custom");
        if (hasCustom) continue;
        return {
          exposureId: ex.id,
          indices: autoIxs.slice(0, 2).map(ix => ix.id),
        };
      }
    }
  }
  throw new Error(
    `No suitable exposure found at ${BACKEND_BASE}. ` +
    `Need an analyzed exposure with ≥2 auto indices and no custom group. ` +
    `See "Pre-conditions on the dev DB" in this file.`);
}

async function navigateToExposure(page: Page, exposureId: number): Promise<void> {
  // The frontend doesn't expose a deep-link to an exposure id; we click
  // through the navigator instead. If the app has an OnboardingFlow, dismiss
  // it first (the test DB may have no users registered).
  await page.goto("/");
  // Try to dismiss onboarding if present (no-op if absent).
  const onboardingCta = page.getByRole("button", { name: /^(continue|skip)/i });
  if (await onboardingCta.isVisible().catch(() => false)) {
    await onboardingCta.click().catch(() => {});
  }
  // Click the exposure tile by data attribute.
  await page.locator(`[data-exposure-id="${exposureId}"]`).first().click();
  // Wait for the index list to render (a card with [data-index-id] appears).
  await expect(page.locator("[data-index-id]").first()).toBeVisible();
}

async function confirmIndex(page: Page, indexId: number): Promise<void> {
  await page.getByRole("button", { name: `Add index ${indexId}` }).click();
  // Wait for the card to flip to active. data-active is set by the optimistic
  // onMutate, so this only confirms the optimistic update fired — assert on
  // the post-settle state separately.
  await expect(page.locator(`[data-index-id="${indexId}"][data-active]`))
    .toBeVisible();
}

test.describe("post-PR#33 extended fixes (issue #37)", () => {
  let exposureId: number;
  let indexA: number;
  let indexB: number;

  test.beforeAll(async () => {
    const t = await findTestExposure();
    exposureId = t.exposureId;
    [indexA, indexB] = t.indices as [number, number];
    console.log(
      `[post-pr37] using exposure=${exposureId} indexA=${indexA} indexB=${indexB}`);
  });

  test("Bug 1a — second confirmation works after first creates the custom group",
    async ({ page }) => {
      await navigateToExposure(page, exposureId);

      // First confirmation: triggers ensure_custom_group! on the backend.
      // Pre-fix: the new custom group's id wasn't in the cache, so the
      // mutator's onSuccess silently no-op'd and the auto group stayed
      // marked active client-side.
      await confirmIndex(page, indexA);
      await page.waitForTimeout(500);  // allow settle

      // Second confirmation: under the post-fix invalidate-then-refetch
      // contract, the cache now has the canonical {auto inactive, custom
      // active} structure, so this op targets the right group.
      await confirmIndex(page, indexB);
      await page.waitForTimeout(500);

      // Both indices should be in the Active set.
      await expect(page.locator(`[data-index-id="${indexA}"][data-active]`))
        .toBeVisible();
      await expect(page.locator(`[data-index-id="${indexB}"][data-active]`))
        .toBeVisible();

      // Backend sanity: a custom group exists with both indices in members.
      const groups: GroupEntry[] = await page.evaluate(async (exId) => {
        return (await fetch(`/api/exposures/${exId}/groups`)).json();
      }, exposureId);
      const custom = groups.find(g => g.kind === "custom");
      expect(custom, "custom group should exist post-confirmation").toBeDefined();
      expect(custom!.active).toBe(true);
      expect(custom!.members.sort()).toEqual([indexA, indexB].sort());
      const auto = groups.find(g => g.kind === "auto");
      expect(auto?.active, "auto group should be demoted").toBe(false);
    });

  test("Bug 1c — foreign tab picks up new custom group via SSE invalidation",
    async ({ browser }) => {
      // Two independent browser contexts simulate two tabs / two users.
      const ctxA = await browser.newContext();
      const ctxB = await browser.newContext();
      const tabA = await ctxA.newPage();
      const tabB = await ctxB.newPage();

      try {
        await navigateToExposure(tabA, exposureId);
        await navigateToExposure(tabB, exposureId);

        // Snapshot initial state in tab B.
        await expect(tabB.locator(`[data-index-id="${indexA}"]:not([data-active])`))
          .toBeVisible();

        // Tab A confirms.
        await confirmIndex(tabA, indexA);

        // Tab B should converge: the SSE for index_confirmed arrives, the
        // applyRemoteToCache fix detects the new custom group_id is not in
        // tab B's cache, invalidates, refetches, and re-renders the index
        // in the Active set. Pre-fix: tab B's surgical splice silently
        // missed the new group_id, so this assertion would time out.
        await expect(
          tabB.locator(`[data-index-id="${indexA}"][data-active]`),
          "tab B should see index in Active set after tab A confirmed",
        ).toBeVisible({ timeout: 8000 });
      } finally {
        await ctxA.close();
        await ctxB.close();
      }
    });

  test("no console errors during reconciliation", async ({ page }) => {
    const errors: string[] = [];
    page.on("pageerror", (e) => errors.push(`pageerror: ${e.message}`));
    page.on("console", (m) => {
      if (m.type() === "error") errors.push(`console.error: ${m.text()}`);
    });
    await navigateToExposure(page, exposureId);
    await confirmIndex(page, indexA);
    await page.waitForTimeout(500);
    // Bug 3 (re-run loop throw safety): if a mutator's onMutate threw and
    // the loop wasn't try/caught, the unhandled rejection would surface
    // here. The dev console should be clean.
    expect(errors, `unexpected errors: ${errors.join(" | ")}`).toHaveLength(0);
  });
});
