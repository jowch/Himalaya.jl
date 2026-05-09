import { test, expect } from "@playwright/test";

// Spec §8.3 — Live integration. Requires backend (port 8090) + Vite (5180)
// per packages/HimalayaUI/frontend/e2e/live/README.md.

test.describe("permalinks — live", () => {
  test("delete sample in tab B → tab A URL invalidates to StaleUrlPage", async ({ browser }) => {
    const ctxA = await browser.newContext();
    const ctxB = await browser.newContext();
    const tabA = await ctxA.newPage();
    const tabB = await ctxB.newPage();

    // Tab A opens an inspect deep link for an existing sample.
    await tabA.goto("/inspect/E2ETest/SeedSample");
    await tabA.waitForTimeout(800);  // SSE handshake (per gotcha)

    // Tab B deletes the sample via direct API.
    await tabB.request.delete("/api/samples/SEED_SAMPLE_ID", {
      headers: { "X-Username": "alice" },
    });

    // Tab A's URL invalidates; StaleUrlPage renders.
    await expect(tabA.locator("[data-testid='stale-url-page']")).toBeVisible({ timeout: 5_000 });
  });

  test("paste URL referencing a deleted sample lands on StaleUrlPage", async ({ page }) => {
    await page.goto("/inspect/E2ETest/AlreadyDeletedSample");
    await expect(page.locator("[data-testid='stale-url-page']")).toBeVisible();
  });

  test("same-user-different-tab: tab A sees tab B's delete (no client_id self-echo)", async ({ browser }) => {
    // Spec §7 — URL invalidation runs on every SSE event regardless of
    // originating client_id, so two tabs of the same username refresh each
    // other's URL state.
    const ctx = await browser.newContext();
    const tabA = await ctx.newPage();
    const tabB = await ctx.newPage();
    // Same browser context = shared sessionStorage for clientId? No — clientId
    // is per-tab via sessionStorage, so the two tabs DO have distinct client_ids
    // even within a single context. The username header is the same on both.
    await tabA.goto("/inspect/E2ETest/SeedSample");
    await tabA.waitForTimeout(800);

    await tabB.request.delete("/api/samples/SEED_SAMPLE_ID", {
      headers: { "X-Username": "alice" },
    });

    await expect(tabA.locator("[data-testid='stale-url-page']")).toBeVisible({ timeout: 5_000 });
  });
});
