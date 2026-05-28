/**
 * Live-server spec for issue #35 Bug 4 — `addSampleTag` SSE-wins shape mismatch.
 *
 * Pre-fix, when the SSE frame for `add_tag` beat the HTTP body, the mutator's
 * `onSuccess` ran with the synth payload (which spreads only the SSE event
 * fields). The synth omitted `id` and `source`, so the cache row landed as
 * `{ id: undefined, key, value, source: undefined }`. Visible: the ×-delete
 * button on a freshly-added tag triggered `onRemoveTag(undefined)` and the
 * HTTP request hit `DELETE /api/samples/:id/tags/undefined` → 404, never
 * removing the row.
 *
 * Post-fix: `synthesizeResponseFromSse` is kind-aware for `add_tag` and maps
 * `tag_id → id` while pinning `source: "manual"`.
 *
 * The user flow is "add a tag, then delete it from the same UI." If the cache
 * carries a real id, deletion succeeds and the tag disappears.
 *
 * Repointed (#207): post-#163 the Inspect page is gone; the corpus contact
 * sheet (`/samples`) is the canonical home for sample-scoped tag editing. Each
 * row carries a `+ tag` invite + inline form, and every chip has a remove ✕.
 * The mutator code path (`addCorpusSampleTagMutator` → backend's `add_tag`
 * route) reuses the same synth shape the original bug class lived in, so this
 * spec still exercises the regression — through the corpus cache key instead
 * of the per-experiment one.
 */
import { test, expect, type Page } from "@playwright/test";

const BACKEND_BASE = process.env["BACKEND_BASE"] ?? "http://127.0.0.1:8090";

type Sample = { id: number; experiment_id: number; tags: { id: number; key: string }[] };

interface Fixture {
  experimentId: number;
  sampleId: number;
}

async function findFixture(): Promise<Fixture> {
  // Any sample works; we test the lifecycle in isolation. Pick the first.
  const exps = await fetch(`${BACKEND_BASE}/api/experiments`).then(r => r.json());
  for (const exp of exps) {
    const samples: Sample[] = await fetch(
      `${BACKEND_BASE}/api/experiments/${exp.id}/samples`).then(r => r.json());
    if (samples.length > 0) {
      return { experimentId: exp.id, sampleId: samples[0]!.id };
    }
  }
  throw new Error(`No sample at ${BACKEND_BASE}.`);
}

async function seedCorpusSheet(page: Page, fx: Fixture): Promise<void> {
  await page.addInitScript((args) => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: {
        username: args.username, firstName: args.username, lastName: "tester",
        // The corpus surface owns its own URL; we seed username only.
        tutorialSeen: true,
        theme: "dark",
      },
      version: 3,
    }));
  }, { username: "tag-tester" });
  // Filter the contact sheet to the fixture's beamtime so only the right rows
  // paint. The corpus surface still mounts every row globally, but a smaller
  // working set makes the row lookup deterministic.
  await page.goto(`/samples?beamtime=${fx.experimentId}`);
  // The row's identity (data-testid="sample-row-<id>") is the canonical
  // anchor for "the corpus sheet is ready and the row I want is in the DOM".
  await expect(page.getByTestId(`sample-row-${fx.sampleId}`)).toBeVisible();
}

async function getSampleTags(sampleId: number): Promise<{ id: number; key: string; value: string }[]> {
  const sample: Sample & { tags: { id: number; key: string; value: string }[] } =
    await fetch(`${BACKEND_BASE}/api/samples/${sampleId}`).then(r => r.json());
  return sample.tags;
}

test.describe("issue #35 sample-tag add-then-delete reconciliation (Bug 4)", () => {
  let fx: Fixture;

  test.beforeAll(async () => {
    fx = await findFixture();
    console.log(`[sample-tag] using experiment=${fx.experimentId} sample=${fx.sampleId}`);
  });

  test("add a tag, then delete it from the same UI", async ({ page }) => {
    await seedCorpusSheet(page, fx);

    // Unique-ish key so re-runs against a non-reset DB don't collide.
    const tagKey = `pr36-test-${Date.now()}`;
    const tagVal = "ok";

    const beforeIds = new Set((await getSampleTags(fx.sampleId)).map(t => t.id));

    // The contact-sheet `+ tag` button lives inside this row's tags cell.
    // (Empty-state rows show the full "+ tag" copy; rows with chips show "+".
    // Either way the testid is "tag-add".) The hover-only `opacity-0` variant
    // is still in the DOM + clickable; Playwright's click flow moves the
    // cursor into the row before dispatching, so no explicit hover is needed.
    const row = page.getByTestId(`sample-row-${fx.sampleId}`);
    await row.getByTestId("tag-add").click();
    const keyInput = row.getByPlaceholder("key");
    const valInput = row.getByPlaceholder("value");
    await expect(keyInput).toBeVisible();
    await keyInput.fill(tagKey);
    await valInput.fill(tagVal);
    await row.getByRole("button", { name: "Add" }).click();

    // The new tag chip should appear, with a × delete button identified by
    // aria-label "Remove ${key} tag".
    const removeBtn = row.getByRole("button", { name: `Remove ${tagKey} tag` });
    await expect(removeBtn).toBeVisible({ timeout: 8000 });

    // Backend sanity: poll for the tag landing — the optimistic update fires
    // before the HTTP response settles, so a single fetch can race ahead of
    // the server-side INSERT.
    let newTagId: number | undefined;
    await expect.poll(async () => {
      const tags = await getSampleTags(fx.sampleId);
      const fresh = tags.filter(t => !beforeIds.has(t.id) && t.key === tagKey);
      if (fresh.length === 1) newTagId = fresh[0]!.id;
      return fresh.length;
    }, { timeout: 8000 }).toBe(1);
    expect(newTagId!, "server assigned a real id").toBeGreaterThan(0);

    // The pre-fix bug: clicking × hit DELETE .../undefined → 404, the chip
    // stayed. Post-fix: the cache row carries the real id from the synth.
    await removeBtn.click();
    await expect(removeBtn).not.toBeVisible({ timeout: 8000 });

    // Backend sanity: poll for server-side removal.
    await expect.poll(async () => {
      const tags = await getSampleTags(fx.sampleId);
      return tags.some(t => t.id === newTagId);
    }, { timeout: 8000 }).toBe(false);
  });
});
