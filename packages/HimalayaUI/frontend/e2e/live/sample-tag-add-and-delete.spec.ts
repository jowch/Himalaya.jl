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
 * Lives on the Inspect page (SampleMetadataCard). Sample-scoped tags exercise
 * the same mutator code path as exposure-scoped tags (the bug class is the
 * synth shape, not the entity), so testing one suffices.
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

async function seedInspectPage(page: Page, fx: Fixture): Promise<void> {
  await page.addInitScript((args) => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: {
        username: args.username, firstName: args.username, lastName: "tester",
        activeExperimentId: args.expId,
        activeSampleId: args.sampId,
        activePage: "inspect",
        tutorialSeen: true,
        theme: "dark",
      },
      version: 3,
    }));
  }, { username: "tag-tester", expId: fx.experimentId, sampId: fx.sampleId });
  await page.goto("/");
  // SampleMetadataCard renders the "+ tag" affordance once the sample loads.
  await expect(page.getByRole("button", { name: "+ tag" })).toBeVisible();
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
    await seedInspectPage(page, fx);

    // Unique-ish key so re-runs against a non-reset DB don't collide.
    const tagKey = `pr36-test-${Date.now()}`;
    const tagVal = "ok";

    const beforeIds = new Set((await getSampleTags(fx.sampleId)).map(t => t.id));

    // Open the inline tag form.
    await page.getByRole("button", { name: "+ tag" }).click();
    const keyInput = page.getByPlaceholder("key");
    const valInput = page.getByPlaceholder("value");
    await expect(keyInput).toBeVisible();
    await keyInput.fill(tagKey);
    await valInput.fill(tagVal);
    await page.getByRole("button", { name: "Add" }).click();

    // The new tag chip should appear, with a × delete button identified by
    // aria-label "Remove ${key} tag".
    const removeBtn = page.getByRole("button", { name: `Remove ${tagKey} tag` });
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
