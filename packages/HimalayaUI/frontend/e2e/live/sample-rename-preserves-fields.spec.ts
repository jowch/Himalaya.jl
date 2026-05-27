/**
 * Live-server spec for issue #35 Bug 5 — `updateSample` SSE-wins clobbers
 * unpatched fields.
 *
 * Pre-fix, when the SSE frame for `update_sample` beat the HTTP body, the
 * mutator's `onSuccess` ran with the synth payload. The synth spread the
 * SSE event's `payload` onto the cache row directly, but the SSE event for
 * `update_sample` carries only the diff (only the fields that were patched).
 * So if the user patched only `name`, the synth produced `{ name }` and the
 * cache row's `notes` and `label` got overwritten to `undefined`. Visible:
 * setting a sample's notes, then renaming the sample, wiped the notes.
 *
 * Post-fix: `updateSampleMutator.onSuccess` builds a partial patch from
 * defined response fields only — `Object.entries(response).filter(([, v]) =>
 * v !== undefined)` — so absent SSE keys leave the cache untouched.
 *
 * The user flow is "set notes, rename, verify notes still there." Pre-fix
 * this fails on at least the SSE-wins path; post-fix it always succeeds.
 *
 * Lives on the Inspect page (SampleMetadataCard).
 */
import { test, expect, type Page } from "@playwright/test";

const BACKEND_BASE = process.env["BACKEND_BASE"] ?? "http://127.0.0.1:8090";

type Sample = {
  id: number; experiment_id: number;
  name: string | null; display_name: string | null; notes: string | null;
};

interface Fixture {
  experimentId: number;
  sampleId: number;
  originalName: string | null;
  originalNotes: string | null;
}

async function findFixture(): Promise<Fixture> {
  const exps = await fetch(`${BACKEND_BASE}/api/experiments`).then(r => r.json());
  for (const exp of exps) {
    const samples: Sample[] = await fetch(
      `${BACKEND_BASE}/api/experiments/${exp.id}/samples`).then(r => r.json());
    if (samples.length > 0) {
      const s = samples[0]!;
      return {
        experimentId: exp.id, sampleId: s.id,
        originalName: s.name, originalNotes: s.notes,
      };
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
        tutorialSeen: true,
        theme: "dark",
      },
      version: 3,
    }));
  }, { username: "rename-tester", expId: fx.experimentId, sampId: fx.sampleId });
  await page.goto("/");
  // SampleMetadataCard's name input is the canonical anchor for "Inspect
  // page is ready".
  await expect(page.getByTestId("sample-name-input")).toBeVisible();
}

async function getSample(sampleId: number): Promise<Sample> {
  return fetch(`${BACKEND_BASE}/api/samples/${sampleId}`).then(r => r.json());
}

test.describe("issue #35 sample partial-patch reconciliation (Bug 5)", () => {
  let fx: Fixture;

  test.beforeAll(async () => {
    fx = await findFixture();
    console.log(
      `[sample-rename] using experiment=${fx.experimentId} sample=${fx.sampleId}`);
  });

  test.afterAll(async () => {
    // Restore the sample's display_name + notes so the dev DB doesn't drift between runs.
    await fetch(`${BACKEND_BASE}/api/samples/${fx.sampleId}`, {
      method: "PATCH",
      headers: { "Content-Type": "application/json", "X-Username": "rename-tester" },
      body: JSON.stringify({ display_name: fx.originalName, notes: fx.originalNotes }),
    }).catch(() => { /* best-effort cleanup */ });
  });

  test("renaming a sample preserves notes set in the same session", async ({ page }) => {
    await seedInspectPage(page, fx);

    // 1. Set a unique notes string and blur — triggers a partial PATCH {notes}.
    const notesText = `pr36-bug5-notes-${Date.now()}`;
    const notesArea = page.getByTestId("sample-notes-textarea");
    await notesArea.fill(notesText);
    await notesArea.blur();

    // Wait for the notes write to land on the server.
    await expect.poll(
      async () => (await getSample(fx.sampleId)).notes,
      { timeout: 8000 },
    ).toBe(notesText);

    // 2. Rename the sample (name-only PATCH). Pre-fix, the SSE-wins synth
    //    spread `{ name }` over the row, wiping notes to undefined; the
    //    optimistic cache view in the UI showed empty notes. Post-fix the
    //    partial-merge keeps notes intact.
    const newName = `pr36-bug5-name-${Date.now()}`;
    const nameInput = page.getByTestId("sample-name-input");
    await nameInput.fill(newName);
    await nameInput.blur();

    // 3. Wait for the display_name write to land.
    await expect.poll(
      async () => (await getSample(fx.sampleId)).display_name,
      { timeout: 8000 },
    ).toBe(newName);

    // 4. The notes textarea should still hold our string. Re-syncing the
    //    field on focus-out is gated by `!notesFocused`, so the post-fix
    //    cache value flows back into the input on next render.
    await expect(notesArea).toHaveValue(notesText);

    // 5. Backend sanity: server still has the notes too (full round-trip).
    const final = await getSample(fx.sampleId);
    expect(final.display_name).toBe(newName);
    expect(final.notes).toBe(notesText);
  });
});
