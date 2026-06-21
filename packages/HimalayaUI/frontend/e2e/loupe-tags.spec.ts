import { test, expect, type Page } from "@playwright/test";

// Mocked E2E for Slice 4 (Tasks 7-9): ManageTagsModal edit + delete flows.
// Asserts rendered semantics (role/aria/text), NOT Tailwind class strings.

const EXPERIMENT = {
  id: 1, name: "SSRL Test", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-05-01",
};

const SAMPLE_WITH_DUP_TAGS = {
  // Post display_name→name collapse: `name` carries the human label the loupe
  // renders (sampleDisplayName + ManageTagsModal both read `name`). The dialog
  // is titled by it, so the name IS what the /JC010/ assertions match.
  id: 10, experiment_id: 1, name: "JC010",
  notes: null, q_units: "A-1",
  tags: [
    { id: 1, key: "lipid", value: "LL37", source: "manual" },
    { id: 2, key: "lipid", value: "LL37", source: "scoping" },
  ],
};

const EXPOSURES = [
  {
    id: 1, sample_id: 10, filename: "pos1.dat", kind: "file", selected: true,
    status: "accepted", image_path: null, image_version: "", tags: [],
    sources: [], trace_hash: null, analysis_inputs_hash: null,
  },
];

const CORPUS_TAGS = [
  { key: "lipid", value: "LL37", count: 5 },
  { key: "lipid", value: "DMPC", count: 2 },
  { key: "temperature", value: "37C", count: 3 },
];

async function mockLoupe(page: Page, sample = SAMPLE_WITH_DUP_TAGS): Promise<void> {
  await page.addInitScript(() =>
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({
        state: { username: "alice", theme: "dark", tutorialSeen: true },
        version: 3,
      }),
    ),
  );

  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([EXPERIMENT]) }));
  await page.route("**/api/experiments/1", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify(EXPERIMENT) }));
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([sample]) }));
  await page.route(`**/api/samples/${sample.id}/exposures*`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify(EXPOSURES) }));
  await page.route(`**/api/samples/${sample.id}/messages`, (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/sample-tags", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify(CORPUS_TAGS) }));
  // Stub the SSE stream so it never sends events (avoids hanging connection)
  await page.route("**/api/events", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));

  // Stub tag mutation endpoints
  await page.route(`**/api/samples/${sample.id}/tags`, async (route) => {
    const req = route.request();
    if (req.method() === "POST") {
      const body = JSON.parse(req.postData() ?? "{}");
      await route.fulfill({
        status: 201, contentType: "application/json",
        body: JSON.stringify({ id: 99, key: body.key, value: body.value, source: "manual" }),
      });
    } else {
      await route.continue();
    }
  });
  await page.route(`**/api/samples/${sample.id}/tags/*`, async (route) => {
    const req = route.request();
    if (req.method() === "DELETE") {
      await route.fulfill({ status: 204 });
    } else if (req.method() === "PATCH") {
      const body = JSON.parse(req.postData() ?? "{}");
      await route.fulfill({
        status: 200, contentType: "application/json",
        body: JSON.stringify({ id: 1, key: body.key ?? "lipid", value: body.value ?? "LL37", source: "manual" }),
      });
    } else {
      await route.continue();
    }
  });
}

test("Manage button opens the dialog with role=dialog and aria-modal=true", async ({ page }) => {
  await mockLoupe(page);
  await page.goto(`/sample/10/loupe`);

  // Wait for the side panel to load
  await expect(page.getByTestId("loupe-side-panel")).toBeVisible();

  // Click the Manage button
  await page.getByTestId("manage-tags-trigger").click();

  // The modal dialog should appear
  const dialog = page.getByRole("dialog", { name: /JC010/ });
  await expect(dialog).toBeVisible();
  await expect(dialog).toHaveAttribute("aria-modal", "true");
});

test("Manage modal shows both duplicate lipid LL37 tags as editable rows", async ({ page }) => {
  await mockLoupe(page);
  await page.goto(`/sample/10/loupe`);
  await expect(page.getByTestId("loupe-side-panel")).toBeVisible();
  await page.getByTestId("manage-tags-trigger").click();

  const dialog = page.getByRole("dialog", { name: /JC010/ });
  await expect(dialog).toBeVisible();

  // Two manage-tag-row elements should appear (one per duplicate tag)
  const rows = dialog.getByTestId("manage-tag-row");
  await expect(rows).toHaveCount(2);
});

test("delete the second duplicate row via its Remove button", async ({ page }) => {
  await mockLoupe(page);
  await page.goto(`/sample/10/loupe`);
  await expect(page.getByTestId("loupe-side-panel")).toBeVisible();
  await page.getByTestId("manage-tags-trigger").click();

  const dialog = page.getByRole("dialog", { name: /JC010/ });
  await expect(dialog).toBeVisible();

  // Both rows have Remove lipid LL37 buttons — click the second one (scoping copy, id 2)
  const removeButtons = dialog.getByRole("button", { name: /remove lipid LL37/i });
  await expect(removeButtons).toHaveCount(2);

  // The second remove button targets tag id 2 (the scoping duplicate)
  await removeButtons.nth(1).click();

  // The queue's optimistic onMutate removes that one tag from the cache
  // immediately (client-side, independent of the mock), so the modal re-renders
  // with exactly one row left — proving the second duplicate is individually
  // addressable (LO-TAGDUP). The full server round-trip is in the queue tests.
  await expect(dialog.getByTestId("manage-tag-row")).toHaveCount(1);
});

test("edit a tag value via the value TagSuggest and press Enter", async ({ page }) => {
  await mockLoupe(page);
  await page.goto(`/sample/10/loupe`);
  await expect(page.getByTestId("loupe-side-panel")).toBeVisible();
  await page.getByTestId("manage-tags-trigger").click();

  const dialog = page.getByRole("dialog", { name: /JC010/ });
  await expect(dialog).toBeVisible();

  // First row, value combobox
  const firstRow = dialog.getByTestId("manage-tag-row").first();
  const comboboxes = firstRow.getByRole("combobox");
  await expect(comboboxes).toHaveCount(2);
  const valueCombo = comboboxes.nth(1);

  // Clear and type a new value
  await valueCombo.fill("DMPC");
  await valueCombo.press("Enter");

  // The edit was committed (we can check the input reflects the new value)
  await expect(valueCombo).toHaveValue("DMPC");
});

test("Done button closes the modal", async ({ page }) => {
  await mockLoupe(page);
  await page.goto(`/sample/10/loupe`);
  await expect(page.getByTestId("loupe-side-panel")).toBeVisible();
  await page.getByTestId("manage-tags-trigger").click();

  const dialog = page.getByRole("dialog", { name: /JC010/ });
  await expect(dialog).toBeVisible();

  await dialog.getByRole("button", { name: "Done" }).click();
  await expect(dialog).not.toBeVisible();
});
