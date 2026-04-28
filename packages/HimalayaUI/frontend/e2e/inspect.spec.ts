import { test, expect, type Page } from "@playwright/test";

const EXPERIMENT = {
  id: 1, name: "SSRL Test", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-05-01",
};
const SAMPLES = [
  { id: 10, experiment_id: 1, label: "D1", name: "cubic_run03", notes: null, tags: [] },
];
const EXPOSURES = [
  {
    id: 1, sample_id: 10, filename: "pos1.dat", kind: "file",
    selected: true, status: "accepted", image_path: "/tmp/pos1.tiff",
    image_version: "v1-1700000001", tags: [], sources: [],
  },
  {
    id: 2, sample_id: 10, filename: "pos2.dat", kind: "file",
    selected: false, status: "rejected", image_path: null,
    image_version: "",
    tags: [{ id: 99, key: "rejection_reason", value: "flare", source: "manual" }],
    sources: [],
  },
  {
    id: 3, sample_id: 10, filename: "pos3.dat", kind: "file",
    selected: false, status: null, image_path: null,
    image_version: "", tags: [], sources: [],
  },
];

async function mockCore(page: Page): Promise<void> {
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([{ id: 1, username: "alice" }]) }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([EXPERIMENT]) }));
  await page.route("**/api/experiments/1", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify(EXPERIMENT) }));
  await page.route("**/api/experiments/1/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify(SAMPLES) }));
  await page.route("**/api/samples/10/exposures**", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify(EXPOSURES) }));
  await page.route("**/api/samples/10/messages", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([
        { id: 1, sample_id: 10, author_id: 1, author: "alice",
          body: "shared sample note", created_at: "2026-05-01T12:00:00Z" },
      ]) }));
  await page.route("**/api/exposures/*/image**", (r) =>
    r.fulfill({ status: 200, contentType: "image/png", body: Buffer.alloc(100) }));
  await page.route("**/api/exposures/*/status", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify({ id: 1, status: "accepted" }) }));
  await page.route("**/api/exposures/*/select", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify({ id: 1, selected: true }) }));
  await page.route("**/api/exposures/*/tags", (r) =>
    r.fulfill({ status: 201, contentType: "application/json",
      body: JSON.stringify({ id: 100, key: "rejection_reason", value: "test", source: "manual" }) }));
}

async function seedState(page: Page): Promise<void> {
  await page.addInitScript((state) => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({ state, version: 3 }));
  }, {
    username: "alice",
    activePage: "inspect",
    activeExperimentId: 1,
    activeSampleId: 10,
    tutorialSeen: true,
    theme: "dark",
  });
}

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() => { localStorage.clear(); });
});

test("Inspect tab appears before Index tab", async ({ page }) => {
  await mockCore(page);
  await page.addInitScript(() => {
    localStorage.setItem("himalaya-ui:state",
      JSON.stringify({ state: { username: "alice", activePage: "index", activeExperimentId: 1, activeSampleId: 10, tutorialSeen: true, theme: "dark" }, version: 3 }));
  });
  await page.goto("/");
  const tabs = page.getByRole("tab");
  await expect(tabs.first()).toHaveText("Inspect");
  await expect(tabs.nth(1)).toHaveText("Index");
});

test("navigating to Inspect shows the page", async ({ page }) => {
  await mockCore(page);
  await seedState(page);
  await page.goto("/");
  await expect(page.getByTestId("inspect-page")).toBeVisible();
});

test("rejected exposure cell is dimmed", async ({ page }) => {
  await mockCore(page);
  await seedState(page);
  await page.goto("/");
  await expect(page.getByTestId("inspect-page")).toBeVisible();
  const rejectedCell = page.getByTestId("thumb-cell-2");
  await expect(rejectedCell).toHaveAttribute("data-rejected", "true");
});

test("Reject → Other reveals custom reason input inline", async ({ page }) => {
  await mockCore(page);
  await seedState(page);
  await page.goto("/");
  await expect(page.getByTestId("inspect-page")).toBeVisible();
  // Click Reject — quick-pick chips appear (Flare, Other), no input yet.
  await page.getByRole("button", { name: /reject/i }).first().click();
  await expect(page.getByPlaceholder(/reason/i)).toHaveCount(0);
  // Click Other to reveal the freeform input.
  await page.getByRole("button", { name: /other/i }).click();
  await expect(page.getByPlaceholder(/reason/i)).toBeVisible();
});

test("clicking thumbnail updates the image card", async ({ page }) => {
  await mockCore(page);
  await seedState(page);
  await page.goto("/");
  await expect(page.getByTestId("inspect-page")).toBeVisible();
  await page.getByTestId("thumb-cell-3").click();
  await expect(page.getByText("pos3.dat")).toBeVisible();
});

test("ChatCard appears on Inspect with sample-scoped message history", async ({ page }) => {
  // Chat is sample-scoped, so the same conversation visible on Index appears
  // here too — gives curators continuity when switching pages.
  await mockCore(page);
  await seedState(page);
  await page.goto("/");
  await expect(page.getByTestId("inspect-page")).toBeVisible();
  const chat = page.getByTestId("chat-card");
  await expect(chat).toBeVisible();
  await expect(chat.getByText("shared sample note")).toBeVisible();
});

test("chat message persists when switching Index → Inspect", async ({ page }) => {
  await mockCore(page);
  await page.addInitScript(() => {
    localStorage.setItem("himalaya-ui:state",
      JSON.stringify({ state: { username: "alice", activePage: "index",
        activeExperimentId: 1, activeSampleId: 10, tutorialSeen: true,
        theme: "dark" }, version: 3 }));
  });
  await page.goto("/");
  // Verify the message is visible on Index
  await expect(page.getByText("shared sample note")).toBeVisible();
  // Switch to Inspect tab
  await page.getByRole("tab", { name: "Inspect" }).click();
  await expect(page.getByTestId("inspect-page")).toBeVisible();
  // Same message still there (same sampleId, same query cache, same component)
  await expect(page.getByText("shared sample note")).toBeVisible();
});
