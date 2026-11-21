import { test, expect } from "./fixtures";

const SAMPLES = ["JC001", "JC002"].map((name, i) => ({
  id: 40 + i, experiment_id: 17, display_name: name, name,
  notes: null, q_units: "A-1",
}));

const NAMES = ["buffer", "37C", "heated twice", "retake", "dilute 1 to 10", "cloudy", "after spin",
  "second capillary", "check flare", "low counts", "repeat tomorrow", "reference"];
const LABELS = NAMES.map((label, i) => ({
  id: i + 1, sample_id: 40, label, state: i === 0 ? "suggested" : "confirmed", source: "user",
}));

test.use({ viewport: { width: 1280, height: 800 } });

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() => {
    localStorage.clear();
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: { username: "alice", tutorialSeen: true },
      version: 5,
    }));
  });
  const json = (body: unknown) => ({ status: 200, contentType: "application/json", body: JSON.stringify(body) });
  const experiment = { id: 17, name: "lipid", path: "", data_dir: "", analysis_dir: "",
    manifest_path: null, created_at: "", q_units: null };
  await page.route("**/api/experiments", (r) => r.fulfill(json([experiment])));
  await page.route("**/api/experiments/17", (r) => r.fulfill(json(experiment)));
  await page.route((url) => url.pathname === "/api/samples", (r) => r.fulfill(json(SAMPLES)));
  await page.route("**/api/users", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/17/bridges", (r) => r.fulfill(json([])));
  await page.route("**/api/label-colors", (r) => r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments/17/scan", (r) =>
    r.fulfill(json({ status: "complete", experiment_id: 17, changed: false })));
  await page.route("**/api/experiments/17/labels", (r) => r.fulfill(json(LABELS)));
  await page.route("**/api/experiments/17/problems", (r) => r.fulfill(json({ ingest_error: null, exposures: [] })));
  await page.route(/\/api\/samples\/\d+\/call$/, (r) => {
    const id = Number(/samples\/(\d+)\//.exec(r.request().url())![1]);
    return r.fulfill(json({ sample_id: id, state: "indexed", phase: null, call_phases: [], members: [] }));
  });
  await page.route("**/api/samples/*/exposures*", (r) => r.fulfill(json([])));
  await page.route("**/api/events", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
});

async function belowLabels(page: import("@playwright/test").Page): Promise<number> {
  return page.getByTestId("workspace-labels").evaluate((el) => {
    const next = el.nextElementSibling as HTMLElement | null;
    return (next ?? el).getBoundingClientRect().top + (next ? 0 : el.getBoundingClientRect().height);
  });
}

test("many labels stay on one line, every label is reachable, and the content below never moves", async ({ page }) => {
  await page.goto("/experiments/17/workspace/41");
  await expect(page.getByRole("button", { name: "+ Label L" })).toBeVisible();
  const emptyTop = await belowLabels(page);

  await page.goto("/experiments/17/workspace/40");
  const chips = page.getByTestId("workspace-label-chips");
  const more = page.getByTestId("workspace-labels-more");
  await expect(more).toBeVisible();
  expect(await belowLabels(page)).toBe(emptyTop);

  const row = await page.getByTestId("workspace-labels").boundingBox();
  const tops = await chips.locator("[data-label]").evaluateAll((els) => els.map((e) => e.getBoundingClientRect().top));
  expect(new Set(tops).size).toBe(1);
  const rights = await chips.locator("[data-label]").evaluateAll((els) => els.map((e) => e.getBoundingClientRect().right));
  const box = await chips.boundingBox();
  for (const r of rights) expect(r).toBeLessThanOrEqual(box!.x + box!.width + 0.5);
  expect(row!.height).toBeLessThan(40);
  expect(await page.locator("main").evaluate((m) => m.scrollWidth - m.clientWidth)).toBeLessThanOrEqual(0);
  expect(await page.locator("[data-label]").count()).toBe(await chips.locator("[data-label]").count());

  const shown = await chips.locator("[data-label]").evaluateAll((els) => els.map((e) => e.getAttribute("data-label")));
  await more.click();
  const overflow = page.getByRole("dialog", { name: "More labels" });
  const hidden = await overflow.locator("[data-label]").evaluateAll((els) => els.map((e) => e.getAttribute("data-label")));
  expect([...shown, ...hidden]).toEqual(NAMES);
  await expect(more).toHaveText(`+${hidden.length}`);
  expect(await belowLabels(page)).toBe(emptyTop);
  await expect(overflow.getByRole("button").first()).toBeFocused();
  await page.keyboard.press("Escape");
  await expect(overflow).toBeHidden();
  await expect(more).toBeFocused();
});

test("L opens the label picker with an empty, focused field, and Enter confirms the suggested buffer", async ({ page }) => {
  const posts: unknown[] = [];
  await page.route(/\/api\/samples\/40\/labels$/, (r) => {
    posts.push(r.request().postDataJSON());
    return r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([{ id: 1, sample_id: 40, label: "Buffer", state: "confirmed", source: "user" }]) });
  });
  await page.goto("/experiments/17/workspace/40");
  await expect(page.getByRole("button", { name: "+ Label L" })).toBeVisible();
  await page.keyboard.press("l");
  const field = page.getByRole("dialog", { name: "Add a label" }).getByLabel("Label name");
  await expect(field).toBeFocused();
  await expect(field).toHaveValue("");
  await page.keyboard.type("buf");
  await page.keyboard.press("Enter");
  await expect(page.getByRole("dialog", { name: "Add a label" })).toBeHidden();
  await expect.poll(() => posts.length).toBe(1);
  expect(posts[0]).toMatchObject({ label: "Buffer" });
});

test("labels share the title row, and the picker opens under + Label", async ({ page }) => {
  await page.goto("/experiments/17/workspace/41");
  const title = (await page.getByRole("heading", { level: 2, name: "JC002" }).boundingBox())!;
  const row = (await page.getByTestId("workspace-labels").boundingBox())!;
  expect(row.y + row.height / 2).toBeGreaterThan(title.y);
  expect(row.y + row.height / 2).toBeLessThan(title.y + title.height);
  expect(row.x).toBeGreaterThan(title.x + title.width);
  const add = (await page.getByRole("button", { name: "+ Label L" }).boundingBox())!;
  await page.keyboard.press("l");
  const picker = (await page.getByRole("dialog", { name: "Add a label" }).boundingBox())!;
  expect(Math.abs(picker.x + picker.width - (add.x + add.width))).toBeLessThan(2);
  expect(picker.y).toBeGreaterThanOrEqual(add.y + add.height);
});

test("Tab and Shift+Tab move the picker's highlight, focus stays in the field, and Enter adds the highlighted row", async ({ page }) => {
  const posts: unknown[] = [];
  await page.route(/\/api\/samples\/41\/labels$/, (r) => {
    posts.push(r.request().postDataJSON());
    return r.fulfill({ status: 200, contentType: "application/json", body: "[]" });
  });
  await page.goto("/experiments/17/workspace/41");
  await expect(page.getByRole("button", { name: "+ Label L" })).toBeVisible();
  await page.keyboard.press("l");
  const picker = page.getByRole("dialog", { name: "Add a label" });
  const field = picker.getByLabel("Label name");
  await field.fill("ca");
  const highlighted = picker.locator("[aria-current=true]");
  const options = await picker.getByRole("button").allTextContents();
  expect(options.length).toBeGreaterThan(1);
  await expect(highlighted).toHaveText(options[0]!);
  await page.keyboard.press("Tab");
  await expect(field).toBeFocused();
  await expect(highlighted).toHaveText(options[1]!);
  await page.keyboard.press("Shift+Tab");
  await page.keyboard.press("Shift+Tab");
  await expect(field).toBeFocused();
  await expect(highlighted).toHaveText(options.at(-1)!);
  await page.keyboard.press("Tab");
  await page.keyboard.press("Tab");
  await page.keyboard.press("Enter");
  await expect(picker).toBeHidden();
  await expect.poll(() => posts.length).toBe(1);
  expect(posts[0]).toMatchObject({ label: options[1] });
});

test("on a narrow window labels fold into +N before the sample name truncates", async ({ page }) => {
  await page.setViewportSize({ width: 1024, height: 700 });
  const json = (body: unknown) => ({ status: 200, contentType: "application/json", body: JSON.stringify(body) });
  const name = "buffer_blank_S2424 (S02P05)";
  await page.route((url) => url.pathname === "/api/samples", (r) =>
    r.fulfill(json([{ id: 42, experiment_id: 17, display_name: name, name, notes: null, q_units: "A-1" }])));
  await page.route("**/api/experiments/17/labels", (r) => r.fulfill(json(
    ["Buffer", "37C", "heated twice", "retake"].map((label, i) => ({ id: i + 1, sample_id: 42, label, state: "confirmed", source: "user" })))));
  await page.goto("/experiments/17/workspace/42");
  const title = page.getByRole("heading", { level: 2, name });
  await expect(title).toBeVisible();
  await expect(page.getByTestId("workspace-labels-more")).toBeVisible();
  expect(await title.evaluate((el) => el.scrollWidth <= el.clientWidth)).toBe(true);
  await expect(page.getByRole("button", { name: "+ Label L" })).toBeVisible();
});
