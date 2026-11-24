import { test, expect, type Page } from "./fixtures";

const json = (body: unknown) => ({ status: 200, contentType: "application/json", body: JSON.stringify(body) });

const GEOMETRY = {
  energy_kev: 9, energy_kev_source: "prp",
  flight_path_m: 1.81, flight_path_m_source: "setup",
  beam_center_x: 421.4, beam_center_x_source: "setup",
  beam_center_y: 836.9, beam_center_y_source: "setup",
  pixel_size_um: 172, pixel_size_um_source: "computed",
};

function setup(name: string, flight_path_m: number, exposures: number): object {
  const root = `/lab/20261103/${name}`;
  return {
    root, name: `20261103/${name}`, data_dir: `${root}/data`,
    analysis_dir: `${root}/analysis/automatic_analysis/tot_files`,
    setup_file: `${root}/analysis/setup_info_20261103_090000.txt`, setup_ambiguous: false,
    layout_detected: true,
    patterns: { image: "{name}_0_001.tif", metadata: "{name}_0_001.prp", integration: "{name}_tot.dat" },
    exposures, without_metadata: 0, without_trace: 0,
    geometry: { ...GEOMETRY, flight_path_m }, discrepancies: [], duplicate_of: null,
  };
}

interface CreatedSetup { root: string; data_dir: string; included: boolean; geometry?: { flight_path_m?: number } }
interface Created { name: string; path: string; setups: CreatedSetup[] }

async function mockFunnel(page: Page, created: Created[]): Promise<void> {
  await page.addInitScript(() => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({ state: { username: "alice", tutorialSeen: true }, version: 6 }));
  });
  await page.route("**/api/users", (r) => r.fulfill(json([])));
  await page.route("**/api/events*", (r) => r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
  await page.route("**/api/label-colors", (r) => r.fulfill(json([])));
  await page.route(/\/api\/fs\/suggest\?/, (r) => {
    const prefix = new URL(r.request().url()).searchParams.get("prefix") ?? "";
    const folders = [
      { path: "/lab/20261103", kind: "beamtime", setups: 2, added: null, described: true },
      { path: "/lab/20261103/notes", kind: null, setups: null, added: null, described: true },
      { path: "/lab/20261103/setup1", kind: "setup", setups: null, added: null, described: true },
      { path: "/lab/20261103/setup2", kind: "setup", setups: null, added: null, described: true },
    ];
    const dir = prefix.slice(0, prefix.lastIndexOf("/") + 1);
    const inDir = folders.filter((f) => f.path.lastIndexOf("/") + 1 === dir.length && f.path.startsWith(prefix));
    return r.fulfill(json({ folders: inDir, more: false }));
  });
  await page.route(/\/api\/fs\/preview\?/, (r) => {
    const path = (new URL(r.request().url()).searchParams.get("path") ?? "").replace(/\/+$/, "");
    return r.fulfill(json(path === "/lab/20261103"
      ? { path, existing: null, setups: [setup("setup1", 1.81, 1220), setup("setup2", 0.5, 340)] }
      : { path, existing: null, setups: [] }));
  });
  await page.route((url) => url.pathname === "/api/experiments", (r) => {
    if (r.request().method() === "POST") {
      const body = r.request().postDataJSON() as Created;
      created.push(body);
      return r.fulfill({ status: 202, contentType: "application/json", body: JSON.stringify({ id: 1, name: body.name, status: "scanning" }) });
    }
    return r.fulfill(json([]));
  });
  await page.route("**/api/experiments/1", (r) => r.fulfill(json({ id: 1, name: "20261103", path: "/lab/20261103", ingest_status: "scanning", setups: [] })));
  await page.route("**/api/experiments/1/**", (r) => r.fulfill(json([])));
  await page.route((url) => url.pathname === "/api/samples", (r) => r.fulfill(json([])));
}

test("a beamtime folder becomes one experiment holding every setup", async ({ page }) => {
  const created: Created[] = [];
  await mockFunnel(page, created);
  await page.goto("/experiments/new");
  await page.getByTestId("dirpicker-input").locator("input").fill("/lab/20261103");
  const rows = page.getByTestId("setup-row");
  await expect(rows).toHaveCount(2);
  await expect(rows.nth(0).getByTestId("setup-summary")).toHaveText("1,220 images · flight path 1.81 m");
  await expect(rows.nth(1).getByTestId("setup-summary")).toHaveText("340 images · flight path 0.50 m");
  await page.getByTestId("dirpicker-input").locator("input").press("ControlOrMeta+Enter");
  await expect(page).toHaveURL(/\/experiments\/1/);
  expect(created).toHaveLength(1);
  expect(created[0]!.name).toBe("20261103");
  expect(created[0]!.path).toBe("/lab/20261103");
  expect(created[0]!.setups.map((s) => [s.data_dir, s.included, s.geometry?.flight_path_m])).toEqual([
    ["/lab/20261103/setup1/data", true, 1.81],
    ["/lab/20261103/setup2/data", true, 0.5],
  ]);
});

test("⌘/Ctrl+Enter on a just-unticked setup sends it unticked and opens the experiment", async ({ page }) => {
  const created: Created[] = [];
  await mockFunnel(page, created);
  await page.goto("/experiments/new");
  await page.getByTestId("dirpicker-input").locator("input").fill("/lab/20261103");
  const second = page.getByTestId("setup-row").nth(1).getByRole("checkbox");
  await second.click();
  await expect(second).toHaveAttribute("aria-checked", "false");
  await page.keyboard.press("ControlOrMeta+Enter");
  await expect(page).toHaveURL(/\/experiments\/1/);
  expect(created).toHaveLength(1);
  expect(created[0]!.setups.map((s) => [s.root, s.included])).toEqual([
    ["/lab/20261103/setup1", true],
    ["/lab/20261103/setup2", false],
  ]);
});

test("typing a path shows the folders below it, and the keyboard picks one", async ({ page }) => {
  const created: Created[] = [];
  await mockFunnel(page, created);
  await page.goto("/experiments/new");
  const field = page.getByRole("combobox", { name: "Data directory" });
  await expect(field).toBeFocused();
  await page.keyboard.type("/lab/2026");
  const list = page.getByRole("listbox", { name: "Folders" });
  await expect(list.getByRole("option")).toHaveCount(1);
  await expect(list.getByRole("option", { name: "20261103, beamtime, 2 setups" })).toBeVisible();
  await page.keyboard.press("Tab");
  await page.keyboard.press("Enter");
  await expect(field).toHaveValue("/lab/20261103/");
  await expect(page.getByTestId("setup-row")).toHaveCount(2);
  await expect(list).toBeHidden();
  await page.keyboard.press("ArrowDown");
  await expect(list.getByRole("option")).toHaveCount(3);
  await expect(list.getByRole("option", { name: "notes" })).toHaveAttribute("aria-selected", "true");
  await page.keyboard.press("Escape");
  await expect(list).toBeHidden();
  await page.keyboard.press("Tab");
  await expect(page.getByRole("textbox", { name: "Experiment name" })).toBeFocused();
  await page.keyboard.press("ControlOrMeta+Enter");
  await expect(page).toHaveURL(/\/experiments\/1/);
  expect(created[0]!.path).toBe("/lab/20261103");
});
