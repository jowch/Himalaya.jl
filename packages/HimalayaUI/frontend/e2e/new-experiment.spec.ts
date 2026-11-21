import { test, expect, type Page } from "./fixtures";

const json = (body: unknown) => ({ status: 200, contentType: "application/json", body: JSON.stringify(body) });

const GEOMETRY = {
  energy_kev: 9, energy_kev_source: "prp",
  flight_path_m: 1.8095, flight_path_m_source: "setup",
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

async function mockFunnel(page: Page, created: unknown[]): Promise<void> {
  await page.addInitScript(() => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({ state: { username: "alice", tutorialSeen: true }, version: 5 }));
  });
  await page.route("**/api/users", (r) => r.fulfill(json([])));
  await page.route("**/api/events*", (r) => r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
  await page.route("**/api/sample-tags", (r) => r.fulfill(json([])));
  await page.route("**/api/tag-colors", (r) => r.fulfill(json([])));
  await page.route(/\/api\/fs\/suggest\?/, (r) => {
    const prefix = new URL(r.request().url()).searchParams.get("prefix") ?? "";
    const folders = ["/lab/20261103/notes", "/lab/20261103/setup1", "/lab/20261103/setup2"];
    return r.fulfill(json({ suggestions: folders.filter((f) => f.startsWith(prefix)) }));
  });
  await page.route(/\/api\/fs\/preview\?/, (r) => {
    const path = new URL(r.request().url()).searchParams.get("path");
    return r.fulfill(json(path === "/lab/20261103"
      ? { path, setups: [setup("setup1", 1.8095, 1220), setup("setup2", 0.5, 340)] }
      : { path, setups: [] }));
  });
  await page.route((url) => url.pathname === "/api/experiments", (r) => {
    if (r.request().method() === "POST") {
      const body = r.request().postDataJSON() as { name: string };
      created.push(body);
      return r.fulfill({ status: 202, contentType: "application/json", body: JSON.stringify({ id: created.length, name: body.name, status: "scanning" }) });
    }
    return r.fulfill(json([]));
  });
}

test("a beamtime folder adds every setup as its own experiment", async ({ page }) => {
  const created: { name: string; data_dir: string; geometry: { flight_path_m: number } }[] = [];
  await mockFunnel(page, created);
  await page.goto("/experiments/new");
  await page.getByTestId("dirpicker-input").locator("input").fill("/lab/20261103");
  const rows = page.getByTestId("setup-row");
  await expect(rows).toHaveCount(2);
  await expect(rows.nth(0).getByTestId("setup-summary")).toHaveText("1,220 exposures · flight path 1.81 m");
  await expect(rows.nth(1).getByTestId("setup-summary")).toHaveText("340 exposures · flight path 0.50 m");
  await page.getByTestId("dirpicker-input").locator("input").press("ControlOrMeta+Enter");
  await expect(page).toHaveURL(/\/experiments$/);
  expect(created.map((c) => [c.name, c.data_dir, c.geometry.flight_path_m])).toEqual([
    ["20261103/setup1", "/lab/20261103/setup1/data", 1.8095],
    ["20261103/setup2", "/lab/20261103/setup2/data", 0.5],
  ]);
});

test("⌘/Ctrl+Enter on a just-unticked setup adds only the ticked one and opens it", async ({ page }) => {
  const created: { name: string }[] = [];
  await mockFunnel(page, created);
  await page.route("**/api/experiments/1", (r) => r.fulfill(json({ id: 1, name: "20261103/setup1", ingest_status: "scanning" })));
  await page.route("**/api/experiments/1/**", (r) => r.fulfill(json([])));
  await page.route((url) => url.pathname === "/api/samples", (r) => r.fulfill(json([])));
  await page.goto("/experiments/new");
  await page.getByTestId("dirpicker-input").locator("input").fill("/lab/20261103");
  const second = page.getByTestId("setup-row").nth(1).getByRole("checkbox");
  await second.click();
  await expect(second).toHaveAttribute("aria-checked", "false");
  await page.keyboard.press("ControlOrMeta+Enter");
  await expect(page).toHaveURL(/\/experiments\/1/);
  expect(created.map((c) => c.name)).toEqual(["20261103/setup1"]);
});
