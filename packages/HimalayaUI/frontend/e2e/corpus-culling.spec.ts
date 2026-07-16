import { test, expect, type Page } from "@playwright/test";

// I1.7 (#163): the Phase-1 mocked spec for the corpus culling flow, ceded here
// by I1.6. Exercises the contact sheet (cull / batch-reject / representative)
// and the loupe (flip), all against mocked /api routes. Replaces the deleted
// Inspect E2E (inspect.spec.ts).

const EXPERIMENT = {
  id: 1, name: "SSRL Test", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-05-01",
};
const SAMPLES = [
  {
    id: 10, experiment_id: 1, display_name: "D1", name: "run03",
    notes: null, tags: [], q_units: "A-1",
  },
];
const EXPOSURES = [
  {
    id: 1, sample_id: 10, filename: "pos1.dat", kind: "file", selected: true,
    status: "accepted", image_path: null, image_version: "", tags: [],
    sources: [], trace_hash: null, analysis_inputs_hash: null,
  },
  {
    id: 2, sample_id: 10, filename: "pos2.dat", kind: "file", selected: false,
    status: "accepted", image_path: null, image_version: "", tags: [],
    sources: [], trace_hash: null, analysis_inputs_hash: null,
  },
  {
    id: 3, sample_id: 10, filename: "pos3.dat", kind: "file", selected: false,
    status: "accepted", image_path: null, image_version: "", tags: [],
    sources: [], trace_hash: null, analysis_inputs_hash: null,
  },
];

async function mockCorpus(page: Page): Promise<void> {
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([EXPERIMENT]) }));
  await page.route("**/api/experiments/1", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXPERIMENT) }));
  await page.route("**/api/experiments/1/loads", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(SAMPLES) }));
  await page.route("**/api/samples/10/exposures*", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXPOSURES) }));
  await page.route("**/api/samples/10/messages", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  // The loupe's Manage-tags modal reads the corpus tag vocabulary; an unstubbed
  // request here would leave the query in-flight and swallow the page's keyboard
  // shortcuts (x/k/arrows). Empty vocabulary is fine for these tests.
  await page.route("**/api/sample-tags", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
}

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() =>
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({
        state: { username: "alice", theme: "dark", tutorialSeen: true },
        version: 3,
      }),
    ),
  );
});

test("representative: picking a representative in the loupe PATCHes select", async ({ page }) => {
  // R2-M11 (#207): the per-thumb rep ⊙ button is gone — representative pick
  // lives in the loupe sidebar's "Set as representative" affordance (and the
  // `R` keystroke). Cover the loupe path here so the contact-sheet spec keeps
  // exercising the underlying selectExposure mutator.
  let selected = false;
  await mockCorpus(page);
  await page.route("**/api/exposures/2/select", async (route) => {
    selected = true;
    await route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({ id: 2, selected: true }),
    });
  });

  await page.goto("/sample/10/loupe");
  await expect(page.getByTestId("loupe-page")).toBeVisible();
  // Open the loupe on exposure 2 by clicking its strip thumbnail (frameNo 2 →
  // nth(1)). Every greenfield strip thumb shares data-testid="thumbnail".
  await page.getByTestId("thumbnail").nth(1).click();
  await page.getByRole("button", { name: "Set as representative" }).click();

  await expect.poll(() => selected).toBe(true);
});

// LO-REPDROP: the backend's Index-stage resolution never consults exposure
// status, so a DROPPED representative still carries forward. Dropping the
// representative in the loupe must surface the honest warning in the
// RepresentativeBox; restoring it clears the warning.
test("loupe: dropping the representative shows the rep-dropped warning; restore clears it", async ({ page }) => {
  await mockCorpus(page);
  // The exposures cache is patched OPTIMISTICALLY by the queue mutator's
  // onMutate (lib/queue/mutators/trivial.ts) — the fulfilled 200 here only
  // prevents the rollback, for drop (rejected) and restore (null) alike.
  await page.route("**/api/exposures/1/status", async (route) => {
    const body = route.request().postDataJSON() as { status: string | null };
    await route.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify({ id: 1, status: body.status }),
    });
  });

  await page.goto("/sample/10/loupe");
  await expect(page.getByTestId("loupe-page")).toBeVisible();
  // The loupe opens on the representative (exposure 1) — kept, so no warning.
  await expect(page.getByTestId("rep-dropped-warning")).toBeHidden();

  await page.keyboard.press("x"); // drop the representative
  await expect(page.getByTestId("rep-dropped-warning")).toBeVisible();
  await expect(page.getByTestId("rep-dropped-warning")).toContainText(
    "Warning. The representative frame is dropped.",
  );

  await page.keyboard.press("x"); // press Drop again to un-drop (toggle)
  await expect(page.getByTestId("rep-dropped-warning")).toBeHidden();
});

test("loupe-flip: arrow keys move between exposures in the loupe", async ({ page }) => {
  await mockCorpus(page);
  await page.goto("/sample/10/loupe");
  await expect(page.getByTestId("loupe-page")).toBeVisible();
  // Loupe opens on the representative (exposure 1, frame 1). The greenfield
  // page has no filename row; assert the active frame via the BigFrame caption.
  await expect(page.locator('[data-role="frame-caption"]')).toContainText("frame 1 of");
  await page.keyboard.press("ArrowRight");
  await expect(page.locator('[data-role="frame-caption"]')).toContainText("frame 2 of");
});

// Regression: a many-exposure loupe must NOT let the filmstrip balloon the grid's
// `1fr` column and shove the side panel off-screen. (Found in the Phase-4 loupe
// walkthrough: a CSS `1fr` track has an implicit `min-width:auto`, so the 20-thumb
// strip forced the column to ~1552px and pushed the side panel past the viewport.
// Fixed via `minmax(0,1fr)` + `min-w-0` on the column + `justify-center-safe` on
// the strip.) The mockup only ever drew a handful of thumbs, so it never exposed this.
test("loupe layout: a many-exposure filmstrip keeps the side panel on-screen", async ({ page }) => {
  await page.setViewportSize({ width: 1200, height: 800 });
  const MANY = Array.from({ length: 16 }, (_, i) => ({
    id: 500 + i, sample_id: 11, filename: `f${i + 1}.dat`, kind: "file",
    selected: i === 0, status: i === 0 ? "accepted" : null,
    image_path: null, image_version: "", tags: [],
    sources: [], trace_hash: null, analysis_inputs_hash: null,
  }));
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([EXPERIMENT]) }));
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify([
      { id: 11, experiment_id: 1, display_name: "D-many", name: "run11", notes: null, tags: [], q_units: "A-1" },
    ]) }));
  await page.route("**/api/samples/11/exposures*", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(MANY) }));
  await page.route("**/api/samples/11/messages", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));

  await page.goto("/sample/11/loupe");
  await expect(page.getByTestId("loupe-page")).toBeVisible();

  // The filmstrip genuinely overflows its column (so the guard is meaningful)…
  const strip = page.getByTestId("thumbnail-gallery");
  expect(await strip.evaluate((g) => g.scrollWidth > g.clientWidth + 1)).toBe(true);

  // …yet the side panel's right edge stays within the viewport.
  const box = await page.getByTestId("loupe-side-panel").boundingBox();
  expect(box).not.toBeNull();
  expect(box!.x + box!.width).toBeLessThanOrEqual(1200 + 1);
});

// SA-RESP (WCAG 1.4.10): the sheet grid's intrinsic min-width (~1054px with the
// checkbox track) exceeds a 1024px viewport. SheetTable wraps header + rows in
// ONE shared horizontal scroller with sticky identity columns, so a narrow
// viewport SCROLLS to the Phase column instead of silently clipping it.
//
// The mock sample carries MANY exposures on purpose: the table wrapper is
// min-w-min, so the scroller's width must stay at the grid's track-min sum
// (~1054px) while the ThumbnailGallery keeps scrolling INTERNALLY. (A
// max-content wrapper regressed this live: a real 8-exposure corpus unwrapped
// the flex-nowrap gallery and blew the scroller out to ~3400px — Status ended
// up 2.3 viewports away. Few-exposure mocks could not see that.)
test("narrow viewport: Phase column is reachable by scrolling; Sample column sticks", async ({ page }) => {
  await page.setViewportSize({ width: 1024, height: 800 });
  await mockCorpus(page);
  // Registered AFTER mockCorpus → takes precedence (Playwright matches the
  // most recently registered route first): sample 10 gets 12 exposures.
  const MANY = Array.from({ length: 12 }, (_, i) => ({
    id: 700 + i, sample_id: 10, filename: `m${i + 1}.dat`, kind: "file",
    selected: i === 0, status: "accepted",
    image_path: null, image_version: "", tags: [],
    sources: [], trace_hash: null, analysis_inputs_hash: null,
  }));
  await page.route("**/api/samples/10/exposures*", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(MANY) }));
  await page.goto("/experiments/1/corpus");
  await expect(page.getByTestId("sample-table-row").first()).toBeVisible();
  await expect(page.getByTestId("thumbnail")).toHaveCount(12);

  // The grid overflows the card at 1024 → the shared container scrolls…
  const scroll = page.getByTestId("sheet-scroll");
  const widths = await scroll.evaluate((el) => ({
    scrollWidth: el.scrollWidth, clientWidth: el.clientWidth,
  }));
  expect(widths.scrollWidth).toBeGreaterThan(widths.clientWidth);
  // …but only to the grid's track-min sum, NOT the gallery's unwrapped
  // max-content width (the min-w-min cap).
  expect(widths.scrollWidth).toBeLessThan(1200);

  // The many-exposure gallery still scrolls internally instead of widening
  // its fr track.
  const gallery = page.getByTestId("thumbnail-gallery");
  expect(
    await gallery.evaluate((g) => g.scrollWidth > g.clientWidth + 1),
  ).toBe(true);

  const status = page.getByRole("columnheader", { name: "Phase" });
  const sample = page.getByRole("columnheader", { name: "Sample" });
  const sampleBefore = await sample.boundingBox();
  expect(sampleBefore).not.toBeNull();

  // Scroll the container to the far right: the Phase column must be REACHABLE
  // (fully inside the viewport) — the 1.4.10 substance.
  await scroll.evaluate((el) => { el.scrollLeft = el.scrollWidth; });
  const statusBox = await status.boundingBox();
  expect(statusBox).not.toBeNull();
  expect(statusBox!.x).toBeGreaterThanOrEqual(0);
  expect(statusBox!.x + statusBox!.width).toBeLessThanOrEqual(1024 + 1);

  // The sticky Sample header did not move with the scroll.
  const sampleAfter = await sample.boundingBox();
  expect(sampleAfter).not.toBeNull();
  expect(Math.abs(sampleAfter!.x - sampleBefore!.x)).toBeLessThanOrEqual(1);
});

// F-STATE (DESIGN.md State taxonomy → Focus-visible): bespoke <button>/<a>
// elements inherit the keyboard-only 2px accent ring from ONE @layer base rule
// in styles.css (the loupe back button is the tracked LO-FOCUS case). Real
// browser only: :focus-visible heuristics (keyboard vs mouse) don't exist in
// JSDOM, and programmatic .focus() does not reliably trigger :focus-visible —
// so this Tab-walks to the button for the keyboard half and mouse-clicks a
// non-navigating bespoke button (a strip thumbnail) for the no-ring half.
test("focus-visible: keyboard focus draws the 2px accent ring on a bespoke button; mouse click draws none", async ({ page }) => {
  await mockCorpus(page);
  await page.goto("/sample/10/loupe");
  await expect(page.getByTestId("loupe-page")).toBeVisible();

  // Tab-walk (bounded) until the bespoke loupe-back button holds focus.
  let reached = false;
  for (let i = 0; i < 30; i++) {
    await page.keyboard.press("Tab");
    reached = await page.evaluate(
      () => document.activeElement?.getAttribute("data-testid") === "loupe-back",
    );
    if (reached) break;
  }
  expect(reached, "Tab order reaches the loupe back button").toBe(true);

  // The keyboard ring is the accent token, not the UA default. Resolve
  // var(--color-accent) through a probe element so both sides go through the
  // same computed-value serialization.
  const accent = await page.evaluate(() => {
    const probe = document.createElement("div");
    probe.style.color = "var(--color-accent)";
    document.body.appendChild(probe);
    const c = getComputedStyle(probe).color;
    probe.remove();
    return c;
  });
  const focused = await page.getByTestId("loupe-back").evaluate((el) => {
    const s = getComputedStyle(el);
    return { style: s.outlineStyle, width: s.outlineWidth, color: s.outlineColor };
  });
  expect(focused.style).toBe("solid");
  expect(focused.width).toBe("2px");
  expect(focused.color).toBe(accent);

  // Mouse clicks never show a ring. loupe-back navigates away, so use another
  // bespoke button that stays put: a filmstrip thumbnail (frame select only).
  const thumb = page.getByTestId("thumbnail").nth(1);
  await thumb.click();
  expect(await thumb.evaluate((el) => getComputedStyle(el).outlineStyle)).toBe("none");
});

test("wide viewport: the sheet does not scroll horizontally (unchanged layout)", async ({ page }) => {
  // Playwright's default 1280×720 viewport is the wide case.
  await mockCorpus(page);
  await page.goto("/experiments/1/corpus");
  await expect(page.getByTestId("sample-table-row").first()).toBeVisible();
  const widths = await page.getByTestId("sheet-scroll").evaluate((el) => ({
    scrollWidth: el.scrollWidth, clientWidth: el.clientWidth,
  }));
  expect(widths.scrollWidth).toBeLessThanOrEqual(widths.clientWidth + 1);
});

// LO-TAGTARGET (WCAG 2.2 §2.5.8 Target Size Minimum): the tag-remove × is a
// bare glyph on a dense pill, and the "+ tag" invite sits only 6px to its
// right (TagList gap-1.5) — too close for the spacing exemption, so the button
// itself must carry a ≥24×24 CSS-px border box. It also must not spill past
// the pill's right border, or it would overlap that neighboring target.
// jsdom has no layout, so this is pinned here in a real browser.
test("loupe tags: the tag-remove × has a ≥24×24 hit target inside the pill's right edge", async ({ page }) => {
  await mockCorpus(page);
  // Re-register /api/samples AFTER mockCorpus so this tagged payload wins
  // (Playwright matches the most recently registered route first).
  await page.route("**/api/samples", (r) =>
    r.fulfill({
      status: 200, contentType: "application/json",
      body: JSON.stringify([{
        ...SAMPLES[0],
        tags: [{ id: 1, key: "temperature", value: "37C", source: "user" }],
      }]),
    }));

  await page.goto("/sample/10/loupe");
  await expect(page.getByTestId("loupe-page")).toBeVisible();

  const pill = page.getByTestId("tag-pill").first();
  const remove = pill.getByRole("button", { name: "Remove" });
  await expect(remove).toBeVisible();

  const pillBox = await pill.boundingBox();
  const removeBox = await remove.boundingBox();
  expect(pillBox).not.toBeNull();
  expect(removeBox).not.toBeNull();
  expect(removeBox!.width).toBeGreaterThanOrEqual(24);
  expect(removeBox!.height).toBeGreaterThanOrEqual(24);
  // Right edge stays inside the pill's border box — the next interactive
  // target (the "+ tag" invite / the next pill's ×) is only 6px away.
  expect(removeBox!.x + removeBox!.width).toBeLessThanOrEqual(
    pillBox!.x + pillBox!.width,
  );
});
