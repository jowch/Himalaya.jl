import { test, expect, type Page } from "@playwright/test";

const EXP = {
  id: 1, name: "TestRun", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-01-01",
};
const SAMPLES = [
  { id: 1, experiment_id: 1, label: "D1", name: "UX1", notes: null,
    tags: [{ id: 1, key: "lipid", value: "DOPC", source: "manual" }] },
  { id: 2, experiment_id: 1, label: "D2", name: "UX2", notes: null, tags: [] },
];

async function mockApi(page: Page, users: { id: number; username: string }[] = []): Promise<void> {
  await page.route("**/api/users", async (route) => {
    const req = route.request();
    if (req.method() === "GET") {
      await route.fulfill({ status: 200, contentType: "application/json",
        body: JSON.stringify(users) });
    } else if (req.method() === "POST") {
      const data = JSON.parse(req.postData() ?? "{}");
      await route.fulfill({ status: 201, contentType: "application/json",
        body: JSON.stringify({ id: users.length + 1, username: data.username }) });
    } else {
      await route.continue();
    }
  });
  await page.route("**/api/experiments/1", async (route) =>
    route.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXP) }));
  await page.route("**/api/experiments/1/samples", async (route) =>
    route.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(SAMPLES) }));
}

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() => { localStorage.clear(); });
});

test("shell loads with navbar, layout placeholders, and sample list", async ({ page }) => {
  await mockApi(page);
  await page.goto("/");

  await expect(page.locator("[data-testid='nav-logo']")).toHaveText("Himalaya");
  // Plan 5 wired MillerPlot + PhasePanel into the right column; no panes
  // remain as placeholders.
  await expect(page.locator("[data-testid='placeholder']")).toHaveCount(0);
  await expect(page.locator("[data-sample-id]")).toHaveCount(2);
  await expect(page.locator('[data-sample-id="1"] [data-testid="sample-label"]')).toHaveText("D1");
});

test("first-visit prompts user modal and submits new user", async ({ page }) => {
  await mockApi(page);
  await page.goto("/");

  await expect(page.locator("[role='dialog']")).toBeVisible();
  await expect(page.locator("[role='dialog'] h2")).toHaveText("Who are you?");

  await page.locator("input[placeholder='Enter username']").fill("alice");
  await page.locator("[role='dialog'] button", { hasText: "Continue" }).click();

  await expect(page.locator("[role='dialog']")).not.toBeVisible();
  await expect(page.locator("[data-testid='nav-user']")).toHaveText("alice");
});

test("clicking a sample updates the active state and breadcrumb", async ({ page }) => {
  await mockApi(page, [{ id: 1, username: "alice" }]);
  await page.addInitScript(() => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({ state: { username: "alice" }, version: 1 }),
    );
  });
  await page.goto("/");

  await expect(page.locator("[role='dialog']")).not.toBeVisible();
  await page.locator('[data-sample-id="2"]').click();
  await expect(page.locator('[data-sample-id="2"]')).toHaveAttribute("data-active", "true");
  await expect(page.locator("[data-testid='nav-breadcrumb']")).toContainText("D2");
});

test("filter narrows the sample list", async ({ page }) => {
  await mockApi(page, [{ id: 1, username: "alice" }]);
  await page.addInitScript(() => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({ state: { username: "alice" }, version: 1 }),
    );
  });
  await page.goto("/");

  await page.locator("input[placeholder='Filter samples\u2026']").fill("DOPC");
  await expect(page.locator("[data-sample-id]")).toHaveCount(1);
  await expect(page.locator("[data-sample-id]")).toHaveAttribute("data-sample-id", "1");
});

test("clicking the trace posts a manual peak", async ({ page }) => {
  let posted: { q: number } | null = null;

  await page.route("**/api/users", (r) => r.fulfill({ status: 200, body: "[]" }));
  await page.route("**/api/experiments/1", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify({
      id: 1, name: "demo", path: "/x", data_dir: "/x/data",
      analysis_dir: "/x/analysis", manifest_path: null,
      created_at: "2026-04-22T00:00:00Z",
    }),
  }));
  await page.route("**/api/experiments/1/samples", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify([
      { id: 10, experiment_id: 1, label: "A1", name: "s1", notes: null, tags: [] },
    ]),
  }));
  await page.route("**/api/samples/10/exposures", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify([
      { id: 100, sample_id: 10, filename: "demo", kind: "file",
        selected: true, tags: [], sources: [] },
    ]),
  }));
  await page.route("**/api/exposures/100/trace", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify({
      q:     [0.05, 0.10, 0.15, 0.20, 0.25, 0.30],
      I:     [ 100,   80,   90,  500,   85,   70],
      sigma: [  5,    4,    5,   25,    4,    3],
    }),
  }));
  await page.route("**/api/exposures/100/peaks", async (r) => {
    if (r.request().method() === "POST") {
      const body = r.request().postDataJSON() as { q: number };
      posted = body;
      return r.fulfill({
        status: 201,
        body: JSON.stringify({
          id: 999, exposure_id: 100, q: body.q, intensity: null,
          prominence: null, sharpness: null, source: "manual",
          stale_indices: 1,
        }),
      });
    }
    return r.fulfill({ status: 200, body: "[]" });
  });
  await page.route("**/api/exposures/100/indices", (r) =>
    r.fulfill({ status: 200, body: "[]" }));

  await page.addInitScript(() => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({
        state: { username: "alice", activeSampleId: 10, activeExposureId: 100 },
        version: 1,
      }),
    );
  });

  await page.goto("/");
  const viewer = page.locator('[data-testid="trace-viewer"]');
  await viewer.waitFor();
  const box = await viewer.boundingBox();
  if (!box) throw new Error("trace viewer has no bounding box");

  // Click somewhere in the middle of the plot — exact coords don't matter,
  // the snap will land on one of the data point q values.
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2);

  await expect.poll(() => posted?.q, { timeout: 2000 }).toBeGreaterThan(0);
});

test("hovering an alternative and clicking + posts to groups", async ({ page }) => {
  let posted: { index_id: number } | null = null;

  await page.route("**/api/users", (r) => r.fulfill({ status: 200, body: "[]" }));
  await page.route("**/api/experiments/1", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify({
      id: 1, name: "demo", path: "/x", data_dir: "/x/data",
      analysis_dir: "/x/analysis", manifest_path: null,
      created_at: "2026-04-22T00:00:00Z",
    }),
  }));
  await page.route("**/api/experiments/1/samples", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify([
      { id: 10, experiment_id: 1, label: "A1", name: "s1", notes: null, tags: [] },
    ]),
  }));
  await page.route("**/api/samples/10/exposures", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify([
      { id: 100, sample_id: 10, filename: "demo", kind: "file",
        selected: true, tags: [], sources: [] },
    ]),
  }));
  await page.route("**/api/exposures/100/trace", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify({
      q: [0.05, 0.1, 0.15, 0.2], I: [100, 80, 90, 70], sigma: [5, 4, 4, 3],
    }),
  }));
  await page.route("**/api/exposures/100/peaks", (r) =>
    r.fulfill({ status: 200, body: "[]" }));
  await page.route("**/api/exposures/100/indices", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify([
      { id: 1, exposure_id: 100, phase: "Pn3m", basis: 0.5, score: 1.0,
        r_squared: 0.99, lattice_d: 12.5, status: "candidate",
        predicted_q: [0.7, 0.86], peaks: [] },
      { id: 2, exposure_id: 100, phase: "Im3m", basis: 0.3, score: 0.6,
        r_squared: 0.71, lattice_d: 9.1, status: "candidate",
        predicted_q: [0.42, 0.6], peaks: [] },
    ]),
  }));
  await page.route("**/api/exposures/100/groups", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify([
      { id: 1, exposure_id: 100, kind: "auto", active: true, members: [1] },
    ]),
  }));
  await page.route("**/api/groups/1/members", async (r) => {
    if (r.request().method() === "POST") {
      posted = r.request().postDataJSON() as { index_id: number };
      return r.fulfill({
        status: 200,
        body: JSON.stringify({
          id: 2, exposure_id: 100, kind: "custom", active: true, members: [1, posted.index_id],
        }),
      });
    }
    return r.fulfill({ status: 404, body: "nope" });
  });

  await page.addInitScript(() => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: { username: "alice", activeSampleId: 10, activeExposureId: 100 },
      version: 1,
    }));
  });

  await page.goto("/");
  const altRow = page.locator('[data-alternative-id="2"]');
  await altRow.waitFor();
  await altRow.hover();
  const addBtn = altRow.getByRole("button", { name: /add index 2/i });
  await addBtn.click();
  await expect.poll(() => posted?.index_id, { timeout: 2000 }).toBe(2);
});

test("switching to Tags tab and adding a tag posts to /api/samples/:id/tags", async ({ page }) => {
  let posted: { key: string; value: string } | null = null;

  await page.route("**/api/users", (r) => r.fulfill({ status: 200, body: "[]" }));
  await page.route("**/api/experiments/1", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify({
      id: 1, name: "demo", path: "/x", data_dir: "/x/data",
      analysis_dir: "/x/analysis", manifest_path: null,
      created_at: "2026-04-22T00:00:00Z",
    }),
  }));
  await page.route("**/api/experiments/1/samples", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify([
      { id: 10, experiment_id: 1, label: "A1", name: "s1", notes: null, tags: [] },
    ]),
  }));
  await page.route("**/api/samples/10/exposures", (r) => r.fulfill({
    status: 200, body: "[]",
  }));
  await page.route("**/api/samples/10/tags", async (r) => {
    if (r.request().method() === "POST") {
      posted = r.request().postDataJSON() as { key: string; value: string };
      return r.fulfill({
        status: 201,
        body: JSON.stringify({
          id: 99, sample_id: 10,
          key: posted.key, value: posted.value, source: "manual",
        }),
      });
    }
    return r.fulfill({ status: 404, body: "nope" });
  });

  await page.addInitScript(() => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: { username: "alice", activeSampleId: 10, activeExposureId: undefined },
      version: 1,
    }));
  });

  await page.goto("/");
  await expect(page.getByRole("tab", { name: /exposures/i }))
    .toHaveAttribute("aria-selected", "true");

  await page.getByRole("tab", { name: /tags/i }).click();
  await expect(page.getByRole("tab", { name: /tags/i }))
    .toHaveAttribute("aria-selected", "true");

  await page.getByPlaceholder(/key/i).fill("lipid");
  await page.getByPlaceholder(/value/i).fill("DOPC");
  await page.getByRole("button", { name: /add tag/i }).click();

  await expect.poll(() => posted?.key, { timeout: 2000 }).toBe("lipid");
  await expect.poll(() => posted?.value, { timeout: 2000 }).toBe("DOPC");
});
