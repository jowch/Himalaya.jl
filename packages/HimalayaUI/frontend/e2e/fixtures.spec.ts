import { test, expect } from "./fixtures";

test("a request the spec did not mock fails the test", async ({ page }) => {
  test.fail();
  await page.goto("/favicon.svg");
  await page.evaluate(() => fetch("/api/not-mocked-anywhere").then(() => undefined, () => undefined));
});

test("a fallback from a spec's handler still reaches the catch-all", async ({ page, unmockedRequests }) => {
  await page.route("**/api/users", (r) => r.fallback());
  await page.goto("/favicon.svg");
  const status = await page.evaluate(() => fetch("/api/users", { method: "DELETE" }).then((r) => r.status));
  expect(status).toBe(501);
  expect(unmockedRequests).toEqual(["DELETE /api/users"]);
  unmockedRequests.length = 0;
});

test("an error the page throws is caught, even after the spec's last assertion", async ({ page, pageCrashes }) => {
  await page.goto("/favicon.svg");
  await page.evaluate(() => { setTimeout(() => { throw new Error("thrown after the last assertion"); }, 0); });
  await expect.poll(() => [...pageCrashes]).toEqual(["thrown after the last assertion"]);
  pageCrashes.length = 0;
});

test("a render that falls into the app's error boundary is caught", async ({ page, pageCrashes }) => {
  const json = (body: unknown) => ({ status: 200, contentType: "application/json", body: JSON.stringify(body) });
  await page.route("**/api/users", (r) => r.fulfill(json([])));
  await page.route("**/api/events*", (r) => r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
  await page.route("**/api/label-colors", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments", (r) => r.fulfill(json([])));
  await page.route("**/api/experiments/1", (r) => r.fulfill(json({ id: 1, name: "x", path: "/lab/x", ingest_status: "complete", setups: [] })));
  await page.route("**/api/experiments/1/**", (r) => r.fulfill(json([])));
  await page.route((url) => url.pathname === "/api/samples", (r) => r.fulfill(json([])));
  await page.goto("/experiments/1/workspace");
  await expect(page.getByRole("alert").getByText("Something went wrong.")).toBeVisible();
  await expect.poll(() => pageCrashes.some((c) => c.startsWith("ErrorBoundary caught"))).toBe(true);
  pageCrashes.length = 0;
});
