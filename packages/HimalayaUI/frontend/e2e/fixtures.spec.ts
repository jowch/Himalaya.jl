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

test("an error the page throws fails the test, even when every assertion passed", async ({ page }) => {
  test.fail();
  await page.goto("/favicon.svg");
  await page.evaluate(() => { setTimeout(() => { throw new Error("thrown after the last assertion"); }, 0); });
  await expect.poll(() => page.evaluate(() => new Promise((done) => setTimeout(() => done(true), 10)))).toBe(true);
});
