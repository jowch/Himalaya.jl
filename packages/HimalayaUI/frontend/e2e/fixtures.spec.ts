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
