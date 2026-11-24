import { test as base, expect } from "@playwright/test";

export const test = base.extend<{ unmockedRequests: string[]; pageCrashes: string[] }>({
  pageCrashes: [
    async ({ page }, use) => {
      const crashes: string[] = [];
      page.on("pageerror", (error) => crashes.push(error.message));
      page.on("console", (message) => {
        if (message.type() === "error" && message.text().startsWith("ErrorBoundary caught")) crashes.push(message.text());
      });
      await use(crashes);
      expect(crashes, "the page threw, or React gave up on a render and showed \"Something went wrong\"").toEqual([]);
    },
    { auto: true },
  ],
  unmockedRequests: [
    async ({ page }, use) => {
      const unmocked: string[] = [];
      await page.route(/\/api\//, (route) => {
        const req = route.request();
        const url = new URL(req.url());
        unmocked.push(`${req.method()} ${url.pathname}${url.search}`);
        return route.fulfill({
          status: 501,
          contentType: "application/json",
          body: JSON.stringify({ error: "no mock for this request in the spec" }),
        });
      });
      await use(unmocked);
      expect(
        [...new Set(unmocked)],
        "every /api request the page makes must be answered by a mock in the spec, not by whatever backend the Vite proxy reaches",
      ).toEqual([]);
    },
    { auto: true },
  ],
});

export { expect };
export type { Locator, Page, Route } from "@playwright/test";
