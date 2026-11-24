import { expect, type Locator, type Page } from "./fixtures";

export async function expectDialogFitsAndScrolls(page: Page, dialog: Locator, lastItem: Locator): Promise<void> {
  const viewport = page.viewportSize()!;
  await expect(dialog).toBeVisible();
  const panel = (await dialog.boundingBox())!;
  expect(panel.y + panel.height, "the dialog ends inside the window").toBeLessThanOrEqual(viewport.height);
  expect(panel.x + panel.width, "the dialog fits the window's width").toBeLessThanOrEqual(viewport.width);
  const body = dialog.getByTestId("overlay-body");
  const overflows = await body.evaluate((el) => el.scrollHeight > el.clientHeight);
  await expect(dialog.getByTestId("overlay-more"), "a fade marks the hidden rows below").toHaveCount(overflows ? 1 : 0);
  const header = dialog.getByTestId("overlay-header");
  const headerBefore = (await header.count()) > 0 ? await header.boundingBox() : null;
  await body.hover();
  await page.mouse.wheel(0, 100_000);
  await expect.poll(async () => {
    const [b, r] = await Promise.all([body.boundingBox(), lastItem.boundingBox()]);
    return r!.y + r!.height <= b!.y + b!.height + 0.5 && r!.y >= b!.y - 0.5;
  }, { message: "scrolling to the end brings the last item fully into view" }).toBe(true);
  await expect(dialog.getByTestId("overlay-more"), "the fade is gone at the end").toHaveCount(0);
  if (headerBefore) expect(await header.boundingBox(), "the header stays put while the body scrolls").toEqual(headerBefore);
}
