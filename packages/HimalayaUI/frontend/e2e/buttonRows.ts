import { expect, type Page } from "./fixtures";

export interface MismatchedRow {
  row: string;
  buttons: { text: string; height: number }[];
}

export async function mismatchedButtonRows(page: Page): Promise<MismatchedRow[]> {
  return page.evaluate(() => {
    const visible = (el: Element): boolean => {
      const r = el.getBoundingClientRect();
      const s = getComputedStyle(el);
      return r.width > 0 && r.height > 0 && s.visibility !== "hidden" && s.display !== "none";
    };
    const describe = (el: Element): string =>
      (el.getAttribute("data-testid") ?? el.getAttribute("aria-label") ?? el.textContent ?? "").trim().replace(/\s+/g, " ").slice(0, 40);
    const isBordered = (b: Element): boolean => {
      const s = getComputedStyle(b);
      return parseFloat(s.borderTopWidth) > 0 && s.borderTopStyle !== "none" && s.borderTopColor !== "rgba(0, 0, 0, 0)";
    };
    const rows = new Map<Element, Element[]>();
    for (const b of Array.from(document.querySelectorAll("button"))) {
      if (!visible(b) || !isBordered(b) || !b.parentElement) continue;
      const row = b.parentElement;
      if (!getComputedStyle(row).display.includes("flex")) continue;
      rows.set(row, [...(rows.get(row) ?? []), b]);
    }
    const out: { row: string; buttons: { text: string; height: number }[] }[] = [];
    for (const [row, buttons] of rows) {
      if (buttons.length < 2) continue;
      const sameLine = buttons.filter((b) => Math.abs(b.getBoundingClientRect().top - buttons[0]!.getBoundingClientRect().top) < 12);
      const heights = sameLine.map((b) => Math.round(b.getBoundingClientRect().height * 4) / 4);
      if (sameLine.length >= 2 && new Set(heights).size > 1) {
        out.push({ row: describe(row), buttons: sameLine.map((b, i) => ({ text: describe(b), height: heights[i]! })) });
      }
    }
    return out;
  });
}

export async function expectButtonRowsAligned(page: Page, where: string): Promise<void> {
  expect(await mismatchedButtonRows(page), `bordered buttons side by side share one height (${where})`).toEqual([]);
}
