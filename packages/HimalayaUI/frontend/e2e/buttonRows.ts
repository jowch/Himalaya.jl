import { expect, type Page } from "./fixtures";

export interface MismatchedPair {
  left: { text: string; height: number };
  right: { text: string; height: number };
}

export async function mismatchedButtonRows(page: Page): Promise<MismatchedPair[]> {
  return page.evaluate(() => {
    const MAX_GAP_PX = 24;
    const SAME_LINE_PX = 6;
    const visible = (el: Element): boolean => {
      const r = el.getBoundingClientRect();
      const s = getComputedStyle(el);
      return r.width > 0 && r.height > 0 && s.visibility !== "hidden" && s.display !== "none";
    };
    const isBordered = (el: Element): boolean => {
      const s = getComputedStyle(el);
      return parseFloat(s.borderTopWidth) > 0 && s.borderTopStyle !== "none" && s.borderTopColor !== "rgba(0, 0, 0, 0)";
    };
    const describe = (el: Element): string =>
      (el.getAttribute("aria-label") ?? el.textContent ?? "").trim().replace(/\s+/g, " ").slice(0, 40);
    const quarterPx = (n: number): number => Math.round(n * 4) / 4;
    const boxes = Array.from(document.querySelectorAll("button"))
      .filter((b) => visible(b) && isBordered(b))
      .map((b) => ({ el: b, r: b.getBoundingClientRect() }))
      .sort((a, b) => a.r.left - b.r.left);
    const out: { left: { text: string; height: number }; right: { text: string; height: number } }[] = [];
    for (const [i, a] of boxes.entries()) {
      for (const b of boxes.slice(i + 1)) {
        const gap = b.r.left - a.r.right;
        const sameLine = Math.abs(a.r.top + a.r.height / 2 - (b.r.top + b.r.height / 2)) <= SAME_LINE_PX;
        if (gap < 0 || gap > MAX_GAP_PX || !sameLine) continue;
        const [ha, hb] = [quarterPx(a.r.height), quarterPx(b.r.height)];
        if (ha !== hb) out.push({ left: { text: describe(a.el), height: ha }, right: { text: describe(b.el), height: hb } });
      }
    }
    return out;
  });
}

export async function expectButtonRowsAligned(page: Page, where: string): Promise<void> {
  expect(await mismatchedButtonRows(page), `bordered buttons side by side share one height (${where})`).toEqual([]);
}
