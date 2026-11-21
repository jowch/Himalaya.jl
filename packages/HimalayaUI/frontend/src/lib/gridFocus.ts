export type GridKey = "ArrowUp" | "ArrowDown" | "ArrowLeft" | "ArrowRight";

type Rect = { left: number; top: number; width: number; height: number };

export function isGridKey(key: string): key is GridKey {
  return key === "ArrowUp" || key === "ArrowDown" || key === "ArrowLeft" || key === "ArrowRight";
}

export function gridStep(rects: readonly Rect[], current: number, key: GridKey): number {
  const last = rects.length - 1;
  if (last < 0) return -1;
  if (current < 0) return 0;
  const linear = (d: 1 | -1): number => Math.min(last, Math.max(0, current + d));
  if (key === "ArrowRight") return linear(1);
  if (key === "ArrowLeft") return linear(-1);
  const here = rects[current]!;
  const down = key === "ArrowDown";
  const noLayout = rects.every((r) => r.width === 0 && r.height === 0);
  if (noLayout) return linear(down ? 1 : -1);
  const centre = here.left + here.width / 2;
  let best = current;
  let bestRow = Infinity;
  let bestDx = Infinity;
  rects.forEach((r, i) => {
    const dy = down ? r.top - here.top : here.top - r.top;
    if (dy <= 0) return;
    const dx = Math.abs(r.left + r.width / 2 - centre);
    if (dy < bestRow || (dy === bestRow && dx < bestDx)) {
      best = i;
      bestRow = dy;
      bestDx = dx;
    }
  });
  return best;
}
