/**
 * `relativeTime("2026-05-14T10:00:00Z")` → "2h ago".
 *
 * Granularity matches typical app UX: "just now" < 1min,
 * "Nm ago" < 1h, "Nh ago" < 24h, "Nd ago" < 30d, ISO date thereafter.
 */
export function relativeTime(iso: string | null, nowMs: number = Date.now()): string | null {
  if (iso === null) return null;
  const t = Date.parse(iso);
  if (Number.isNaN(t)) return null;
  const diffMs = nowMs - t;
  const s = Math.floor(diffMs / 1000);
  if (s < 60) return "just now";
  const m = Math.floor(s / 60);
  if (m < 60) return `${m}m ago`;
  const h = Math.floor(m / 60);
  if (h < 24) return `${h}h ago`;
  const d = Math.floor(h / 24);
  if (d <= 30) return `${d}d ago`;
  return new Date(t).toISOString().slice(0, 10);
}
