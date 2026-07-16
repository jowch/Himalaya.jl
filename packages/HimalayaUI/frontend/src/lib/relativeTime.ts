/**
 * Coarse "N min ago" relative-time label for status lines (e.g. the corpus
 * "scanned 2 min ago" header, §5.6). Coarse on purpose — a scan-freshness cue,
 * not a precise duration. `now` is injectable so callers/tests don't mock Date.
 */
export function formatRelativeTime(iso: string, now: number = Date.now()): string {
  const then = new Date(iso).getTime();
  if (Number.isNaN(then)) return "";
  const sec = Math.max(0, Math.round((now - then) / 1000));
  if (sec < 45) return "just now";
  const min = Math.round(sec / 60);
  if (min < 60) return `${min} min ago`;
  const hr = Math.round(min / 60);
  if (hr < 24) return `${hr} hr ago`;
  const days = Math.round(hr / 24);
  return `${days} day${days === 1 ? "" : "s"} ago`;
}
