const collator = new Intl.Collator("en", { sensitivity: "base" });
const NUMBER = /(\d+(?:\.\d+)?)/;

function chunks(s: string): string[] {
  return s.split(NUMBER);
}

function compareChunks(a: readonly string[], b: readonly string[]): number {
  const n = Math.min(a.length, b.length);
  for (let i = 0; i < n; i++) {
    const x = a[i]!;
    const y = b[i]!;
    const c = i % 2 === 1 ? Number(x) - Number(y) : collator.compare(x, y);
    if (c !== 0) return Math.sign(c);
  }
  return a.length - b.length;
}

export function naturalCompare(a: string, b: string): number {
  return compareChunks(chunks(a), chunks(b)) || (a < b ? -1 : a > b ? 1 : 0);
}

export function naturalSortBy<T>(items: readonly T[], keyOf: (item: T) => string): T[] {
  return [...items].sort((x, y) => naturalCompare(keyOf(x), keyOf(y)));
}
