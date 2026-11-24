const collator = new Intl.Collator("en", { numeric: true, sensitivity: "base" });

export function naturalCompare(a: string, b: string): number {
  return collator.compare(a, b) || (a < b ? -1 : a > b ? 1 : 0);
}

export function naturalSortBy<T>(items: readonly T[], keyOf: (item: T) => string): T[] {
  return [...items].sort((x, y) => naturalCompare(keyOf(x), keyOf(y)));
}
