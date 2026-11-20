/**
 * Plain-decimal axis label formatter for SAXS log scales. Avoids Observable
 * Plot's default SI-suffix formatter, which renders 0.04 as "40m" (milli) —
 * unidiomatic for SAXS q values and actively confusing on first read. Uses
 * scientific notation only at the extremes where decimals would be unreadable.
 *
 * Used by `TraceViewer` (Index page) and `MultiTracePlot` (Compare page) so
 * both views render q ticks with the same vocabulary.
 */
export function formatAxis(d: number): string {
	const ad = Math.abs(d);
	if (ad === 0) return "0";
	if (ad < 1e-3 || ad >= 1e4)
		return d
			.toExponential(0)
			.replace("e+", "e")
			.replace("e-0", "e-")
			.replace("e0", "e");
	if (ad < 0.01) return d.toFixed(3);
	if (ad < 1) return d.toFixed(2);
	if (ad < 100) return d.toFixed(ad < 10 ? 1 : 0);
	return d.toFixed(0);
}

export function formatTicks(values: readonly number[]): string[] {
	if (values.length > 1 && values.every((v) => Number.isInteger(v) && Math.abs(v) < 1e4)) return values.map((v) => v.toFixed(0));
	const base = values.map(formatAxis);
	if (new Set(base).size === base.length) return base;
	for (let digits = 1; digits <= 6; digits++) {
		const labels = values.map((v) => (Math.abs(v) < 1e-3 || Math.abs(v) >= 1e4 ? formatAxis(v) : v.toFixed(digits)));
		if (new Set(labels).size === labels.length) return labels;
	}
	return base;
}
