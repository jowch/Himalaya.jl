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
function faithful(label: string, v: number): boolean {
	return Math.abs(Number(label) - v) <= 1e-9 * Math.max(Math.abs(v), 1e-300);
}
function exponential(v: number): string {
	for (let digits = 0; digits <= 10; digits++) {
		const label = v.toExponential(digits).replace("e+", "e");
		if (faithful(label, v)) return label;
	}
	return formatAxis(v);
}

function isExponential(v: number): boolean {
	const av = Math.abs(v);
	return av !== 0 && (av < 1e-3 || av >= 1e4);
}
export function formatTicks(values: readonly number[]): string[] {
	if (values.every((v) => Number.isInteger(v) && Math.abs(v) < 1e4)) return values.map((v) => v.toFixed(0));
	const base = values.map(formatAxis);
	if (new Set(base).size === base.length && base.every((label, i) => faithful(label, values[i]!))) return base;
	const fixed = values.filter((v) => !isExponential(v));
	let digits = 0;
	while (digits < 10 && !fixed.every((v) => faithful(v.toFixed(digits), v))) digits++;
	return values.map((v) => (isExponential(v) ? exponential(v) : v.toFixed(digits)));
}