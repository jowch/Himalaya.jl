export function sampleEyebrow(setup: string | undefined, index: number, count: number): string | undefined {
  if (index < 0) return setup;
  if (setup === undefined) return `Sample ${index + 1} of ${count}`;
  return `${setup} · sample ${index + 1} of ${count}`;
}
