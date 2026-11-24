export function sampleEyebrow(setup: string | undefined, index: number, count: number): string | undefined {
  if (index < 0) return setup;
  const place = `sample ${index + 1} of ${count}`;
  return setup ? `${setup} · ${place}` : `S${place.slice(1)}`;
}
