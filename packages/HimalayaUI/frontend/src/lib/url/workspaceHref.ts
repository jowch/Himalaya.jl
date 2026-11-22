export function workspaceHref(experimentId: number, sampleId: number, exposureId?: number | null): string {
  const base = `/experiments/${experimentId}/workspace/${sampleId}`;
  return exposureId ? `${base}?exposure=${exposureId}` : base;
}
