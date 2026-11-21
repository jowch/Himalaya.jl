export const INGEST_STAGES = ["discovery", "analyzing"] as const;

export type IngestStage = (typeof INGEST_STAGES)[number];

const STAGE_LABELS: Record<IngestStage, string> = {
  discovery: "Reading files",
  analyzing: "Analyzing",
};

export function isIngestStage(v: string | undefined): v is IngestStage {
  return v !== undefined && (INGEST_STAGES as readonly string[]).includes(v);
}

export function stageLabel(stage: IngestStage): string {
  return STAGE_LABELS[stage];
}
