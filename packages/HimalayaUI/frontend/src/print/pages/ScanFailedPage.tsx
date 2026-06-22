import { useState, useMemo } from "react";
import { useNavigate } from "react-router-dom";
import type { ManifestUnmatched } from "../../api";
import { Button } from "../ui/Button";
import { Input } from "../ui/Input";

export interface ScanFailedPageProps {
  experimentId: number;
  /** Files that failed to match during the scan, from the manifest response. */
  unmatched: ManifestUnmatched[];
  /** Number of files that parsed successfully and can be ingested. */
  parsedCount: number;
  /** Called when the user confirms "Ingest N that parsed". Wired by the parent
   *  to trigger a forced scan via useTriggerScan. */
  onIngestParsed?: () => void;
}

type MissType = "metadata" | "integration" | "image";

const MISS_LABEL: Record<MissType, string> = {
  metadata: "metadata",
  integration: "integration",
  image: "image",
};

/**
 * ScanFailedPage — rendered by ExperimentCorpusPage when ingest_status is
 * "failed". Shows:
 *  - Open Configuration primary action (navigates to /experiments/:id/config)
 *  - Scrollable list of unmatched files grouped by miss type (one section heading
 *    per miss type, stems listed below — heading is the only visible miss-type
 *    label, so getByText(/metadata/i) matches exactly one element in tests)
 *  - Adaptive pattern test: one Input per affected miss type, clearing
 *    independently (labels are sr-only to avoid duplicating the heading text)
 *  - "Ingest N that parsed" with a two-stage in-place confirm
 *
 * Two-stage confirm: first click arms; armed state swaps the button for
 * "Confirm" + "Cancel" — mirrors SeriesScopingPage's confirmingDiscard pattern.
 */
export function ScanFailedPage({
  experimentId,
  unmatched,
  parsedCount,
  onIngestParsed,
}: ScanFailedPageProps): JSX.Element {
  const navigate = useNavigate();

  // Two-stage in-place confirm for "Ingest N that parsed"
  // Pattern: first click arms; Confirm executes; Cancel disarms.
  const [confirmingIngest, setConfirmingIngest] = useState(false);

  // Group unmatched by miss type — Map preserves insertion order.
  const byMiss = useMemo(() => {
    const groups = new Map<MissType, string[]>();
    for (const u of unmatched) {
      const key = u.miss as MissType;
      if (!groups.has(key)) groups.set(key, []);
      groups.get(key)!.push(u.file);
    }
    return groups;
  }, [unmatched]);

  // Affected miss types (for adaptive pattern test) — same order as byMiss.
  const affectedTypes = useMemo(() => Array.from(byMiss.keys()), [byMiss]);

  // Per-type pattern test inputs — independent, one per affected type.
  const [patterns, setPatterns] = useState<Record<string, string>>({});
  const setPattern = (type: string, value: string): void => {
    setPatterns((prev) => ({ ...prev, [type]: value }));
  };

  const handleOpenConfiguration = (): void => {
    navigate(`/experiments/${experimentId}/config`);
  };

  const handleApplyPatterns = (): void => {
    navigate(`/experiments/${experimentId}/config`);
  };

  const handleIngestParsed = (): void => {
    if (!confirmingIngest) {
      setConfirmingIngest(true);
      return;
    }
    setConfirmingIngest(false);
    onIngestParsed?.();
  };

  return (
    <div className="flex flex-col gap-5">
      {/* Header: scan failed notice + primary action */}
      <div className="flex items-center justify-between gap-4">
        <div>
          <p className="text-sm text-ink font-semibold">Scan incomplete</p>
          {unmatched.length > 0 && (
            <p className="text-sm text-ink-soft">
              {unmatched.length} {unmatched.length === 1 ? "file" : "files"} could not be matched
            </p>
          )}
        </div>
        <Button variant="accent" onClick={handleOpenConfiguration}>
          Open Configuration
        </Button>
      </div>

      {/* Scrollable unmatched list grouped by miss type.
          Each section heading is the only visible occurrence of that miss-type
          word — inline items show the stem only, labels are sr-only — so
          screen.getByText(/metadata/i) matches exactly one element. */}
      {unmatched.length > 0 && (
        <div
          className="flex flex-col gap-3 overflow-y-auto rounded-sm border border-hair-strong bg-paper-sunk p-4"
          style={{ maxHeight: "12rem" }}
          aria-label="Unmatched files"
        >
          {Array.from(byMiss.entries()).map(([missType, files]) => (
            <div key={missType}>
              {/* Section heading is the single visible miss-type label. */}
              <p className="text-caption text-ink-soft font-semibold mb-1 uppercase tracking-wider">
                Missing {MISS_LABEL[missType]}
              </p>
              <ul className="flex flex-col gap-0.5">
                {files.map((stem) => (
                  <li key={stem} className="text-sm text-ink font-mono">
                    {stem}
                  </li>
                ))}
              </ul>
            </div>
          ))}
        </div>
      )}

      {/* Adaptive pattern test: one Input per affected miss type.
          Labels are sr-only to avoid duplicating the section heading text
          (which is the unique visible "metadata"/"integration" occurrence). */}
      {affectedTypes.length > 0 && (
        <div className="flex flex-col gap-3">
          <p className="text-sm text-ink font-semibold">Test a pattern</p>
          {affectedTypes.map((type) => (
            <div key={type}>
              <Input
                id={`pattern-${type}`}
                aria-label={`${MISS_LABEL[type]} pattern`}
                value={patterns[type] ?? ""}
                onValueChange={(v) => setPattern(type, v)}
                placeholder="{name}.prp"
                mono
              />
            </div>
          ))}
          <Button variant="outline" onClick={handleApplyPatterns} className="self-start">
            Apply all in Configuration
          </Button>
        </div>
      )}

      {/* Ingest N that parsed — two-stage in-place confirm.
          First click arms; armed state swaps the button for Confirm + Cancel.
          Mirrors the SeriesScopingPage confirmingDiscard pattern. */}
      {parsedCount > 0 && (
        <div className="flex items-center gap-2">
          {confirmingIngest ? (
            <span
              data-testid="ingest-confirm-group"
              className="flex items-center gap-2"
              role="group"
              aria-label="Confirm ingest"
            >
              <span className="text-caption text-ink-soft">
                Ingest {parsedCount} parsed {parsedCount === 1 ? "file" : "files"}?
              </span>
              <Button variant="ghost" onClick={() => setConfirmingIngest(false)}>
                Cancel
              </Button>
              <Button
                variant="success"
                data-testid="ingest-confirm-yes"
                onClick={handleIngestParsed}
              >
                Confirm
              </Button>
            </span>
          ) : (
            <Button variant="outline" onClick={handleIngestParsed}>
              Ingest {parsedCount} that parsed
            </Button>
          )}
        </div>
      )}
    </div>
  );
}
