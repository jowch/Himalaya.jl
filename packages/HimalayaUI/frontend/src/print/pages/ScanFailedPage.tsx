import { useState, useMemo, useEffect } from "react";
import { useNavigate } from "react-router-dom";
import type { ManifestUnmatched } from "../../api";
import * as api from "../../api";
import { Button } from "../ui/Button";
import { Input } from "../ui/Input";
import { Card } from "../ui/Card";
import { Kicker } from "../ui/Kicker";
import { KbKey } from "../ui/KbKey";

type Patterns = { image?: string; metadata?: string; integration?: string };

export interface ScanFailedPageProps {
  experimentId: number;
  /** Files that failed to match during the scan, from the manifest response. */
  unmatched: ManifestUnmatched[];
  /** Number of files that parsed successfully and can be ingested. */
  parsedCount: number;
  /** Called when the user confirms "Ingest N that parsed". Wired by the parent
   *  to trigger a forced scan via useTriggerScan. */
  onIngestParsed?: () => void;
  /** Data dir for the live trial-pattern preview. When absent (the prop-less
   *  unit render), the readout does not fetch. */
  dataDir?: string;
  /** Analysis subtree for the live trial-pattern preview (where .dat lives). */
  analysisDir?: string;
  /** The experiment's current stored patterns, used both as the trial baseline
   *  and to show the missing-file `expected <stem><suffix>` annotation. */
  patterns?: Patterns;
}

type MissType = "metadata" | "integration" | "image";

// Heading per miss type. NOTE: the image group is "Image mismatch" (a near-miss
// extension), not "Missing image".
const MISS_HEADING: Record<MissType, string> = {
  metadata: "Missing metadata",
  integration: "Missing integration",
  image: "Image mismatch",
};

/** Reduce a glob/template pattern to its trailing extension-ish suffix so the
 *  "expected <stem><suffix>" annotation reads `expected HA_5_044_S1991.dat`.
 *  `*.dat` → `.dat`, `{name}.prp` → `.prp`, `*_total.dat` → `_total.dat`. */
function patternSuffix(pattern: string | undefined): string {
  if (!pattern) return "";
  // Drop everything up to and including the last `*` or `{name}` placeholder.
  const star = pattern.lastIndexOf("*");
  const brace = pattern.lastIndexOf("}");
  const cut = Math.max(star, brace);
  return cut >= 0 ? pattern.slice(cut + 1) : pattern;
}

/**
 * ScanFailedPage — rendered by ExperimentCorpusPage when ingest_status is
 * "failed" (spec §5.5, p1-failed). Shows:
 *  - A warn "Scan incomplete" kicker + a count-derived summary sentence.
 *  - Two cards: a scrollable grouped list of what did not parse (each row with
 *    its nearest existing file), and an adaptive "Test patterns" card with one
 *    trial input + live ✓ N/N per affected miss type.
 *  - An action bar: a config-note, the two-stage in-place partial-ingest
 *    confirm, and an "Open configuration" accent primary with an ↵ chip.
 *
 * Two-stage confirm: first click arms; armed state swaps the button for a
 * Cancel + success Confirm pair — mirrors SeriesScopingPage's confirmingDiscard.
 */
export function ScanFailedPage({
  experimentId,
  unmatched,
  parsedCount,
  onIngestParsed,
  dataDir,
  analysisDir,
  patterns,
}: ScanFailedPageProps): JSX.Element {
  const navigate = useNavigate();

  // Two-stage in-place confirm for "Ingest N that parsed".
  const [confirmingIngest, setConfirmingIngest] = useState(false);

  // Group unmatched by miss type — Map preserves insertion order.
  const byMiss = useMemo(() => {
    const groups = new Map<MissType, ManifestUnmatched[]>();
    for (const u of unmatched) {
      const key = u.miss as MissType;
      if (!groups.has(key)) groups.set(key, []);
      groups.get(key)!.push(u);
    }
    return groups;
  }, [unmatched]);

  // Affected miss types (for adaptive pattern test) — same order as byMiss.
  const affectedTypes = useMemo(() => Array.from(byMiss.keys()), [byMiss]);

  const total = parsedCount + unmatched.length;
  const imagePattern = patterns?.image ?? "*.tif";

  const handleOpenConfiguration = (): void => {
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

  // Summary breakdown clause: "mostly missing integration traces, a few metadata
  // mismatches." — list only affected types, dominant (largest) first.
  const breakdown = useMemo(() => {
    const entries = affectedTypes
      .map((t) => ({ t, n: byMiss.get(t)!.length }))
      .sort((a, b) => b.n - a.n);
    if (entries.length === 0) return "";
    const phrase = (t: MissType): string => {
      if (t === "integration") return "missing integration traces";
      if (t === "metadata") return "metadata mismatches";
      return "image mismatches";
    };
    const parts = entries.map((e, i) =>
      i === 0 ? `mostly ${phrase(e.t)}` : `a few ${phrase(e.t)}`,
    );
    return parts.join(", ");
  }, [affectedTypes, byMiss]);

  return (
    <div className="flex flex-col gap-5">
      {/* Header: warn kicker + count-derived summary sentence. */}
      <div className="flex flex-col gap-2">
        <Kicker tone="warning" as="h2">
          Scan incomplete
        </Kicker>
        <p className="text-sm text-ink-soft">
          <span className="font-semibold text-ink">
            {parsedCount} of {total} parsed.
          </span>{" "}
          {unmatched.length} couldn't be paired
          {breakdown ? `, ${breakdown}` : ""}. Each row shows the nearest file
          that does exist, so you can see what to match.
        </p>
      </div>

      {/* Two cards side by side: what didn't parse + adaptive trial patterns. */}
      <div className="grid grid-cols-2 gap-5">
        {/* Card A — Didn't parse list, grouped by miss type, scrollable. */}
        <Card padding="lg" border="strong" className="flex flex-col gap-3">
          <p className="text-caption text-ink-soft font-semibold uppercase tracking-wider">
            Didn't parse, {unmatched.length}{" "}
            {unmatched.length === 1 ? "file" : "files"}
          </p>
          <div
            className="flex flex-col gap-3 overflow-y-auto"
            style={{ maxHeight: "14rem" }}
            aria-label="Unmatched files"
          >
            {Array.from(byMiss.entries()).map(([missType, items]) => (
              <div key={missType} className="flex flex-col gap-1">
                <p
                  data-testid={`miss-heading-${missType}`}
                  className="text-caption text-ink-soft font-semibold uppercase tracking-wider"
                >
                  {MISS_HEADING[missType]} · {items.length}
                </p>
                <ul className="flex flex-col gap-1.5">
                  {items.map((u) => (
                    <li key={u.file} className="flex flex-col">
                      <span className="text-sm text-ink font-mono">{u.file}</span>
                      {missType === "image" && u.near ? (
                        <span className="text-caption text-ink-soft font-mono">
                          nearest{" "}
                          <span className="text-warning">{u.near}</span> (vs{" "}
                          {imagePattern})
                        </span>
                      ) : (
                        <span className="text-caption text-ink-faint font-mono">
                          expected {u.file}
                          {patternSuffix(patterns?.[missType])}
                        </span>
                      )}
                    </li>
                  ))}
                </ul>
              </div>
            ))}
          </div>
        </Card>

        {/* Card B — adaptive Test patterns: one trial input + live ✓ per type. */}
        <Card padding="lg" border="strong" className="flex flex-col gap-3">
          <div className="flex flex-col gap-1">
            <p className="text-caption text-ink-soft font-semibold uppercase tracking-wider">
              Test patterns
            </p>
            <p className="text-caption text-ink-soft">
              One field per file type with misses. Edit until each clears.
            </p>
          </div>
          {affectedTypes.length > 0 ? (
            <div className="flex flex-col gap-3">
              {affectedTypes.map((type) => (
                <TrialPattern
                  key={type}
                  type={type}
                  basePattern={patterns?.[type] ?? ""}
                  patterns={patterns}
                  dataDir={dataDir}
                  analysisDir={analysisDir}
                />
              ))}
              <Button
                variant="accent"
                onClick={handleOpenConfiguration}
                className="w-full"
              >
                Apply all in Configuration →
              </Button>
            </div>
          ) : (
            <p className="text-caption text-ink-faint">No patterns to test.</p>
          )}
        </Card>
      </div>

      {/* Action bar: config-note · spacer · two-stage ingest · primary CTA. */}
      <div className="flex items-center gap-4 border-t border-hair pt-4">
        <p className="text-caption text-ink-soft">
          A failed scan is usually a configuration issue.
        </p>
        <div className="ml-auto flex items-center gap-2">
          {parsedCount > 0 && (
            <>
              {confirmingIngest ? (
                <span
                  data-testid="ingest-confirm-group"
                  className="flex items-center gap-2"
                  role="group"
                  aria-label="Confirm ingest"
                >
                  <span className="text-caption text-ink-soft">
                    Ingest {parsedCount}, skip {unmatched.length}?
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
            </>
          )}
          <Button variant="accent" onClick={handleOpenConfiguration}>
            Open configuration
            <KbKey className="ml-1.5">↵</KbKey>
          </Button>
        </div>
      </div>
    </div>
  );
}

interface TrialPatternProps {
  type: MissType;
  basePattern: string;
  patterns: Patterns | undefined;
  dataDir: string | undefined;
  analysisDir: string | undefined;
}

/** One adaptive trial-pattern row: a mono Input + a live `✓ N / N` readout.
 *  When `dataDir` is provided, typing debounces a fetchManifest call swapping in
 *  the trial pattern for this miss type and reports matched[type] / matched.image
 *  (exposures). When absent, the readout simply does not fetch (a static hint),
 *  and never renders an extra textbox. */
function TrialPattern({
  type,
  basePattern,
  patterns,
  dataDir,
  analysisDir,
}: TrialPatternProps): JSX.Element {
  const [value, setValue] = useState(basePattern);
  const [result, setResult] = useState<{ matched: number; total: number } | null>(
    null,
  );

  useEffect(() => {
    if (!dataDir || value.trim() === "") {   // "" (unresolved experiment) no-ops too
      setResult(null);
      return;
    }
    let cancelled = false;
    const handle = setTimeout(() => {
      const trial: Patterns = { ...patterns, [type]: value };
      void api
        .fetchManifest(dataDir, trial, undefined, analysisDir)
        .then((r) => {
          if (cancelled) return;
          setResult({ matched: r.matched[type], total: r.matched.image });
        })
        .catch(() => {
          if (!cancelled) setResult(null);
        });
    }, 300);
    return () => {
      cancelled = true;
      clearTimeout(handle);
    };
  }, [value, dataDir, analysisDir, type, patterns]);

  const cleared = result != null && result.matched === result.total;

  return (
    <div className="flex items-center gap-3">
      <span className="text-meta text-ink-soft w-24 shrink-0">{type}</span>
      <Input
        id={`pattern-${type}`}
        aria-label={`${type} pattern`}
        value={value}
        onValueChange={setValue}
        placeholder="{name}.prp"
        mono
        className="flex-1"
      />
      {result != null && (
        <span
          className={`text-data shrink-0 ${cleared ? "text-success" : "text-warning"}`}
          role="status"
        >
          {cleared ? "✓ " : ""}
          {result.matched} / {result.total}
        </span>
      )}
    </div>
  );
}
