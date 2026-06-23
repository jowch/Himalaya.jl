import { useEffect, useState } from "react";
import type { JSX } from "react";
import { useParams, useNavigate } from "react-router-dom";
import { useQuery, useMutation } from "@tanstack/react-query";
import { ConfigurationBody } from "../components/ConfigurationBody";
import { DirectoryPickerField } from "../components/DirectoryPickerField";
import { useDraftExperiment } from "../../lib/draftExperiment";
import * as api from "../../api";
import { ProgressBar } from "../ui/ProgressBar";
import { Button } from "../ui/Button";
import { Card } from "../ui/Card";
import { Kicker } from "../ui/Kicker";
import { Input } from "../ui/Input";
import { IconButton } from "../ui/IconButton";
import { Chip } from "../ui/Chip";

/**
 * ConfigurationPage -- the Configuration tab shell (spec sec 8.1). Extracts
 * the experiment id from the route and delegates all body rendering to
 * ConfigurationBody (Task 19, Phase E2).
 *
 * First-run mode: absent `:id` (route `/experiments/new/config`) reads
 * path+patterns from `useDraftExperiment`, fetches the phase-1 manifest
 * (`GET /api/fs/manifest`), shows the indexing progress + file-pattern editor
 * + auto-derived geometry + latest files scanned. Gates Approve until the
 * manifest resolves. On Approve calls `createExperiment`, then `clear()`s the
 * draft and navigates to the new experiment's corpus.
 *
 * Later-edit mode: present `:id` (route `/experiments/:id/config`) keeps
 * the existing behaviour including the geometry/sources table.
 *
 * The two wrapper divs with data-testid="config-geometry-region" and
 * "config-sources-region" preserve the E1 test contract (ConfigurationPage
 * Phase E1 shell test) while the real content is now composed by
 * ConfigurationBody.
 */
export function ConfigurationPage(): JSX.Element {
  const { id } = useParams<{ id: string }>();

  if (!id) {
    return <ConfigurationFirstRun />;
  }

  const expId = Number(id);
  return (
    <div data-testid="config-geometry-region">
      {/* data-testid="config-sources-region" is also required by the E1 test.
          It lives inside the ConfigurationBody, but we need BOTH testids to be
          present. Wrap the body in a div that satisfies the geometry region test,
          then render ConfigurationBody which emits config-sources-region. */}
      <div data-testid="config-sources-region" className="contents">
        <ConfigurationBody experimentId={expId} />
      </div>
    </div>
  );
}

/** Default patterns used when the draft has no overrides -- mirrors the backend defaults. */
const DEFAULT_PATTERNS = {
  image: "{name}.tif",
  metadata: "{name}.prp",
  integration: "{name}.dat",
} as const;

/**
 * CoverageLine — truthful per-pattern match counts for the indexing preview.
 * The exposure count is the number of image matches; metadata and integration
 * are shown as `n / exposures`, flagged amber when a pattern covers fewer than
 * the exposures (or zero). A shortfall raises an explicit warning instead of
 * being hidden — integration (.dat) is matched against the analysis subtree, so
 * a 0 here means the pattern or the analysis directory is off, not a near-miss.
 */
function CoverageLine({ m }: { m: api.ManifestResponse }): JSX.Element {
  const exp = m.matched.image;
  const cov = [
    { key: "image", label: "Image", n: m.matched.image },
    { key: "metadata", label: "Metadata", n: m.matched.metadata },
    { key: "integration", label: "Integration", n: m.matched.integration },
  ] as const;
  const short = cov.filter((c) => c.key !== "image" && c.n < exp);
  return (
    <>
      <div className="flex flex-wrap items-center gap-x-3 gap-y-1" data-testid="manifest-coverage">
        {cov.map((c) => {
          const under = c.key !== "image" && c.n < exp;
          // Denominator hidden when there are no exposures (image === 0): a `/0`
          // reads like a real ratio when it is just "nothing to divide" (§6, D4).
          const showDenom = c.key !== "image" && exp > 0;
          return (
            <span key={c.key} className="text-caption">
              <span className="text-ink-soft">{c.label} </span>
              <span className={c.n === 0 || under ? "text-data text-warning" : "text-data text-ink"}>{c.n}</span>
              {showDenom && <span className="text-ink-faint">/{exp}</span>}
            </span>
          );
        })}
      </div>
      {short.length > 0 && exp > 0 && (
        <p className="text-caption text-warning" role="status" data-testid="manifest-warning">
          {short
            .map((c) =>
              c.key === "integration"
                ? `Integration matched ${c.n} of ${exp}. Check the integration pattern or the analysis directory.`
                : `Metadata matched ${c.n} of ${exp}. Check the metadata pattern.`,
            )
            .join(" ")}
        </p>
      )}
    </>
  );
}

/**
 * D4 zero-coverage block (§6). An exposure is the triple (image, metadata,
 * integration) and all three legs are essential "for now": a whole type matching
 * zero everywhere means the triple cannot be formed, so the experiment cannot be
 * ingested. Returns the first zero leg, or null. Single block tier (no advisory).
 */
function zeroCoverageType(
  m: api.ManifestResponse,
): "image" | "metadata" | "integration" | null {
  if (m.matched.image === 0) return "image";
  if (m.matched.metadata === 0) return "metadata";
  if (m.matched.integration === 0) return "integration";
  return null;
}

/**
 * EditableRow — a read-only value with an inline Edit/Done toggle (the same
 * click-to-edit idiom as the Analysis-directory field). Read mode shows
 * `display` (defaults to the raw value) + an "Edit" link; clicking it swaps in a
 * text Input that commits on Done / Enter / blur and reverts on Escape. The
 * committed string is trimmed. Used for the pattern fields (and, where editable,
 * geometry).
 */
function EditableRow({
  label,
  value,
  display,
  mono = false,
  inputTestId,
  valueTestId,
  onCommit,
}: {
  label: string;
  value: string;
  display?: JSX.Element | string;
  mono?: boolean;
  inputTestId?: string;
  valueTestId?: string;
  onCommit: (next: string) => void;
}): JSX.Element {
  const [editing, setEditing] = useState(false);
  const [draft, setDraft] = useState(value);
  // Re-seed the draft from upstream while not actively editing (e.g. a resolve
  // re-fills the value), so the field never shows a stale draft on next open.
  useEffect(() => { if (!editing) setDraft(value); }, [value, editing]);

  const commit = (): void => { onCommit(draft.trim()); setEditing(false); };
  const cancel = (): void => { setDraft(value); setEditing(false); };

  if (editing) {
    return (
      <div className="flex items-center justify-between gap-4">
        <span className="text-meta shrink-0">{label}</span>
        <span className="flex items-center gap-2">
          <Input
            value={draft}
            onValueChange={setDraft}
            inputSize="sm"
            mono={mono}
            autoFocus
            className="w-56"
            aria-label={label}
            {...(inputTestId ? { testId: inputTestId } : {})}
            onKeyDown={(e) => {
              if (e.key === "Enter") { e.preventDefault(); commit(); }
              else if (e.key === "Escape") { e.preventDefault(); cancel(); }
            }}
            onBlur={commit}
          />
          {/* preventDefault on mousedown keeps focus so the click commits before
              the input's blur fires a redundant first commit. */}
          <button
            type="button"
            className="text-caption text-accent shrink-0 hover:underline"
            onMouseDown={(e) => e.preventDefault()}
            onClick={commit}
          >
            Done
          </button>
        </span>
      </div>
    );
  }
  return (
    <div className="flex items-center justify-between gap-4">
      <span className="text-meta shrink-0">{label}</span>
      <span className="flex items-center gap-2 min-w-0">
        <span
          className={`text-data text-ink truncate text-right${mono ? " font-mono" : ""}`}
          {...(valueTestId ? { "data-testid": valueTestId } : {})}
        >
          {display ?? value}
        </span>
        <button
          type="button"
          className="text-caption text-accent shrink-0 hover:underline"
          onClick={() => setEditing(true)}
        >
          Edit
        </button>
      </span>
    </div>
  );
}

/** Maps a geometry source string to the short label shown in the chip. */
function sourceLabel(source: string | undefined | null): string {
  if (!source || source === "default") return "unset";
  if (source === "prp") return "prp";
  if (source === "setup") return "setup";
  if (source === "computed") return "computed";
  if (source === "user") return "edited";
  return "unset";
}

/**
 * GeometryEditRow — a single-value geometry field (flight path / pixel size /
 * energy) with the Edit/Done mechanism. Read mode shows the formatted value +
 * source chip + Edit; editing swaps in a numeric input (raw number, no unit)
 * committing on Done / Enter / blur, reverting on Escape. `onCommit` gets the
 * raw string; the parent parses + stores the override (blank reverts to derived).
 */
function GeometryEditRow({
  label, rawValue, displayValue, sourceLbl, onCommit,
}: {
  label: string;
  rawValue: number | null;
  displayValue: string;
  sourceLbl: string;
  onCommit: (raw: string) => void;
}): JSX.Element {
  const [editing, setEditing] = useState(false);
  const seed = rawValue == null ? "" : String(rawValue);
  const [draft, setDraft] = useState(seed);
  useEffect(() => { if (!editing) setDraft(rawValue == null ? "" : String(rawValue)); }, [rawValue, editing]);
  const commit = (): void => { onCommit(draft); setEditing(false); };
  const cancel = (): void => { setDraft(seed); setEditing(false); };
  return (
    <div className="flex items-center justify-between gap-4 py-1.5">
      <span className="text-meta text-ink-soft">{label}</span>
      {editing ? (
        <span className="flex items-center gap-2">
          <Input
            value={draft} onValueChange={setDraft} inputSize="sm" mono autoFocus
            className="w-28" aria-label={label}
            onKeyDown={(e) => {
              if (e.key === "Enter") { e.preventDefault(); commit(); }
              else if (e.key === "Escape") { e.preventDefault(); cancel(); }
            }}
            onBlur={commit}
          />
          <button type="button" className="text-caption text-accent shrink-0 hover:underline"
            onMouseDown={(e) => e.preventDefault()} onClick={commit}>Done</button>
        </span>
      ) : (
        <span className="flex items-center gap-2">
          <span className="text-data text-ink">{displayValue}</span>
          <Chip variant="static" size="sm">{sourceLbl}</Chip>
          <button type="button" className="text-caption text-accent shrink-0 hover:underline"
            onClick={() => setEditing(true)}>Edit</button>
        </span>
      )}
    </div>
  );
}

/** BeamCenterRow — the beam-center (x, y) pair with the Edit/Done mechanism: read
 *  shows "x, y px" + chip + Edit; editing shows two numeric inputs. `onCommit`
 *  gets both raw strings (committed together so neither half clobbers the other). */
function BeamCenterRow({
  x, y, sourceLbl, onCommit,
}: {
  x: number | null;
  y: number | null;
  sourceLbl: string;
  onCommit: (rx: string, ry: string) => void;
}): JSX.Element {
  const [editing, setEditing] = useState(false);
  const seedX = x == null ? "" : String(x);
  const seedY = y == null ? "" : String(y);
  const [dx, setDx] = useState(seedX);
  const [dy, setDy] = useState(seedY);
  useEffect(() => {
    if (!editing) { setDx(x == null ? "" : String(x)); setDy(y == null ? "" : String(y)); }
  }, [x, y, editing]);
  const commit = (): void => { onCommit(dx, dy); setEditing(false); };
  const cancel = (): void => { setDx(seedX); setDy(seedY); setEditing(false); };
  const keys = (e: React.KeyboardEvent): void => {
    if (e.key === "Enter") { e.preventDefault(); commit(); }
    else if (e.key === "Escape") { e.preventDefault(); cancel(); }
  };
  return (
    <div className="flex items-center justify-between gap-4 py-1.5">
      <span className="text-meta text-ink-soft">Beam center</span>
      {editing ? (
        /* Commit when focus leaves the WHOLE pair (relatedTarget outside this
           span) — not on the internal X→Y tab. A per-input onBlur on X would
           commit on that tab; with no onBlur on X, an outside-click from X would
           strand the edit. The Done button's onMouseDown-preventDefault keeps
           focus so its click commits before this blur fires. */
        <span
          className="flex items-center gap-1.5"
          onBlur={(e) => { if (!e.currentTarget.contains(e.relatedTarget as Node)) commit(); }}
        >
          <Input value={dx} onValueChange={setDx} inputSize="sm" mono autoFocus
            className="w-20" aria-label="Beam center X" onKeyDown={keys} />
          <span className="text-caption text-ink-faint">,</span>
          <Input value={dy} onValueChange={setDy} inputSize="sm" mono
            className="w-20" aria-label="Beam center Y" onKeyDown={keys} />
          <button type="button" className="text-caption text-accent shrink-0 hover:underline"
            onMouseDown={(e) => e.preventDefault()} onClick={commit}>Done</button>
        </span>
      ) : (
        <span className="flex items-center gap-2">
          <span className="text-data text-ink">
            {x != null && y != null ? `${x.toFixed(1)}, ${y.toFixed(1)} px` : "—"}
          </span>
          <Chip variant="static" size="sm">{sourceLbl}</Chip>
          <button type="button" className="text-caption text-accent shrink-0 hover:underline"
            onClick={() => setEditing(true)}>Edit</button>
        </span>
      )}
    </div>
  );
}

/**
 * First-run mode rendered at `/experiments/new/config`. Reads the picked ROOT
 * from the draft store, structurally RESOLVES it (`/api/fs/resolve`) into
 * auto-discovered defaults (name / data_dir / analysis_dir / setup_file) the
 * user CORRECTS in place, then fetches the phase-1 manifest against the resolved
 * data_dir + setup_file (geometry source). Gates Approve until indexing resolves.
 * On Approve creates the experiment with the CONFIRMED values + navigates to its corpus.
 *
 * Geometry is shown first-run (derived at index time from the resolved setup
 * file + sidecar .prp files -- no need to defer to after creation).
 */
function ConfigurationFirstRun(): JSX.Element {
  const navigate = useNavigate();
  const {
    root, resolved, name, data_dir, analysis_dir, setup_file, setup_ambiguous,
    patterns, geometry: geomOverride, applyResolved, patch, clear,
  } = useDraftExperiment();

  // 1. Resolve the picked root into correctable defaults (once per root).
  const resolveQuery = useQuery({
    queryKey: ["resolve", root],
    queryFn: () => api.resolveLayout(root),
    enabled: root !== "" && !resolved,
  });
  useEffect(() => {
    if (resolveQuery.data && !resolved) applyResolved(resolveQuery.data);
  }, [resolveQuery.data, resolved, applyResolved]);

  // 2. Manifest (geometry/files) against the resolved + corrected data_dir +
  //    setup_file. Re-keys when the user edits data_dir / patterns / setup_file.
  const manifestQuery = useQuery({
    queryKey: ["manifest", data_dir, analysis_dir, patterns, setup_file],
    queryFn: () => api.fetchManifest(data_dir, patterns, setup_file || undefined, analysis_dir || undefined),
    enabled: resolved && data_dir !== "",
  });

  const createMutation = useMutation({
    mutationFn: (body: api.CreateExperimentBody) => api.createExperiment(body),
    onSuccess: (created) => {
      clear();
      // Land on the combined scan + grouping-review surface (p1-grouping): loads
      // unfold live, and Confirm-groups gates on the scan finishing. It persists
      // through completion (the standalone route doesn't switch away like the
      // corpus state machine would).
      navigate(`/experiments/${created.id}/grouping`);
    },
  });

  const handleApprove = (): void => {
    // Send the CONFIRMED values (the create route uses them verbatim — no guessing).
    // Any geometry the user edited is sent as an override (persisted source='user'
    // at create; the scan derives the rest and never clobbers these).
    createMutation.mutate({
      name,
      data_dir,
      path: data_dir,
      ...(analysis_dir ? { analysis_dir } : {}),
      patterns,
      ...(Object.keys(geomOverride).length > 0 ? { geometry: geomOverride } : {}),
    });
  };

  // Set or clear one geometry override (a blank/NaN value reverts to auto-derived).
  const setGeomOverride = (key: keyof typeof geomOverride, raw: string): void => {
    const n = parseFloat(raw);
    const next = { ...geomOverride };
    if (raw.trim() === "" || Number.isNaN(n)) delete next[key];
    else next[key] = n;
    patch({ geometry: next });
  };
  // Effective value = user override (if set) wins over the derived manifest value.
  const gEff = (key: keyof typeof geomOverride, derived: number | null | undefined): number | null =>
    geomOverride[key] ?? derived ?? null;
  // Effective source label: an override reads as "edited"; else the derived source.
  const gSrc = (key: keyof typeof geomOverride, derivedSource: string | undefined): string =>
    geomOverride[key] !== undefined ? "edited" : sourceLabel(derivedSource);
  // Beam center commits BOTH halves in one patch — two setGeomOverride calls would
  // each read the same stale closure and the second would clobber the first.
  const setBeamCenter = (rx: string, ry: string): void => {
    const next = { ...geomOverride };
    const nx = parseFloat(rx); const ny = parseFloat(ry);
    if (rx.trim() === "" || Number.isNaN(nx)) delete next.beam_center_x; else next.beam_center_x = nx;
    if (ry.trim() === "" || Number.isNaN(ny)) delete next.beam_center_y; else next.beam_center_y = ny;
    patch({ geometry: next });
  };

  const manifest = manifestQuery.data;
  const indexing = manifestQuery.isFetching;
  const resolving = resolveQuery.isFetching && !resolved;
  // D4 hard block: a whole essential type matched zero everywhere (§6). Disables
  // Approve and swaps the headline; the partial-shortfall path stays non-blocking.
  const blockType = manifest && !indexing ? zeroCoverageType(manifest) : null;

  // Local editable field state, synced from the draft once it resolves; commits
  // on blur/Enter via patch() (which re-keys the manifest query where relevant).
  const [nameLocal, setNameLocal] = useState(name);
  // The title shows as text at rest; the pencil opens a wide rename field.
  const [editingName, setEditingName] = useState(false);
  const commitName = (): void => {
    const trimmed = nameLocal.trim();
    setNameLocal(trimmed);          // keep the rest-state display in sync with what's saved
    patch({ name: trimmed });
    setEditingName(false);
  };
  const [analysisDirLocal, setAnalysisDirLocal] = useState(analysis_dir);
  const [setupFileLocal, setSetupFileLocal] = useState(setup_file);
  const [imagePattern, setImagePattern] = useState(patterns.image ?? DEFAULT_PATTERNS.image);
  const [metadataPattern, setMetadataPattern] = useState(patterns.metadata ?? DEFAULT_PATTERNS.metadata);
  const [integrationPattern, setIntegrationPattern] = useState(patterns.integration ?? DEFAULT_PATTERNS.integration);

  // Analysis directory is shown read-only ("what we found"), RELATIVE to the
  // experiment root; the user reveals an editable relative-suffix field only to
  // CORRECT it (note 10). analysisDirLocal stays absolute (the backend contract);
  // the field shows/edits the suffix under the root.
  const [editingAnalysis, setEditingAnalysis] = useState(false);

  // Directory autocomplete for the analysis-dir field (same source as the
  // primary picker: GET /api/fs/suggest). Debounced; only while editing.
  const [analysisSuggest, setAnalysisSuggest] = useState<ReadonlyArray<string>>([]);
  useEffect(() => {
    if (!editingAnalysis || analysisDirLocal.trim() === "") { setAnalysisSuggest([]); return; }
    let live = true;
    const t = setTimeout(() => {
      void api.suggestPaths(analysisDirLocal)
        .then((r) => { if (live) setAnalysisSuggest(r.suggestions); })
        .catch(() => {});
    }, 200);
    return () => { live = false; clearTimeout(t); };
  }, [editingAnalysis, analysisDirLocal]);

  // When resolve populates the draft, seed the editable fields from it.
  useEffect(() => {
    if (resolved) {
      setNameLocal(name);
      setAnalysisDirLocal(analysis_dir);
      setSetupFileLocal(setup_file);
      // Re-seed pattern fields from any patterns the resolver detected (e.g. the
      // SSRL tot_files convention); fall back to the defaults when absent.
      setImagePattern(patterns.image ?? DEFAULT_PATTERNS.image);
      setMetadataPattern(patterns.metadata ?? DEFAULT_PATTERNS.metadata);
      setIntegrationPattern(patterns.integration ?? DEFAULT_PATTERNS.integration);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [resolved]);

  const commitPattern = (key: "image" | "metadata" | "integration", value: string): void => {
    patch({ patterns: { ...patterns, [key]: value } });
  };

  // Root-relative path display: the root is absolute (shown under the title);
  // the data/analysis locations read RELATIVE to it (e.g. `data`, `analysis`).
  // An absolute value outside the root degrades gracefully to showing absolute.
  const rootTrim = root.replace(/\/+$/, "");
  const stripRoot = (abs: string): string => {
    if (!abs) return "";
    if (abs === rootTrim) return ".";
    return abs.startsWith(rootTrim + "/") ? abs.slice(rootTrim.length + 1) : abs;
  };

  const geo = manifest?.geometry;

  // Pre-resolve: the root is set but discovery hasn't returned yet.
  if (!resolved) {
    return (
      <div data-testid="config-resolving" className="flex flex-col gap-3 px-8 py-8">
        <Kicker tone="accent">New experiment · Review before scan</Kicker>
        <p className="text-sm text-ink-soft">
          {resolving || root !== "" ? "Locating experiment files…" : "No experiment selected."}
        </p>
        {(resolving || root !== "") && <ProgressBar value={0} total={1} label="Locating files" />}
      </div>
    );
  }

  return (
    <>
      <div className="flex flex-col gap-6 px-8 py-8 pb-28">
        {/* 1. Header — editable Name + root path (left) · task framing (right) */}
        <div data-testid="config-header" className="flex items-start justify-between gap-8">
          <div className="min-w-0">
            <Kicker tone="accent" className="mb-1">
              New experiment · Review before scan
            </Kicker>
            {/* The title reads as text; the pencil opens a wide rename field. The
                "rename" cue lives on the pencil (its label/tooltip), not inline. */}
            {editingName ? (
              <Input
                variant="title"
                testId="config-name"
                value={nameLocal}
                onValueChange={setNameLocal}
                className="w-full"
                autoFocus
                onBlur={commitName}
                onKeyDown={(e) => {
                  if (e.key === "Enter") { e.preventDefault(); commitName(); }
                  else if (e.key === "Escape") { e.preventDefault(); setNameLocal(name); setEditingName(false); }
                }}
                aria-label="Experiment name"
              />
            ) : (
              <div className="flex items-baseline gap-2">
                <h1 data-testid="config-name" className="text-display text-ink truncate">
                  {nameLocal || "Untitled"}
                </h1>
                <IconButton
                  label="Rename experiment"
                  tone="ghost"
                  className="shrink-0 self-center"
                  onClick={() => setEditingName(true)}
                >
                  ✎
                </IconButton>
              </div>
            )}
            <p className="text-caption font-mono text-ink-soft mt-1 truncate">{root}</p>
          </div>
          {/* Task framing — what this screen is for (note 6). */}
          <p className="text-body text-ink-soft max-w-[36ch] shrink-0">
            We scanned this folder and filled in what we found. Check the
            locations and geometry below, correct anything off, then Approve to
            create the experiment and run the full scan.
          </p>
        </div>

        {/* 2. Indexing progress + truthful per-pattern coverage */}
        <div className="flex flex-col gap-2">
          <div className="flex flex-wrap items-baseline gap-x-4 gap-y-1">
            <span className="text-meta text-ink-soft">
              {indexing ? "Indexing directory…" : "Indexed"}
            </span>
            {manifest && (
              blockType ? (
                <span
                  className="text-meta text-warning font-semibold"
                  role="status"
                  data-testid="manifest-block"
                >
                  No {blockType} matched
                </span>
              ) : (
                <span className="text-meta text-ink-soft">
                  <span className="text-data text-ink">{manifest.matched.image}</span>
                  {manifest.matched.image === 1 ? " exposure found" : " exposures found"}
                </span>
              )
            )}
          </div>
          {manifest && <CoverageLine m={manifest} />}
          {blockType && (
            <p className="text-caption text-warning">
              Every exposure needs an image, metadata, and integration file. Fix the{" "}
              {blockType} pattern or directory, then re-index, to enable Approve.
            </p>
          )}
          {indexing && <ProgressBar value={0} total={1} label="Indexing directory" />}
        </div>

        {/* 3. Two cards side by side: Sources (locations + patterns) and Geometry */}
        <div className="grid grid-cols-2 gap-5">
          {/* 3a. SOURCES card — resolved locations (read-only) + file patterns.
              Labels are sentence-case ink (.text-meta); hints are dim caption,
              so label and hint read as distinct roles (note 9). */}
          <Card padding="lg">
            <Kicker tone="soft" className="mb-1">Sources · locations and patterns</Kicker>
            <p className="text-caption text-ink-soft mb-4">
              Locations are relative to the experiment root above.
            </p>
            <div className="flex flex-col gap-3">
              {/* Data directory — read-only confirmation, relative to the root
                  (note 7). You pointed at the root; this is what we resolved. */}
              <div className="flex items-center justify-between gap-4">
                <span className="text-meta shrink-0">Data directory</span>
                <span data-testid="config-data-dir" className="text-data text-ink-soft truncate text-right">
                  {stripRoot(data_dir) || "—"}
                </span>
              </div>
              {/* Analysis directory — resolved, read-only, relative; reveal a
                  relative-suffix field to correct it (note 10). */}
              <div className="flex flex-col gap-1.5">
                <div className="flex items-center justify-between gap-4">
                  <span className="text-meta shrink-0">Analysis directory</span>
                  {!editingAnalysis && (
                    <span className="flex items-center gap-2 min-w-0">
                      <span data-testid="config-analysis-dir" className="text-data text-ink-soft truncate text-right">
                        {stripRoot(analysisDirLocal) || "not found"}
                      </span>
                      <button
                        type="button"
                        className="text-caption text-accent shrink-0 hover:underline"
                        onClick={() => setEditingAnalysis(true)}
                      >
                        Edit
                      </button>
                    </span>
                  )}
                </div>
                {editingAnalysis && (
                  <div className="flex items-center gap-2">
                    {/* Absolute path with Tab/↑↓/↵ autocomplete — same picker the
                        data directory uses (note 10). The read-only display above
                        stays relative to the root. */}
                    <DirectoryPickerField
                      value={analysisDirLocal}
                      onChange={(abs) => { setAnalysisDirLocal(abs); patch({ analysis_dir: abs }); }}
                      suggestions={analysisSuggest}
                      validation={null}
                      ariaLabel="Analysis directory"
                      className="flex-1"
                    />
                    <button
                      type="button"
                      className="text-caption text-accent shrink-0 hover:underline"
                      onClick={() => { patch({ analysis_dir: analysisDirLocal.trim() }); setEditingAnalysis(false); }}
                    >
                      Done
                    </button>
                  </div>
                )}
              </div>
              <EditableRow
                label="Exposure pattern"
                value={imagePattern}
                mono
                inputTestId="config-image-pattern"
                onCommit={(v) => { setImagePattern(v); commitPattern("image", v); }}
              />
              <EditableRow
                label="Metadata pattern"
                value={metadataPattern}
                mono
                inputTestId="config-metadata-pattern"
                onCommit={(v) => { setMetadataPattern(v); commitPattern("metadata", v); }}
              />
              <EditableRow
                label="Integration pattern"
                value={integrationPattern}
                mono
                inputTestId="config-integration-pattern"
                onCommit={(v) => { setIntegrationPattern(v); commitPattern("integration", v); }}
              />
            </div>
            <p className="text-caption text-ink-soft mt-4">
              Editing a pattern re-runs the index.
            </p>
          </Card>

          {/* 3b. GEOMETRY card */}
          <Card padding="lg">
            <Kicker tone="soft" className="mb-4">Geometry · auto-derived</Kicker>
            {geo ? (
              <div className="flex flex-col divide-y divide-hair">
                {(() => {
                  const bx = gEff("beam_center_x", geo.beam_center_x);
                  const by = gEff("beam_center_y", geo.beam_center_y);
                  const fp = gEff("flight_path_m", geo.flight_path_m);
                  const px = gEff("pixel_size_um", geo.pixel_size_um);
                  const en = gEff("energy_kev", geo.energy_kev);
                  return (
                    <>
                      <BeamCenterRow
                        x={bx} y={by}
                        sourceLbl={
                          geomOverride.beam_center_x !== undefined || geomOverride.beam_center_y !== undefined
                            ? "edited" : sourceLabel(geo.beam_center_x_source)
                        }
                        onCommit={setBeamCenter}
                      />
                      <GeometryEditRow
                        label="Flight path"
                        rawValue={fp}
                        displayValue={fp != null ? `${fp.toFixed(4)} m` : "—"}
                        sourceLbl={gSrc("flight_path_m", geo.flight_path_m_source)}
                        onCommit={(raw) => setGeomOverride("flight_path_m", raw)}
                      />
                      <GeometryEditRow
                        label="Pixel size"
                        rawValue={px}
                        displayValue={px != null ? `${px.toFixed(1)} µm` : "—"}
                        sourceLbl={gSrc("pixel_size_um", geo.pixel_size_um_source)}
                        onCommit={(raw) => setGeomOverride("pixel_size_um", raw)}
                      />
                      <GeometryEditRow
                        label="Energy"
                        rawValue={en}
                        displayValue={en != null ? `${en.toFixed(1)} keV` : "—"}
                        sourceLbl={gSrc("energy_kev", geo.energy_kev_source)}
                        onCommit={(raw) => setGeomOverride("energy_kev", raw)}
                      />
                    </>
                  );
                })()}
              </div>
            ) : (
              <p className="text-meta text-ink-soft">
                Geometry will be derived after indexing.
              </p>
            )}
            {manifest?.discrepancies && manifest.discrepancies.length > 0 && (
              <div className="mt-3 flex flex-col gap-1">
                {manifest.discrepancies.map((d) => (
                  <p key={d.field} className="text-meta text-ink-soft">
                    {d.field}: {d.message}
                  </p>
                ))}
              </div>
            )}
            {/* Setup file (geometry source): exposed for correction only when
                discovery was ambiguous (none / multiple). Otherwise a quiet
                caption naming which file fed the geometry. */}
            {setup_ambiguous ? (
              <div className="mt-3 flex flex-col gap-1.5">
                <span className="text-meta text-warning">
                  We could not find the beamline setup file on its own.
                </span>
                <span className="text-caption text-ink-soft">
                  Point us to the setup_info file so we can read the beam center and flight path.
                </span>
                <Input
                  value={setupFileLocal}
                  onValueChange={setSetupFileLocal}
                  inputSize="sm" mono testId="config-setup-file"
                  placeholder="…/analysis/setup_info_*.txt"
                  onBlur={() => patch({ setup_file: setupFileLocal.trim() })}
                  onKeyDown={(e) => { if (e.key === "Enter") patch({ setup_file: setupFileLocal.trim() }); }}
                  className="w-full mt-0.5"
                />
              </div>
            ) : setup_file ? (
              <p className="text-caption text-ink-soft mt-3" data-testid="config-setup-caption">
                Geometry read from {setup_file.split("/").pop()}
              </p>
            ) : null}
            <p className="text-caption text-ink-soft mt-4">
              Derived automatically: prp = file sidecar, setup = beamline log. You can override later.
            </p>
          </Card>
        </div>

        {/* 4. LATEST FILES SCANNED card -- full width */}
        {manifest && (
          <Card padding="lg">
            <div className="flex items-center justify-between mb-3">
              <Kicker tone="soft">Latest files scanned</Kicker>
              <span className="text-meta text-ink-soft">
                {manifest.matched.image} matched
              </span>
            </div>
            {manifest.matched_files && manifest.matched_files.length > 0 ? (
              <div className="grid grid-cols-2 gap-x-8 gap-y-1">
                {manifest.matched_files.map((f) => (
                  <span key={f} className="text-data text-ink-soft">
                    {f}
                  </span>
                ))}
              </div>
            ) : (
              <p className="text-meta text-ink-soft">No files matched yet.</p>
            )}
          </Card>
        )}
      </div>

      {/* 5. Sticky funnel footer -- mirrors NewExperimentPage pattern */}
      <footer
        data-testid="funnel-footer"
        className="fixed bottom-0 left-0 right-0 z-30 flex items-center justify-between gap-4 border-t border-hair bg-paper px-8 py-3"
      >
        <span className="text-meta text-ink-soft">
          Approve creates the experiment and starts the full scan. Enabled when indexing finishes.
        </span>
        <div className="flex items-center gap-3">
          <Button
            variant="ghost"
            onClick={() => {
              clear();
              navigate("/experiments/new");
            }}
          >
            Cancel
          </Button>
          <Button
            variant="accent"
            data-testid="approve-button"
            disabled={indexing || !manifest || createMutation.isPending || blockType != null}
            onClick={handleApprove}
          >
            Approve
          </Button>
        </div>
      </footer>
    </>
  );
}
