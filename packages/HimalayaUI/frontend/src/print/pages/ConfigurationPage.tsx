import { useEffect, useState } from "react";
import type { JSX } from "react";
import { useParams, useNavigate } from "react-router-dom";
import { useQuery, useMutation } from "@tanstack/react-query";
import { ConfigurationBody } from "../components/ConfigurationBody";
import { useDraftExperiment } from "../../lib/draftExperiment";
import * as api from "../../api";
import { ProgressBar } from "../ui/ProgressBar";
import { Button } from "../ui/Button";
import { Card } from "../ui/Card";
import { Kicker } from "../ui/Kicker";
import { Input } from "../ui/Input";
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
          return (
            <span key={c.key} className="text-caption">
              <span className="text-ink-soft">{c.label} </span>
              <span className={under ? "text-data text-warning" : "text-data text-ink"}>{c.n}</span>
              {c.key !== "image" && <span className="text-ink-faint">/{exp}</span>}
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

/** Maps a geometry source string to the short label shown in the chip. */
function sourceLabel(source: string | undefined | null): string {
  if (!source || source === "default") return "unset";
  if (source === "prp") return "prp";
  if (source === "setup") return "setup";
  if (source === "computed") return "computed";
  return "unset";
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
    patterns, applyResolved, patch, clear,
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
    createMutation.mutate({
      name,
      data_dir,
      path: data_dir,
      ...(analysis_dir ? { analysis_dir } : {}),
      patterns,
    });
  };

  const manifest = manifestQuery.data;
  const indexing = manifestQuery.isFetching;
  const resolving = resolveQuery.isFetching && !resolved;

  // Local editable field state, synced from the draft once it resolves; commits
  // on blur/Enter via patch() (which re-keys the manifest query where relevant).
  const [nameLocal, setNameLocal] = useState(name);
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
  const rootBase = rootTrim.split("/").filter(Boolean).pop() ?? "";
  const stripRoot = (abs: string): string => {
    if (!abs) return "";
    if (abs === rootTrim) return ".";
    return abs.startsWith(rootTrim + "/") ? abs.slice(rootTrim.length + 1) : abs;
  };
  const joinRoot = (rel: string): string => {
    const v = rel.trim();
    if (v === "") return "";
    if (v.startsWith("/")) return v; // absolute override (e.g. a shared analysis dir outside root)
    return `${rootTrim}/${v.replace(/^\.\//, "")}`;
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
            {/* The title is editable in place; a quiet cue says so (note 6). */}
            <div className="flex items-baseline gap-2">
              <Input
                variant="title"
                testId="config-name"
                value={nameLocal}
                onValueChange={setNameLocal}
                onBlur={() => patch({ name: nameLocal.trim() })}
                onKeyDown={(e) => {
                  if (e.key === "Enter") { patch({ name: nameLocal.trim() }); (e.target as HTMLElement).blur(); }
                }}
                aria-label="Experiment name"
              />
              <span className="text-caption text-ink-faint shrink-0">click to rename</span>
            </div>
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
              <span className="text-meta text-ink-soft">
                <span className="text-data text-ink">{manifest.matched.image}</span>
                {manifest.matched.image === 1 ? " exposure found" : " exposures found"}
              </span>
            )}
          </div>
          {manifest && <CoverageLine m={manifest} />}
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
                    {/* The root is fixed (dim prefix); the user types the suffix.
                        A leading "/" still accepts an absolute path. */}
                    <Input
                      value={stripRoot(analysisDirLocal)}
                      onValueChange={(rel) => {
                        const abs = joinRoot(rel);
                        setAnalysisDirLocal(abs);
                        patch({ analysis_dir: abs });
                      }}
                      inputSize="sm" mono
                      aria-label="Analysis directory (relative to root)"
                      leading={<span className="text-data text-ink-faint shrink-0">…/{rootBase}/</span>}
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
              <div className="flex items-center justify-between gap-4">
                <span className="text-meta shrink-0">Exposure pattern</span>
                <Input
                  value={imagePattern}
                  onValueChange={setImagePattern}
                  inputSize="sm" mono
                  onBlur={() => commitPattern("image", imagePattern)}
                  onKeyDown={(e) => { if (e.key === "Enter") commitPattern("image", imagePattern); }}
                  className="w-56"
                />
              </div>
              <div className="flex items-center justify-between gap-4">
                <span className="text-meta shrink-0">Metadata pattern</span>
                <Input
                  value={metadataPattern}
                  onValueChange={setMetadataPattern}
                  inputSize="sm" mono
                  onBlur={() => commitPattern("metadata", metadataPattern)}
                  onKeyDown={(e) => { if (e.key === "Enter") commitPattern("metadata", metadataPattern); }}
                  className="w-56"
                />
              </div>
              <div className="flex items-center justify-between gap-4">
                <span className="text-meta shrink-0">Integration pattern</span>
                <Input
                  value={integrationPattern}
                  onValueChange={setIntegrationPattern}
                  inputSize="sm" mono
                  onBlur={() => commitPattern("integration", integrationPattern)}
                  onKeyDown={(e) => { if (e.key === "Enter") commitPattern("integration", integrationPattern); }}
                  className="w-56"
                />
              </div>
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
                <div className="flex items-center justify-between gap-4 py-1.5">
                  <span className="text-meta text-ink-soft">Beam center</span>
                  <span className="flex items-center gap-2">
                    <span className="text-data text-ink">
                      {geo.beam_center_x != null && geo.beam_center_y != null
                        ? `${geo.beam_center_x.toFixed(1)}, ${geo.beam_center_y.toFixed(1)} px`
                        : "—"}
                    </span>
                    <Chip variant="static" size="sm">
                      {sourceLabel(geo.beam_center_x_source)}
                    </Chip>
                  </span>
                </div>
                <div className="flex items-center justify-between gap-4 py-1.5">
                  <span className="text-meta text-ink-soft">Flight path</span>
                  <span className="flex items-center gap-2">
                    <span className="text-data text-ink">
                      {geo.flight_path_m != null ? `${geo.flight_path_m.toFixed(4)} m` : "—"}
                    </span>
                    <Chip variant="static" size="sm">
                      {sourceLabel(geo.flight_path_m_source)}
                    </Chip>
                  </span>
                </div>
                <div className="flex items-center justify-between gap-4 py-1.5">
                  <span className="text-meta text-ink-soft">Pixel size</span>
                  <span className="flex items-center gap-2">
                    <span className="text-data text-ink">
                      {geo.pixel_size_um != null ? `${geo.pixel_size_um.toFixed(1)} µm` : "—"}
                    </span>
                    <Chip variant="static" size="sm">
                      {sourceLabel(geo.pixel_size_um_source)}
                    </Chip>
                  </span>
                </div>
                <div className="flex items-center justify-between gap-4 py-1.5">
                  <span className="text-meta text-ink-soft">Energy</span>
                  <span className="flex items-center gap-2">
                    <span className="text-data text-ink">
                      {geo.energy_kev != null ? `${geo.energy_kev.toFixed(1)} keV` : "—"}
                    </span>
                    <Chip variant="static" size="sm">
                      {sourceLabel(geo.energy_kev_source)}
                    </Chip>
                  </span>
                </div>
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
            disabled={indexing || !manifest || createMutation.isPending}
            onClick={handleApprove}
          >
            Approve
          </Button>
        </div>
      </footer>
    </>
  );
}
