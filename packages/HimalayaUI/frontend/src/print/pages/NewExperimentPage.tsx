import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useNavigate } from "react-router-dom";
import { useQuery, useQueryClient } from "@tanstack/react-query";
import * as api from "../../api";
import { queryKeys } from "../../queries";
import { PageFrame } from "../components/PageFrame";
import { DirectoryPickerField } from "../components/DirectoryPickerField";
import { Button } from "../ui/Button";
import { Checkbox } from "../ui/Checkbox";
import { Input } from "../ui/Input";
import { KbKey } from "../ui/KbKey";
import { modCombo } from "../../lib/modKey";
import { Kicker } from "../ui/Kicker";

type Patterns = api.SetupPreview["patterns"];
type GeometryKey = "flight_path_m" | "beam_center_x" | "beam_center_y" | "energy_kev" | "pixel_size_um";
type CreateGeometry = NonNullable<api.CreateExperimentBody["geometry"]>;

interface RowEdit {
  include: boolean;
  name: string;
  data_dir: string;
  analysis_dir: string;
  setup_file: string;
  patterns: Patterns;
  geometry: Partial<Record<GeometryKey, string>>;
}

const GEOMETRY_FIELDS: ReadonlyArray<{ key: GeometryKey; label: string }> = [
  { key: "flight_path_m", label: "Flight path (m)" },
  { key: "beam_center_x", label: "Beam center x (px)" },
  { key: "beam_center_y", label: "Beam center y (px)" },
  { key: "energy_kev", label: "Energy (keV)" },
  { key: "pixel_size_um", label: "Pixel size (µm)" },
];

const count = (n: number): string => n.toLocaleString("en-US");

function initialEdit(s: api.SetupPreview): RowEdit {
  return {
    include: s.duplicate_of === null && s.exposures > 0 && s.layout_detected,
    name: s.name,
    data_dir: s.data_dir,
    analysis_dir: s.analysis_dir ?? "",
    setup_file: s.setup_file ?? "",
    patterns: s.patterns,
    geometry: {},
  };
}

function layoutEdited(s: api.SetupPreview, e: RowEdit): boolean {
  return e.data_dir !== s.data_dir
    || e.analysis_dir !== (s.analysis_dir ?? "")
    || e.setup_file !== (s.setup_file ?? "")
    || e.patterns.image !== s.patterns.image
    || e.patterns.metadata !== s.patterns.metadata
    || e.patterns.integration !== s.patterns.integration;
}

interface Summary {
  exposures: number;
  without_metadata: number;
  without_trace: number;
  geometry: api.ManifestGeometry;
}

function summaryFromManifest(m: api.ManifestResponse, fallback: api.ManifestGeometry): Summary {
  return {
    exposures: m.matched.image,
    without_metadata: m.unmatched.filter((u) => u.miss === "metadata").length,
    without_trace: m.unmatched.filter((u) => u.miss === "integration").length,
    geometry: m.geometry ?? fallback,
  };
}

function parsedGeometry(e: RowEdit): Partial<Record<GeometryKey, number>> {
  const out: Partial<Record<GeometryKey, number>> = {};
  for (const { key } of GEOMETRY_FIELDS) {
    const raw = e.geometry[key];
    if (raw === undefined || raw.trim() === "") continue;
    const v = Number(raw);
    if (Number.isFinite(v)) out[key] = v;
  }
  return out;
}

function createGeometry(g: api.ManifestGeometry, edits: Partial<Record<GeometryKey, number>>): CreateGeometry | undefined {
  const out: CreateGeometry = {};
  for (const { key } of GEOMETRY_FIELDS) {
    const edited = edits[key];
    const value = edited ?? g[key];
    if (value == null) continue;
    out[key] = value;
    out[`${key}_source`] = edited !== undefined ? "user" : g[`${key}_source`];
  }
  return Object.keys(out).length > 0 ? out : undefined;
}

function summaryLine(s: api.SetupPreview, sum: Summary, setupFile: string): string {
  const parts = [`${count(sum.exposures)} ${sum.exposures === 1 ? "exposure" : "exposures"}`];
  const fp = sum.geometry.flight_path_m;
  if (fp != null) parts.push(`flight path ${fp.toFixed(2)} m`);
  if (!s.layout_detected) parts.push("integration layout not recognised");
  if (setupFile === "") parts.push("no setup_info, so geometry needs checking");
  if (sum.without_metadata > 0) parts.push(`${count(sum.without_metadata)} without a .prp`);
  if (sum.without_trace > 0) parts.push(`${count(sum.without_trace)} without a trace yet`);
  return parts.join(" · ");
}

function manifestQuery(e: RowEdit) {
  return {
    queryKey: queryKeys.manifest(e.data_dir, e.patterns, e.setup_file, e.analysis_dir),
    queryFn: () => api.fetchManifest(e.data_dir, e.patterns, e.setup_file || undefined, e.analysis_dir || undefined),
  };
}

function DetailsField({ label, value, onChange, placeholder, testId }: {
  label: string; value: string; onChange: (v: string) => void; placeholder?: string; testId: string;
}): JSX.Element {
  return (
    <label className="flex flex-col gap-1">
      <span className="text-caption text-ink-soft">{label}</span>
      <Input inputSize="sm" value={value} onValueChange={onChange} testId={testId}
        {...(placeholder !== undefined ? { placeholder } : {})} />
    </label>
  );
}

function SetupRow({ setup, edit, onEdit, error }: {
  setup: api.SetupPreview;
  edit: RowEdit;
  onEdit: (patch: Partial<RowEdit>) => void;
  error: string | undefined;
}): JSX.Element {
  const [open, setOpen] = useState(false);
  const settled = useDebounced(edit, 300);
  const edited = layoutEdited(setup, settled);
  const manifest = useQuery({
    ...manifestQuery(settled),
    enabled: edited && settled.data_dir !== "",
  });
  const sum: Summary = edited && manifest.data
    ? summaryFromManifest(manifest.data, setup.geometry)
    : setup;
  const g = sum.geometry;
  const shown = (key: GeometryKey): string => {
    const v = g[key];
    return v == null ? "not found" : String(Number(v.toFixed(4)));
  };

  return (
    <li data-testid="setup-row" data-root={setup.root} className="flex flex-col gap-2 border-b border-hair py-3">
      <div className="flex items-center gap-3">
        <Checkbox
          checked={edit.include}
          disabled={setup.duplicate_of !== null}
          onChange={(include) => onEdit({ include })}
          aria-label={`Add ${edit.name}`}
        />
        <Input value={edit.name} onValueChange={(name) => onEdit({ name })} testId="setup-name"
          aria-label="Experiment name" className="min-w-0 flex-1" />
        <Button variant="ghost" onClick={() => setOpen((o) => !o)} aria-expanded={open}>
          {open ? "Hide details" : "Details"}
        </Button>
      </div>
      <p data-testid="setup-summary" className="text-meta text-ink-soft">
        {setup.duplicate_of !== null
          ? `Already added as "${setup.duplicate_of.name}".`
          : summaryLine(setup, sum, edit.setup_file)}
      </p>
      {error !== undefined && <p role="alert" className="text-meta text-error">{error}</p>}
      {open && (
        <div data-testid="setup-details" className="grid grid-cols-1 gap-3 md:grid-cols-2">
          <DetailsField label="Data folder" value={edit.data_dir} testId="details-data-dir"
            onChange={(data_dir) => onEdit({ data_dir })} />
          <DetailsField label="Integration folder" value={edit.analysis_dir} testId="details-analysis-dir"
            onChange={(analysis_dir) => onEdit({ analysis_dir })} />
          <DetailsField label="Setup file" value={edit.setup_file} testId="details-setup-file"
            onChange={(setup_file) => onEdit({ setup_file })} />
          <DetailsField label="Image files" value={edit.patterns.image} testId="details-image-pattern"
            onChange={(image) => onEdit({ patterns: { ...edit.patterns, image } })} />
          <DetailsField label="Metadata files" value={edit.patterns.metadata} testId="details-metadata-pattern"
            onChange={(metadata) => onEdit({ patterns: { ...edit.patterns, metadata } })} />
          <DetailsField label="Integration files" value={edit.patterns.integration} testId="details-integration-pattern"
            onChange={(integration) => onEdit({ patterns: { ...edit.patterns, integration } })} />
          {GEOMETRY_FIELDS.map(({ key, label }) => (
            <DetailsField key={key} label={label} value={edit.geometry[key] ?? ""} testId={`details-${key}`}
              placeholder={shown(key)}
              onChange={(v) => onEdit({ geometry: { ...edit.geometry, [key]: v } })} />
          ))}
        </div>
      )}
    </li>
  );
}

function useDebounced<T>(value: T, ms: number): T {
  const [v, setV] = useState(value);
  useEffect(() => {
    const t = setTimeout(() => setV(value), ms);
    return () => clearTimeout(t);
  }, [value, ms]);
  return v;
}

export function NewExperimentPage(): JSX.Element {
  const navigate = useNavigate();
  const queryClient = useQueryClient();
  const [path, setPath] = useState("");
  const [suggestions, setSuggestions] = useState<string[]>([]);
  const [edits, setEdits] = useState<Record<string, RowEdit>>({});
  const [errors, setErrors] = useState<Record<string, string>>({});
  const [submitting, setSubmitting] = useState(false);
  const picked = useDebounced(path.trim(), 200);

  useEffect(() => {
    if (picked === "") { setSuggestions([]); return; }
    let live = true;
    void api.suggestPaths(picked).then((r) => { if (live) setSuggestions(r.suggestions); }).catch(() => {});
    return () => { live = false; };
  }, [picked]);

  const preview = useQuery({
    queryKey: queryKeys.folderPreview(picked),
    queryFn: () => api.previewFolder(picked),
    enabled: picked !== "",
    retry: false,
  });
  const setups = useMemo(() => preview.data?.setups ?? [], [preview.data]);

  useEffect(() => {
    setEdits((cur) => Object.fromEntries(setups.map((s) => {
      const kept = cur[s.root];
      if (kept === undefined) return [s.root, initialEdit(s)];
      return [s.root, s.duplicate_of !== null ? { ...kept, include: false } : kept];
    })));
    setErrors((cur) => Object.fromEntries(Object.entries(cur).filter(([root]) => setups.some((s) => s.root === root))));
  }, [setups]);

  const chosen = setups.filter((s) => edits[s.root]?.include === true);
  const current = path.trim() === picked && !preview.isFetching && preview.data !== undefined && setups.every((s) => edits[s.root] !== undefined);
  const canSubmit = current && chosen.length > 0 && !submitting && chosen.every((s) => edits[s.root]!.name.trim() !== "");
  const [submitWhenReady, setSubmitWhenReady] = useState<string | undefined>(undefined);
  const inFlight = useRef(false);

  const submit = useCallback(async (): Promise<void> => {
    if (!canSubmit || inFlight.current) return;
    inFlight.current = true;
    setSubmitting(true);
    const created: number[] = [];
    const failed: Record<string, string> = {};
    for (const s of chosen) {
      const e = edits[s.root]!;
      let manifest: api.ManifestResponse | undefined;
      if (layoutEdited(s, e)) {
        try {
          manifest = await queryClient.fetchQuery(manifestQuery(e));
        } catch (err) {
          failed[s.root] = `Could not read this setup's files: ${err instanceof Error ? err.message : String(err)}`;
          continue;
        }
      }
      const geometry = createGeometry(manifest?.geometry ?? s.geometry, parsedGeometry(e));
      try {
        const exp = await api.createExperiment({
          name: e.name.trim(),
          data_dir: e.data_dir,
          path: e.data_dir,
          ...(e.analysis_dir !== "" ? { analysis_dir: e.analysis_dir } : {}),
          patterns: e.patterns,
          ...(geometry ? { geometry } : {}),
        });
        created.push(exp.id);
        setEdits((cur) => ({ ...cur, [s.root]: { ...cur[s.root]!, include: false } }));
      } catch (err) {
        failed[s.root] = err instanceof Error ? err.message : "Could not add this setup.";
      }
    }
    setErrors(failed);
    setSubmitting(false);
    inFlight.current = false;
    void queryClient.invalidateQueries({ queryKey: queryKeys.experiments });
    void queryClient.invalidateQueries({ queryKey: queryKeys.folderPreviews });
    if (Object.keys(failed).length > 0) return;
    navigate(created.length === 1 ? `/experiments/${created[0]}` : "/experiments");
  }, [canSubmit, chosen, edits, navigate, queryClient]);

  useEffect(() => {
    if (submitWhenReady === undefined) return;
    if (submitWhenReady !== path.trim()) setSubmitWhenReady(undefined);
    else if (current) {
      setSubmitWhenReady(undefined);
      if (canSubmit) void submit();
    }
  }, [submitWhenReady, path, current, canSubmit, submit]);

  useEffect(() => {
    const onKey = (ev: KeyboardEvent): void => {
      if (ev.key !== "Enter" || !(ev.metaKey || ev.ctrlKey)) return;
      if (ev.isComposing || ev.repeat || ev.shiftKey || ev.altKey || !canSubmit) return;
      const target = ev.target instanceof Element ? ev.target : null;
      if (target?.closest('[role="dialog"], [data-keys-trap]')) return;
      if (document.querySelector('[role="dialog"][aria-modal="true"]')) return;
      ev.preventDefault();
      ev.stopPropagation();
      void submit();
    };
    window.addEventListener("keydown", onKey, true);
    return () => window.removeEventListener("keydown", onKey, true);
  }, [submit, canSubmit]);

  const status = picked === ""
    ? "Pick a beamtime folder to add every setup in it, or one setup folder."
    : preview.isError
      ? preview.error instanceof api.ApiError && preview.error.status === 400
        ? "That folder does not exist."
        : `Could not read that folder: ${preview.error instanceof Error ? preview.error.message : "unknown error"}`
      : preview.isFetching && setups.length === 0
        ? "Looking for setups…"
        : setups.length === 0
          ? "No images found here or in the folders directly inside it."
          : setups.length === 1
            ? "Found 1 setup."
            : `Found ${setups.length} setups.`;

  return (
    <>
      <PageFrame width="home" className="px-6 py-8 pb-28">
        <Kicker tone="accent">New experiment</Kicker>
        <h1 className="text-display text-ink">Add a beamtime</h1>
        <div className="mt-6">
          <DirectoryPickerField value={path} onChange={setPath} onSubmit={() => (canSubmit ? void submit() : setSubmitWhenReady(path.trim()))} suggestions={suggestions} autoFocus />
          <p data-testid="preview-status" role="status" className="text-meta text-ink-soft mt-2">{status}</p>
        </div>
        {setups.length > 0 && (
          <ul data-testid="setup-list" className="mt-4 flex flex-col">
            {setups.map((s) => edits[s.root] && (
              <SetupRow key={s.root} setup={s} edit={edits[s.root]!} error={errors[s.root]}
                onEdit={(patch) => setEdits((cur) => ({ ...cur, [s.root]: { ...cur[s.root]!, ...patch } }))} />
            ))}
          </ul>
        )}
      </PageFrame>
      <footer
        data-testid="funnel-footer"
        className="fixed bottom-0 left-0 right-0 z-30 flex items-center justify-end gap-3 border-t border-hair bg-paper px-8 py-3"
      >
        <Button variant="ghost" onClick={() => navigate("/experiments")}>Cancel</Button>
        <Button variant="accent" data-testid="create-submit" disabled={!canSubmit} onClick={() => { void submit(); }}>
          {chosen.length > 1 ? `Add ${chosen.length} experiments` : "Add experiment"}
          <KbKey variant="frost" className="ml-2">{modCombo("↵")}</KbKey>
        </Button>
      </footer>
    </>
  );
}
