import { useEffect, useState, type JSX } from "react";
import { useExperiment, useLoads, useUpdateExperiment, useTriggerScan } from "../../queries";
import { useUndoStack } from "../../hooks/useUndoStack";
import { suppressGlobalKeys } from "../../lib/keys";
import { Card } from "../ui/Card";
import { GeometryLedger, type GeometryRow } from "./GeometryLedger";
import { AcquisitionChart, type AcqSession } from "../plot/AcquisitionChart";
import { SourcesCard, type SourceRow } from "./SourcesCard";
import type { ExperimentPatch } from "../../api";

/** Undo entry: records the field key, its previous raw value (+ source) so
 *  Revert/Undo can replay the inverse write, and its new value so Redo can
 *  replay the forward write. */
interface UndoEntry {
  key: string;
  prevValue: number | string;
  prevSource: string;
  newValue: number | string;
}

/** Pattern field keys that trigger a forced rescan on edit (the set of file
 *  glob patterns: changing them changes which files get discovered). */
const PATTERN_KEYS = new Set<string>(["image_pattern", "metadata_pattern", "integration_pattern"]);

/** Geometry key -> patch key mapping. */
const GEOM_PATCH_KEY: Record<string, keyof ExperimentPatch> = {
  energy_kev:     "energy_kev",
  flight_path_m:  "flight_path_m",
  beam_center_x:  "beam_center_x",
  beam_center_y:  "beam_center_y",
  pixel_size_um:  "pixel_size_um",
  q_units:        "q_units",
};

/** Return the raw (number | string) value of a geometry field from the
 *  experiment, for seeding the inline editor and comparing on commit. */
function rawGeoValue(
  exp: NonNullable<ReturnType<typeof useExperiment>["data"]>,
  key: string,
): number | string | undefined {
  switch (key) {
    case "energy_kev":    return exp.energy_kev    ?? undefined;
    case "flight_path_m": return exp.flight_path_m ?? undefined;
    case "beam_center_x": return exp.beam_center_x ?? undefined;
    case "beam_center_y": return exp.beam_center_y ?? undefined;
    case "pixel_size_um": return exp.pixel_size_um ?? undefined;
    case "q_units":       return exp.q_units       ?? undefined;
    default:              return undefined;
  }
}

/** Convert a raw geometry value to a draft string for the inline editor.
 *  Numeric fields: String(n) (plain number, no units).
 *  String fields (q_units): the string itself. */
function rawToDraft(raw: number | string): string {
  return String(raw);
}

/** Parse the user's draft string back to a raw value (number | string).
 *  Returns undefined when the string cannot be parsed as a number for
 *  numeric fields. */
function parseDraft(key: string, draft: string): number | string | undefined {
  if (key === "q_units") return draft.trim() || undefined;
  const numeric = parseFloat(draft);
  if (Number.isNaN(numeric)) return undefined;
  return numeric;
}

export interface ConfigurationBodyProps {
  experimentId: number;
}

/**
 * ConfigurationBody -- composes the internals of the Configuration tab.
 *
 * Renders, in order:
 *   1. Editable description textarea (E1-owned column; E2 renders it here).
 *   2. Two-column grid: GeometryLedger (left) + Acquisition card (right).
 *   3. SourcesCard (full width).
 *
 * Owns the geometry override/undo state via useUndoStack. Override opens an
 * inline editor (no PATCH yet). On commit, if the parsed new value differs
 * from the current raw value, pushes undo + calls updateMutate. If the value
 * is unchanged, the edit is discarded silently (no PATCH, no mis-stamp of
 * source='user'). Escape always discards.
 *
 * TODO(Phase-D/E1): source discrepancyCount from geometry_discrepancy when
 * the backend field lands (spec sec 9.6 detection lives in the scan path;
 * Phase A-C carries no such column).
 */
export function ConfigurationBody({ experimentId }: ConfigurationBodyProps): JSX.Element {
  const { data: exp } = useExperiment(experimentId);
  const { data: loads } = useLoads(experimentId);
  const { mutate: updateMutate } = useUpdateExperiment(experimentId);
  const { mutate: rescanMutate } = useTriggerScan(experimentId);

  const undoStack = useUndoStack<UndoEntry>();

  // Description editing state.
  const [editingDesc, setEditingDesc] = useState(false);
  const [descDraft, setDescDraft] = useState("");

  // Geometry inline-edit state.
  const [editingKey, setEditingKey] = useState<string | undefined>(undefined);
  const [editDraft, setEditDraft] = useState("");

  // --- Build geometry rows ---
  const geometryRows: GeometryRow[] = exp
    ? [
        {
          key: "energy_kev",
          label: "Beam energy",
          value: exp.energy_kev != null ? `${exp.energy_kev.toFixed(2)} keV` : "—",
          source: exp.energy_kev_source,
        },
        {
          key: "flight_path_m",
          label: "Flight path",
          value: exp.flight_path_m != null ? `${exp.flight_path_m.toFixed(4)} m` : "—",
          source: exp.flight_path_m_source,
        },
        {
          key: "beam_center_x",
          label: "Beam center X",
          value: exp.beam_center_x != null ? `${exp.beam_center_x.toFixed(1)} px` : "—",
          source: exp.beam_center_x_source,
        },
        {
          key: "beam_center_y",
          label: "Beam center Y",
          value: exp.beam_center_y != null ? `${exp.beam_center_y.toFixed(1)} px` : "—",
          source: exp.beam_center_y_source,
        },
        {
          key: "pixel_size_um",
          label: "Pixel pitch",
          value: exp.pixel_size_um != null ? `${exp.pixel_size_um} µm` : "—",
          source: exp.pixel_size_um_source,
        },
        {
          key: "q_units",
          label: "q units",
          value: exp.q_units ?? "—",
          source: exp.q_units_source,
        },
      ]
    : [];

  // --- Build acquisition sessions from loads ---
  const acquisitionSessions: AcqSession[] = [];
  if (loads && loads.length > 0) {
    // Group loads by the date portion of start_time.
    const byDate = new Map<string, number[]>();
    for (const load of loads) {
      const label = load.start_time
        ? load.start_time.slice(0, 10)
        : `Load ${load.load_index}`;
      const existing = byDate.get(label) ?? [];
      existing.push(load.frame_count);
      byDate.set(label, existing);
    }
    for (const [label, frameCounts] of byDate) {
      acquisitionSessions.push({ label, loadFrameCounts: frameCounts });
    }
  }

  // --- Build sources rows ---
  const sourcesRows: SourceRow[] = exp
    ? [
        // Read-only directory rows -- directories are fixed at create.
        { key: "data_dir",      label: "Data directory",     value: exp.data_dir,      editable: false },
        { key: "analysis_dir",  label: "Analysis directory",  value: exp.analysis_dir,  editable: false },
        // Editable pattern rows -- editing a pattern triggers a backend rescan.
        { key: "image_pattern",       label: "Image pattern",       value: exp.image_pattern       ?? "", editable: true },
        { key: "metadata_pattern",    label: "Metadata pattern",    value: exp.metadata_pattern    ?? "", editable: true },
        { key: "integration_pattern", label: "Integration pattern", value: exp.integration_pattern ?? "", editable: true },
      ]
    : [];

  // --- Geometry override handlers ---

  /** Activate inline editing for a row. Seeds the draft with the raw numeric
   *  (or string) value -- no units suffix -- so the user edits a clean number. */
  const handleOverride = (key: string) => {
    if (!exp) return;
    const raw = rawGeoValue(exp, key);
    if (raw === undefined) return;
    setEditingKey(key);
    setEditDraft(rawToDraft(raw));
  };

  /** Commit an inline edit. Parses the draft; if the value is unchanged (or
   *  unparseable), exits edit mode without a PATCH. Only when the value
   *  actually changes: push undo + PATCH (backend stamps *_source='user'). */
  const handleEditCommit = () => {
    if (!exp || !editingKey) {
      setEditingKey(undefined);
      return;
    }
    const newRaw = parseDraft(editingKey, editDraft);
    const prevRaw = rawGeoValue(exp, editingKey);
    const row = geometryRows.find((r) => r.key === editingKey);
    // Exit edit mode first to prevent double-commit on blur after Enter.
    setEditingKey(undefined);
    if (newRaw === undefined || newRaw === prevRaw || !row) return;
    const patchKey = GEOM_PATCH_KEY[editingKey];
    if (!patchKey) return;
    // Record undo entry BEFORE the PATCH so Undo can restore the old value and
    // Redo can replay the new one.
    if (prevRaw !== undefined) {
      undoStack.push({ key: editingKey, prevValue: prevRaw, prevSource: row.source, newValue: newRaw });
    }
    updateMutate({ [patchKey]: newRaw });
  };

  /** Cancel inline edit -- no PATCH. */
  const handleEditCancel = () => {
    setEditingKey(undefined);
  };

  const handleUndo = () => {
    const entry = undoStack.pop();
    if (!entry) return;
    const patchKey = GEOM_PATCH_KEY[entry.key];
    if (patchKey) {
      updateMutate({ [patchKey]: entry.prevValue });
    }
  };

  const handleRedo = () => {
    const entry = undoStack.popRedo();
    if (!entry) return;
    const patchKey = GEOM_PATCH_KEY[entry.key];
    if (patchKey) {
      updateMutate({ [patchKey]: entry.newValue });
    }
  };

  // ⌘Z / ⌘⇧Z geometry undo/redo. Guarded by suppressGlobalKeys so a chord
  // inside the description textarea or an inline geometry Input stays native
  // text-undo (the field owns the keyboard there).
  useEffect(() => {
    function onKey(e: KeyboardEvent): void {
      if (!(e.metaKey || e.ctrlKey) || e.key.toLowerCase() !== "z") return;
      if (suppressGlobalKeys(e)) return;
      e.preventDefault();
      if (e.shiftKey) handleRedo();
      else handleUndo();
    }
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
    // handleUndo/handleRedo close only over stable mutate/undoStack refs.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const handleSourceEdit = (key: string, value: string) => {
    const patchKey = key as keyof ExperimentPatch;
    if (PATTERN_KEYS.has(key)) {
      // Pattern field: PATCH then force rescan so the new glob is applied immediately.
      updateMutate({ [patchKey]: value }, { onSuccess: () => rescanMutate(true) });
    } else {
      updateMutate({ [patchKey]: value });
    }
  };

  const handleRescan = () => {
    rescanMutate();
  };

  // --- Description edit handlers ---
  const beginDescEdit = () => {
    setDescDraft(exp?.description ?? "");
    setEditingDesc(true);
  };

  const commitDesc = () => {
    const trimmed = descDraft.trim();
    updateMutate({ description: trimmed || null });
    setEditingDesc(false);
  };

  const cancelDescEdit = () => setEditingDesc(false);

  return (
    <div className="flex flex-col gap-6">
      {/* 1. Editable description */}
      <Card padding="md">
        {editingDesc ? (
          <textarea
            className="w-full resize-none font-sans text-sm text-ink"
            rows={3}
            autoFocus
            value={descDraft}
            onChange={(e) => setDescDraft(e.target.value)}
            onBlur={commitDesc}
            onKeyDown={(e) => {
              if (e.key === "Enter" && (e.metaKey || e.ctrlKey)) commitDesc();
              if (e.key === "Escape") cancelDescEdit();
            }}
            aria-label="Experiment description"
          />
        ) : (
          <button
            type="button"
            className="w-full text-left font-sans text-sm"
            onClick={beginDescEdit}
            aria-label="Edit experiment description"
          >
            {exp?.description
              ? <span className="text-ink">{exp.description}</span>
              : <span className="text-ink-faint">Add a description...</span>
            }
          </button>
        )}
      </Card>

      {/* 2. Two-column grid: GeometryLedger + Acquisition card */}
      <div className="grid grid-cols-2 gap-4">
        <GeometryLedger
          rows={geometryRows}
          onOverride={handleOverride}
          onRevert={(_key) => handleUndo()}
          onUndo={handleUndo}
          canUndo={undoStack.canUndo}
          onRedo={handleRedo}
          canRedo={undoStack.canRedo}
          discrepancyCount={0}
          {...(editingKey !== undefined ? { editingKey, editDraft } : {})}
          onEditDraftChange={setEditDraft}
          onEditCommit={handleEditCommit}
          onEditCancel={handleEditCancel}
          // TODO(Phase-D/E1): source from geometry_discrepancy when field lands
        />

        <Card>
          <div className="flex items-center border-b border-hair px-4 py-3.5">
            <h3 className="text-headline text-ink">Acquisition</h3>
          </div>
          <div className="px-4 py-3">
            {acquisitionSessions.length === 0 ? (
              <p data-testid="acquisition-empty" className="text-sm text-ink-soft">
                No acquisition timeline yet. Rescan this experiment to record session timing.
              </p>
            ) : (
              <div data-testid="acquisition-timeline">
                <AcquisitionChart sessions={acquisitionSessions} />
              </div>
            )}
          </div>
        </Card>
      </div>

      {/* 3. SourcesCard full width */}
      <SourcesCard
        rows={sourcesRows}
        onEdit={handleSourceEdit}
        onRescan={handleRescan}
      />
    </div>
  );
}
