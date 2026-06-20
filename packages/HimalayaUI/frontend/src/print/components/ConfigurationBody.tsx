import { useState, type JSX } from "react";
import { useExperiment, useLoads, useUpdateExperiment } from "../../queries";
import { useUndoStack } from "../../hooks/useUndoStack";
import { Card } from "../ui/Card";
import { GeometryLedger, type GeometryRow } from "./GeometryLedger";
import { AcquisitionTimeline, type AcqSession } from "./AcquisitionTimeline";
import { SourcesCard, type SourceRow } from "./SourcesCard";
import type { ExperimentPatch } from "../../api";

/** Undo entry: records the field key, its previous display value, and its
 *  previous source so Revert/Undo can replay the inverse write. */
interface UndoEntry {
  key: string;
  prevValue: number | string;
  prevSource: string;
}

/** Geometry key -> patch key mapping. */
const GEOM_PATCH_KEY: Record<string, keyof ExperimentPatch> = {
  energy_kev:     "energy_kev",
  flight_path_m:  "flight_path_m",
  beam_center_x:  "beam_center_x",
  beam_center_y:  "beam_center_y",
  pixel_size_um:  "pixel_size_um",
  q_units:        "q_units",
};

/** Parse a display string like "9.00 keV" or "1.81 m" back to a number (or
 *  string for q_units). Returns undefined when unparseable. */
function parseDisplayValue(key: string, display: string): number | string | undefined {
  if (key === "q_units") return display.trim();
  // Strip trailing unit suffix and parse the numeric part.
  const numeric = parseFloat(display);
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
 * Owns the geometry override/undo state via useUndoStack. On Override commit
 * it pushes {key, prevValue, prevSource} and calls updateMutate. On Revert
 * or Undo it pops the entry and re-calls updateMutate with the previous value.
 *
 * TODO(Phase-D/E1): source discrepancyCount from geometry_discrepancy when
 * the backend field lands (spec sec 9.6 detection lives in the scan path;
 * Phase A-C carries no such column).
 */
export function ConfigurationBody({ experimentId }: ConfigurationBodyProps): JSX.Element {
  const { data: exp } = useExperiment(experimentId);
  const { data: loads } = useLoads(experimentId);
  const { mutate: updateMutate } = useUpdateExperiment(experimentId);

  const undoStack = useUndoStack<UndoEntry>();

  // Description editing state.
  const [editingDesc, setEditingDesc] = useState(false);
  const [descDraft, setDescDraft] = useState("");

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
          value: exp.pixel_size_um != null ? `${exp.pixel_size_um} um` : "—",
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

  // --- Handlers ---
  const handleOverride = (key: string) => {
    if (!exp) return;
    const row = geometryRows.find((r) => r.key === key);
    if (!row) return;
    // For simplicity, prompt-less override: the GeometryLedger calls onOverride
    // when the user clicks Override. This body opens an inline editing flow
    // using the existing row value as the starting draft.
    // NOTE: the GeometryLedger does not yet render an inline input; the
    // Override button triggers this callback and a future E2 task will add the
    // inline input. For now we record the intent in the undo stack so the
    // callback contract is satisfied.
    //
    // Derive the raw previous value from the row value string.
    const prevRaw = parseDisplayValue(key, row.value);
    if (prevRaw === undefined) return;

    undoStack.push({ key, prevValue: prevRaw, prevSource: row.source });
    // Mutate using the SAME value (no-op geometry change) to stamp source=user
    // so the source chip updates correctly. A real override will pass the user's
    // edited value once inline-editing lands.
    const patchKey = GEOM_PATCH_KEY[key];
    if (patchKey) {
      updateMutate({ [patchKey]: prevRaw } as ExperimentPatch);
    }
  };

  const handleRevert = (key: string) => {
    const entry = undoStack.pop();
    if (!entry) return;
    const patchKey = GEOM_PATCH_KEY[entry.key];
    if (patchKey) {
      updateMutate({ [patchKey]: entry.prevValue } as ExperimentPatch);
    }
    void key; // key identifies the row; entry.key matches it
  };

  const handleUndo = () => {
    const entry = undoStack.pop();
    if (!entry) return;
    const patchKey = GEOM_PATCH_KEY[entry.key];
    if (patchKey) {
      updateMutate({ [patchKey]: entry.prevValue } as ExperimentPatch);
    }
  };

  const handleSourceEdit = (key: string, value: string) => {
    const patchKey = key as keyof ExperimentPatch;
    updateMutate({ [patchKey]: value } as ExperimentPatch);
  };

  const handleRescan = () => {
    // Rescan is triggered via the ExperimentShell toolbar (scan button). The
    // SourcesCard "Rescan now" button is a convenience shortcut -- it calls
    // the same triggerScan API via the shell. For now this is a no-op
    // placeholder; wiring to useTriggerScan is a future Phase C/E task.
  };

  // --- Description edit handlers ---
  const beginDescEdit = () => {
    setDescDraft(exp?.description ?? "");
    setEditingDesc(true);
  };

  const commitDesc = () => {
    const trimmed = descDraft.trim();
    updateMutate({ description: trimmed || null } as ExperimentPatch);
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
          onRevert={handleRevert}
          onUndo={handleUndo}
          canUndo={undoStack.canUndo}
          discrepancyCount={0}
          // TODO(Phase-D/E1): source from geometry_discrepancy when field lands
        />

        <Card>
          <div className="flex items-center border-b border-hair px-4 py-3.5">
            <h3 className="text-headline text-ink">Acquisition</h3>
          </div>
          <div className="px-4 py-3">
            <AcquisitionTimeline sessions={acquisitionSessions} />
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
