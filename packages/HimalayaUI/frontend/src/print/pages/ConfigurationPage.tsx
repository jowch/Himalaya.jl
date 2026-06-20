import { useParams } from "react-router-dom";
import { useExperiment } from "../../queries";
import { Card } from "../ui/Card";
import { Kicker } from "../ui/Kicker";

/**
 * ConfigurationPage — the Configuration tab body shell (spec §8.1). Lays out
 * three regions: editable description, the Geometry ledger, and Sources. The
 * region INTERNALS are E2 components — GeometryLedger (per-field value +
 * provenance chip + Override + discrepancy banner), AcquisitionTimeline, and
 * SourcesCard. E1 lays out + slots them as labelled placeholders.
 */
export function ConfigurationPage(): JSX.Element {
  const { id } = useParams<{ id: string }>();
  const expId = id ? Number(id) : 0;
  const exp = useExperiment(expId);

  return (
    <div className="flex flex-col gap-6">
      <Card padding="md" data-testid="config-geometry-region">
        <Kicker>Geometry</Kicker>
        {/* E2 GeometryLedger: per-field value + prp/setup/user chip + Override
            + the multi-setup discrepancy banner. E1 shows the derived values. */}
        <dl className="mt-3 grid grid-cols-[auto_1fr] gap-x-6 gap-y-1 text-sm">
          <dt className="text-ink-soft">Energy</dt>
          <dd className="font-mono text-ink">{exp.data?.energy_kev ?? "—"} keV</dd>
          <dt className="text-ink-soft">Flight path</dt>
          <dd className="font-mono text-ink">{exp.data?.flight_path_m ?? "—"} m</dd>
          <dt className="text-ink-soft">Beam center</dt>
          <dd className="font-mono text-ink">
            {exp.data ? `${exp.data.beam_center_x ?? "—"}, ${exp.data.beam_center_y ?? "—"}` : "—"}
          </dd>
        </dl>
      </Card>

      <Card padding="md" data-testid="config-sources-region">
        <Kicker>Sources</Kicker>
        {/* E2 SourcesCard: editable pattern rows + read-only directory rows.
            DECISION (E1→E2 contract):
            - image_pattern / metadata_pattern / integration_pattern: EDITABLE
              edit-in-place rows wired to updateExperiment (rescan-on-change).
            - data_dir / analysis_dir: READ-ONLY display rows (directory is
              fixed at creation; no DirectoryPickerField / path-validation here).
            E1 shows the static read-only dirs as a placeholder. */}
        <dl className="mt-3 grid grid-cols-[auto_1fr] gap-x-6 gap-y-1 text-sm">
          <dt className="text-ink-soft">Data</dt>
          <dd className="font-mono text-ink truncate">{exp.data?.data_dir ?? "—"}</dd>
          <dt className="text-ink-soft">Analysis</dt>
          <dd className="font-mono text-ink truncate">{exp.data?.analysis_dir ?? "—"}</dd>
        </dl>
      </Card>
    </div>
  );
}
