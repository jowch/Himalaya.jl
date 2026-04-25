import { useState } from "react";
import { useAppState } from "../state";
import { useExperiment, useSamples } from "../queries";

/**
 * TitleButton — page-scoped title element.
 *
 * Resting state:    plain large sans-serif text — sample name only.
 * Hovered:          subtle button-like surface (background + border) so the user
 *                   knows it's interactive; tooltip shows the experiment.
 * Empty state:      faint italic prompt; clicking opens the nav modal.
 */
export function TitleButton(): JSX.Element {
  const experimentId = useAppState((s) => s.activeExperimentId);
  const sampleId     = useAppState((s) => s.activeSampleId);
  const openModal    = useAppState((s) => s.openNavModal);

  const experimentQ = useExperiment(experimentId ?? 0);
  const samplesQ    = useSamples(experimentId ?? 0);

  const experimentName = experimentId !== undefined
    ? (experimentQ.data?.name ?? `Experiment ${experimentId}`)
    : undefined;
  const sampleObj = sampleId !== undefined
    ? samplesQ.data?.find((s) => s.id === sampleId)
    : undefined;
  const sampleName = sampleObj
    ? (sampleObj.name ?? sampleObj.label ?? `Sample ${sampleId}`)
    : undefined;

  const step: "experiment" | "sample" =
    experimentId === undefined ? "experiment" : "sample";

  const [hovered, setHovered] = useState(false);

  const onClick = (): void => openModal(step);

  // Decide which label is the visible "title."
  const titleText = sampleName ?? (experimentName ? "pick a sample" : "pick an experiment");
  const isEmpty = !sampleName;
  // Below-title context line — visible while hovered.
  const subtitle = experimentName ?? "";

  return (
    <button
      type="button"
      data-testid="title-button"
      onClick={onClick}
      onMouseEnter={() => setHovered(true)}
      onMouseLeave={() => setHovered(false)}
      onFocus={() => setHovered(true)}
      onBlur={() => setHovered(false)}
      aria-label={
        sampleName
          ? `${experimentName ?? ""} ${sampleName}. Change`.trim()
          : "Open experiment / sample picker"
      }
      className={
        "relative inline-flex flex-col items-center gap-0.5 px-5 py-1.5 rounded-md " +
        "transition-colors " +
        "focus-visible:outline focus-visible:outline-2 focus-visible:outline-accent " +
        (hovered
          ? "bg-bg-elevated/70 border border-border"
          : "border border-transparent")
      }
    >
      <span
        className={
          "font-sans font-semibold text-[19px] leading-tight " +
          (isEmpty ? "text-fg-muted italic" : "text-fg")
        }
      >
        {titleText}
      </span>
      {/* Detail row: experiment + shortcut hint, only while hovered or empty */}
      <span
        className={
          "flex items-center gap-2 text-[10.5px] text-fg-dim " +
          "transition-opacity " +
          (hovered || isEmpty ? "opacity-100" : "opacity-0")
        }
        aria-hidden={!(hovered || isEmpty)}
      >
        {subtitle && <span>{subtitle}</span>}
        {subtitle && <span className="text-fg-dim/50">·</span>}
        <kbd
          data-testid="title-button-kbd"
          className="text-[10px] text-fg-dim
                     border border-border rounded px-1 py-px"
        >/</kbd>
      </span>
    </button>
  );
}
