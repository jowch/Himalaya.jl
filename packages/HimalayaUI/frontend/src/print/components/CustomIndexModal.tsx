import { ModalShell, SegmentedControl, Button } from "../ui";
import { CustomPreview } from "../comb";
import type { CombSeries } from "../comb";
import { ModalHead } from "./ModalHead";
import { ModalFieldRow } from "./ModalFieldRow";
import { LatticeParamControl } from "./LatticeParamControl";
import { FitMetadata } from "./FitMetadata";
import { ModalFooter } from "./ModalFooter";

export interface CustomIndexFit {
  landed: number;
  total: number;
  snapped?: boolean;
}

export interface CustomIndexModalProps {
  open: boolean;
  onClose: () => void;
  onCancel: () => void;
  onAdd: () => void;
  symmetries: readonly string[];
  symmetry: string;
  onSymmetryChange: (s: string) => void;
  paramName: string;
  paramValue: string;
  paramMin: number;
  paramMax: number;
  paramStep?: number;
  onParamChange: (v: string) => void;
  unit?: string;
  previewSeries: CombSeries;
  observed: number[];
  /** Click-to-snap: clicking an observed peak in the preview emits its q so the
   *  consumer can set the lattice for the first reflection to land on it. */
  onSelectObserved?: (q: number) => void;
  fit: CustomIndexFit;
  /** Disables the "Add to assignment" action — set when the lattice parameter is
   *  empty, non-finite, or out of the symmetry's range, so a bad value can never
   *  round-trip to a server 400 (mirrors the trace q-add field's `disabled`). */
  addDisabled?: boolean;
  /** Flags the lattice number field `aria-invalid` to explain a disabled Add. */
  paramInvalid?: boolean;
  className?: string;
}

export function CustomIndexModal(props: CustomIndexModalProps): JSX.Element | null {
  const {
    open, onClose, onCancel, onAdd,
    symmetries, symmetry, onSymmetryChange,
    paramName, paramValue, paramMin, paramMax, paramStep, onParamChange, unit = "Å",
    previewSeries, observed, onSelectObserved, fit, addDisabled, paramInvalid, className,
  } = props;

  return (
    <ModalShell
      open={open}
      onClose={onClose}
      size="md"
      testId="custom-index-modal"
      aria-labelledby="cix-title"
      {...(className ? { className } : {})}
    >
      <ModalHead kicker="Speculative" title="Custom index" titleId="cix-title" onClose={onClose} />

      <div className="px-5 pt-4 pb-2 flex flex-col gap-4">
        <ModalFieldRow label="Symmetry">
          <SegmentedControl
            stretch
            aria-label="Crystal symmetry"
            options={symmetries.map((s) => ({ value: s, label: s }))}
            value={symmetry}
            onChange={onSymmetryChange}
          />
        </ModalFieldRow>

        <ModalFieldRow label="Lattice" labelSuffix={paramName}>
          <LatticeParamControl
            value={paramValue}
            min={paramMin}
            max={paramMax}
            {...(paramStep !== undefined ? { step: paramStep } : {})}
            onValueChange={onParamChange}
            unit={unit}
            invalid={paramInvalid ?? false}
          />
        </ModalFieldRow>

        <div className="bg-paper-sunk border border-hair rounded-sm px-2 py-1.5">
          <CustomPreview
            series={previewSeries}
            observed={observed}
            {...(onSelectObserved ? { onSelectObserved } : {})}
            className="block w-full"
          />
        </div>

        <FitMetadata
          landed={fit.landed}
          total={fit.total}
          paramName={paramName}
          paramValue={paramValue}
          unit={unit}
          {...(fit.snapped !== undefined ? { snapped: fit.snapped } : {})}
        />
      </div>

      <ModalFooter
        note="Drag the lattice until the teeth land on your peaks."
        actions={
          <>
            <Button variant="outline" onClick={onCancel}>Cancel</Button>
            <Button variant="accent" onClick={onAdd} disabled={addDisabled ?? false}>Add to assignment</Button>
          </>
        }
      />
    </ModalShell>
  );
}
