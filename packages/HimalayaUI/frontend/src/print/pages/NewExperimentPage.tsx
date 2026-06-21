import { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import { useExperiments } from "../../queries";
import * as api from "../../api";
import { useDraftExperiment } from "../../lib/draftExperiment";
import { PageFrame } from "../components/PageFrame";
import { DirectoryPickerField } from "../components/DirectoryPickerField";
import { Button } from "../ui/Button";
import { Kicker } from "../ui/Kicker";
import type { ValidatePathResponse } from "../../api";

/**
 * NewExperimentPage — /experiments/new (spec §8.7). Directory picker +
 * suggestions + validation; primary action commits the validated path to a
 * client-side draft and navigates to /experiments/new/config (first-run
 * Configuration). No DB row is created here — creation happens at Approve.
 *
 * T4.0: two-phase funnel handoff. "Review →" replaces the old "Scan and
 * create" action. The dup-dir guard (inline error + disabled submit) is
 * preserved from the pre-T4.0 implementation.
 */
export function NewExperimentPage(): JSX.Element {
  const navigate = useNavigate();
  const { setDraft, clear } = useDraftExperiment();

  const { data: experiments } = useExperiments();

  const [path, setPath] = useState("");
  const [suggestions, setSuggestions] = useState<string[]>([]);
  const [validation, setValidation] = useState<ValidatePathResponse | null>(null);
  const [submitting, setSubmitting] = useState(false);

  // One experiment per directory: warn (and block) if this dir is already taken.
  const trimmedPath = path.trim();
  const duplicateOf =
    trimmedPath === ""
      ? undefined
      : (experiments ?? []).find((e) => e.data_dir.trim() === trimmedPath);

  // Debounced suggestion + validation fetch on path change.
  useEffect(() => {
    if (path.trim() === "") {
      setSuggestions([]);
      setValidation(null);
      return;
    }
    let live = true;
    const t = setTimeout(() => {
      void api.suggestPaths(path).then((r) => { if (live) setSuggestions(r.suggestions); }).catch(() => {});
      void api.validatePath(path).then((r) => { if (live) setValidation(r); }).catch(() => {});
    }, 200);
    return () => { live = false; clearTimeout(t); };
  }, [path]);

  const canSubmit = validation?.ok === true && duplicateOf === undefined && !submitting;

  const handleReview = async (): Promise<void> => {
    if (!canSubmit) return;
    setSubmitting(true);
    try {
      setDraft({ path: trimmedPath });
      navigate("/experiments/new/config");
    } finally {
      setSubmitting(false);
    }
  };

  const handleCancel = (): void => {
    clear();
    navigate("/experiments");
  };

  return (
    <PageFrame width="home" className="px-6 py-8">
      <button
        type="button"
        onClick={() => navigate("/experiments")}
        className="text-sm text-ink-soft hover:text-ink mb-4"
      >
        ← Experiments
      </button>
      <Kicker>New experiment</Kicker>
      <h1 className="text-display text-ink mb-1">Point at an experiment directory</h1>
      <p className="text-base text-ink-soft mb-6">
        Choose the directory of exposures. Himalaya reads the PRP and setup files,
        groups the frames into samples, and derives the geometry.
      </p>

      <div className="flex flex-col gap-5 max-w-[680px]">
        <div>
          <label htmlFor="dirpicker" className="block text-xs font-bold uppercase tracking-wide text-ink-faint mb-1.5">
            Data directory
          </label>
          <DirectoryPickerField
            value={path}
            onChange={setPath}
            suggestions={suggestions}
            validation={validation}
          />
          {duplicateOf !== undefined && (
            <p className="text-sm text-error mt-1.5" role="alert">
              This directory is already an experiment
              {duplicateOf.name ? ` ("${duplicateOf.name}")` : ""}. Each experiment is one directory.
            </p>
          )}
        </div>

        <div className="flex items-center gap-3">
          <Button variant="ghost" onClick={handleCancel}>Cancel</Button>
          <Button
            variant="accent"
            data-testid="create-submit"
            disabled={!canSubmit}
            onClick={() => { void handleReview(); }}
          >
            Review →
          </Button>
        </div>
      </div>
    </PageFrame>
  );
}
