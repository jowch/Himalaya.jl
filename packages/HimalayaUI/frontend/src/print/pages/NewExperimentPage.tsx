import { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import { useAppState } from "../../state";
import { useExperiments } from "../../queries";
import * as api from "../../api";
import { authOpts } from "../../lib/authOpts";
import { getClientId } from "../../lib/clientId";
import { PageFrame } from "../components/PageFrame";
import { DirectoryPickerField } from "../components/DirectoryPickerField";
import { Button } from "../ui/Button";
import { Kicker } from "../ui/Kicker";
import { Input } from "../ui/Input";
import type { ValidatePathResponse } from "../../api";

const CLIENT_ID = getClientId();

/**
 * NewExperimentPage — /experiments/new (spec §8.7). Directory picker +
 * suggestions + validation + optional name + a File-patterns advanced
 * disclosure; submit creates the experiment and routes to its corpus, where the
 * live-ingest unfold runs. No file browser.
 */
export function NewExperimentPage(): JSX.Element {
  const navigate = useNavigate();
  const username = useAppState((s) => s.username);
  const setActiveExperiment = useAppState((s) => s.setActiveExperiment);

  const { data: experiments } = useExperiments();

  const [path, setPath] = useState("");
  const [name, setName] = useState("");
  const [suggestions, setSuggestions] = useState<string[]>([]);
  const [validation, setValidation] = useState<ValidatePathResponse | null>(null);
  const [submitting, setSubmitting] = useState(false);
  const [serverError, setServerError] = useState<string | null>(null);

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

  const submit = async (): Promise<void> => {
    if (!canSubmit) return;
    setSubmitting(true);
    setServerError(null);
    try {
      const exp = await api.createExperiment(
        { path, ...(name.trim() ? { name: name.trim() } : {}) },
        authOpts(username, CLIENT_ID),
      );
      setActiveExperiment(exp.id);
      navigate(`/experiments/${exp.id}/corpus`);
    } catch (err) {
      // Server-side fallback for the one-experiment-per-directory rule (409).
      setServerError(
        err instanceof Error && /409|already uses this directory/i.test(err.message)
          ? "An experiment already uses this directory."
          : "Could not create the experiment. Please try again.",
      );
    } finally {
      setSubmitting(false);
    }
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
              {duplicateOf.name ? ` (“${duplicateOf.name}”)` : ""}. Each experiment is one directory.
            </p>
          )}
          {serverError !== null && (
            <p className="text-sm text-error mt-1.5" role="alert">{serverError}</p>
          )}
        </div>

        <div>
          <label htmlFor="exp-name" className="block text-xs font-bold uppercase tracking-wide text-ink-faint mb-1.5">
            Experiment name <span className="font-normal normal-case text-ink-soft">(optional)</span>
          </label>
          <Input
            testId="exp-name"
            value={name}
            onValueChange={setName}
            placeholder="Defaults to the directory name"
            aria-label="Experiment name"
          />
        </div>

        <div className="flex items-center gap-3">
          <Button variant="ghost" onClick={() => navigate("/experiments")}>Cancel</Button>
          <Button
            variant="accent"
            data-testid="create-submit"
            disabled={!canSubmit}
            onClick={() => { void submit(); }}
          >
            Scan and create
          </Button>
        </div>
      </div>
    </PageFrame>
  );
}
