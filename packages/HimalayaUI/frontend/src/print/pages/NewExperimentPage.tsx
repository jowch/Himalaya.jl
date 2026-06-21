import { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import { useExperiments } from "../../queries";
import * as api from "../../api";
import { useDraftExperiment } from "../../lib/draftExperiment";
import { PageFrame } from "../components/PageFrame";
import { DirectoryPickerField } from "../components/DirectoryPickerField";
import { Card } from "../ui/Card";
import { Button } from "../ui/Button";
import { Kicker } from "../ui/Kicker";
import type { ValidatePathResponse } from "../../api";

/** One inline pre-flight check in the dual ✓ row (mockup p1-new): a green ✓
 *  when satisfied, a red ✗ + reason when not, faint while pending. */
function PreflightCheck({
  state,
  pass,
  fail,
}: {
  state: "pass" | "fail" | "pending";
  pass: string;
  fail: string;
}): JSX.Element {
  if (state === "pending") {
    return <span className="text-meta text-ink-faint">◦ {pass}</span>;
  }
  if (state === "pass") {
    return <span className="text-meta font-semibold text-success" role="status">✓ {pass}</span>;
  }
  return <span className="text-meta font-semibold text-error" role="alert">✗ {fail}</span>;
}

/**
 * NewExperimentPage — /experiments/new (spec §8.7, mockup p1-new). Directory
 * picker + suggestions + two pre-flight checks (directory exists · not already
 * an experiment). The primary action commits the validated path to a
 * client-side draft and navigates to /experiments/new/config (first-run
 * Configuration). No DB row is created here — creation happens at Approve.
 *
 * Two-phase funnel handoff (T4.0): "Review →" replaces the old "Scan and
 * create". The dup-dir guard (inline ✗ + disabled submit) is preserved.
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

  // Pre-flight check states. Both are "pending" until a path is entered.
  const existsState: "pass" | "fail" | "pending" =
    trimmedPath === "" || validation == null ? "pending" : validation.ok ? "pass" : "fail";
  const uniqueState: "pass" | "fail" | "pending" =
    trimmedPath === "" ? "pending" : duplicateOf === undefined ? "pass" : "fail";

  return (
    <>
      <PageFrame width="home" className="px-6 py-8 pb-28">
        <div className="max-w-[760px]">
          <Kicker tone="accent">New experiment</Kicker>
          <h1 className="text-display text-ink">Point at a directory</h1>
          <p className="text-body text-ink-soft mt-2 mb-6 max-w-[60ch]">
            Pick the folder for this experiment. The next step indexes it and lets
            you review before anything is created.
          </p>

          <Card padding="lg">
            <Kicker tone="soft" className="mb-2">Directory</Kicker>
            <DirectoryPickerField
              value={path}
              onChange={setPath}
              suggestions={suggestions}
              validation={null}
            />
            <div
              className="mt-3 flex flex-wrap items-center gap-x-6 gap-y-1.5"
              data-testid="dirpicker-checks"
            >
              <PreflightCheck
                state={existsState}
                pass="directory exists"
                fail={validation?.message ?? "directory not found"}
              />
              <PreflightCheck
                state={uniqueState}
                pass="not already an experiment"
                fail={
                  duplicateOf?.name
                    ? `already an experiment ("${duplicateOf.name}")`
                    : "already an experiment"
                }
              />
            </div>
          </Card>
        </div>
      </PageFrame>

      {/* Sticky funnel footer (p1-new): reassurance + Cancel/Review. */}
      <footer
        data-testid="funnel-footer"
        className="fixed bottom-0 left-0 right-0 z-30 flex items-center justify-between gap-4 border-t border-hair bg-paper px-8 py-3"
      >
        <span className="text-meta text-ink-soft">
          Nothing is created yet. The next step indexes and lets you review.
        </span>
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
      </footer>
    </>
  );
}
