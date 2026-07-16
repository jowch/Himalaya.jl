import { useEffect } from "react";
import { useNavigate, type NavigateFunction } from "react-router-dom";
import { useAppState } from "../../state";
import { suppressGlobalKeys } from "../../lib/keys";
import type { StaleUrlContext } from "../../state";

interface Props {
  staleUrlContext: StaleUrlContext;
}

interface VariantUI {
  dataMissing: string;
  header: string;
  ctaLabel: string;
  onPick: () => void;
}

function uiFor(
  ctx: StaleUrlContext,
  store: ReturnType<typeof useAppState.getState>,
  navigate: NavigateFunction,
): VariantUI {
  if (ctx.kind === "unknown_path") {
    return {
      dataMissing: "path",
      header: "Page not found.",
      ctaLabel: "Go to Experiments",
      // T3.2: unknown-path escape hatch now lands on /experiments (the home).
      // Replace mode so the back button doesn't loop back to /foo/bar.
      onPick: () => {
        navigate("/experiments", { replace: true });
      },
    };
  }
  if (ctx.missing === "experiment") {
    return {
      dataMissing: "experiment",
      header: `Experiment '${ctx.missing_value}' not found.`,
      ctaLabel: "Select an experiment",
      onPick: () => store.recoverFromStaleUrl({
        step: "experiment",
        experimentId: undefined,
        sampleId: undefined,
      }),
    };
  }
  if (ctx.missing === "sample") {
    const expName = ctx.experiment_resolved?.name ?? "?";
    const expId = ctx.experiment_resolved?.id;
    return {
      dataMissing: "sample",
      header: `Sample '${ctx.missing_value}' not found in '${expName}'.`,
      ctaLabel: "Select another sample",
      onPick: () => store.recoverFromStaleUrl({
        step: "sample",
        experimentId: expId,
        sampleId: undefined,
      }),
    };
  }
  // missing === "exposure"
  const sampleName = ctx.sample_resolved?.name ?? "?";
  const expId = ctx.experiment_resolved?.id;
  const sampleId = ctx.sample_resolved?.id;
  return {
    dataMissing: "exposure",
    header: `Exposure '${ctx.missing_value}' not found in '${sampleName}'.`,
    ctaLabel: "Back to sample",
    onPick: () => store.recoverFromStaleUrl({
      step: "sample",
      experimentId: expId,
      sampleId: sampleId,
      openModal: false,
    }),
  };
}

export function StaleUrlPage({ staleUrlContext }: Props): JSX.Element {
  const store = useAppState.getState();
  const navigate = useNavigate();
  const ui = uiFor(staleUrlContext, store, navigate);

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (suppressGlobalKeys(e)) return;
      if (e.key === "/" && !e.metaKey && !e.ctrlKey && !e.altKey) {
        e.preventDefault();
        ui.onPick();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [ui]);

  return (
    <div
      role="alert"
      data-testid="stale-url-page"
      data-missing={ui.dataMissing}
      className="flex-1 min-h-0 flex flex-col items-center justify-center p-8 text-center"
    >
      <h2 className="text-xl mb-2 text-ink">{ui.header}</h2>
      <p className="text-ink-soft mb-6">It may have been renamed or removed.</p>
      <button
        onClick={ui.onPick}
        data-testid="stale-url-cta"
        className="px-4 py-2 rounded border border-accent hover:bg-accent/10"
      >
        {ui.ctaLabel}
        <kbd className="ml-2 text-xs opacity-60">/</kbd>
      </button>
    </div>
  );
}
