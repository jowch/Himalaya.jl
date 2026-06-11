import { describe, it, expect, beforeEach } from "vitest";
import { renderHook, render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter, useLocation } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { useGlobalShortcuts } from "../src/hooks/useGlobalShortcuts";
import { useAppState } from "../src/state";
import { queryKeys } from "../src/queries";
import type { CorpusSample, Sample } from "../src/api";

// staleTime: Infinity so a setQueryData-seeded corpus list is served without a
// mount refetch (no fetch mock needed); retry: false keeps cold-cache tests
// from spinning.
function makeClient(): QueryClient {
  return new QueryClient({
    defaultOptions: { queries: { retry: false, staleTime: Infinity, gcTime: Infinity } },
  });
}

function wrapperAt(path: string, client: QueryClient = makeClient()) {
  return function Wrapper({ children }: { children: ReactNode }): JSX.Element {
    return (
      <QueryClientProvider client={client}>
        <MemoryRouter initialEntries={[path]}>{children}</MemoryRouter>
      </QueryClientProvider>
    );
  };
}

function sample(id: number, experiment_id = 1): Sample {
  return { id, experiment_id, name: `S${id}`, display_name: null,
    notes: null, tags: [], q_units: "A-1" } as Sample;
}

// The corpus list the focus-route step derives its siblings from (F5): same
// rows as the per-experiment list, in corpus order, as CorpusSample.
function corpusOf(...samples: Sample[]): CorpusSample[] {
  return samples as CorpusSample[];
}

function clientWithCorpus(corpus: CorpusSample[]): QueryClient {
  const client = makeClient();
  client.setQueryData(queryKeys.corpusSamples, corpus);
  return client;
}

// Renders the hook alongside a pathname probe so a test can observe whether the
// ,/. shortcut navigated the URL (focus route) vs. only set store state.
function Harness({ samples }: { samples: Sample[] | undefined }): JSX.Element {
  useGlobalShortcuts(samples);
  const loc = useLocation();
  return <div data-testid="pathname">{loc.pathname}</div>;
}

// I5.1 (#182): the dual-nav `activePage` model + its ArrowLeft/Right page-tab
// step are deleted. The arrow keys are no longer bound by useGlobalShortcuts,
// so corpus surfaces (e.g. the loupe at /samples/loupe/:id) own them outright.
// These tests guard against a regression that re-binds them globally.
describe("useGlobalShortcuts — arrow keys are unbound (no page-tab step)", () => {
  beforeEach(() => {
    useAppState.setState({ activeSampleId: undefined });
  });

  it("does not mutate any nav state on a corpus route (loupe owns the arrows)", () => {
    const before = useAppState.getState();
    renderHook(() => useGlobalShortcuts(undefined), {
      wrapper: wrapperAt("/samples/loupe/7"),
    });
    fireEvent.keyDown(document.body, { key: "ArrowRight" });
    fireEvent.keyDown(document.body, { key: "ArrowLeft" });
    const after = useAppState.getState();
    // The hook must not touch the active ids on an arrow press.
    expect(after.activeSampleId).toBe(before.activeSampleId);
    expect(after.activeExperimentId).toBe(before.activeExperimentId);
  });

  it("is a no-op on a non-corpus route too (no global arrow binding remains)", () => {
    renderHook(() => useGlobalShortcuts(undefined), {
      wrapper: wrapperAt("/totally/unknown"),
    });
    // Should neither throw nor mutate state.
    expect(() => {
      fireEvent.keyDown(document.body, { key: "ArrowRight" });
      fireEvent.keyDown(document.body, { key: "ArrowLeft" });
    }).not.toThrow();
  });

  // R0a (#221): the `T` theme-toggle shortcut is retired with the dark theme.
  // Pressing T must no longer mutate any store state (no `theme` field exists).
  it("does not bind T (theme toggle retired with the dark theme)", () => {
    renderHook(() => useGlobalShortcuts(undefined), {
      wrapper: wrapperAt("/samples"),
    });
    const before = useAppState.getState();
    expect("theme" in (before as Record<string, unknown>)).toBe(false);
    expect(() => {
      fireEvent.keyDown(document.body, { key: "T" });
      fireEvent.keyDown(document.body, { key: "t" });
    }).not.toThrow();
    // No store mutation from T.
    expect(useAppState.getState().activeSampleId).toBe(before.activeSampleId);
  });
});

// F5: on the focus route the step derives the active sample's
// experiment-siblings from the CORPUS cache (the same shared derivation the
// topbar stepper uses) — NOT from the `samplesInExperiment` parameter, which
// is undefined whenever `activeExperimentId` was never set (direct visit /
// door entry; only the NavModal picker and recoverFromStaleUrl write it).
describe("useGlobalShortcuts — ,/. sample step (focus route, corpus-derived)", () => {
  const CORPUS = corpusOf(sample(10), sample(11), sample(12), sample(90, 2));

  beforeEach(() => {
    useAppState.setState({ activeSampleId: 11 });
  });

  it("navigates to the next sibling URL on '.' from /sample/:id with NO samplesInExperiment", () => {
    render(<Harness samples={undefined} />, {
      wrapper: wrapperAt("/sample/11", clientWithCorpus(CORPUS)),
    });
    fireEvent.keyDown(document.body, { key: "." });
    expect(screen.getByTestId("pathname")).toHaveTextContent("/sample/12");
  });

  it("navigates to the previous sibling URL on ',' from /sample/:id with NO samplesInExperiment", () => {
    render(<Harness samples={undefined} />, {
      wrapper: wrapperAt("/sample/11", clientWithCorpus(CORPUS)),
    });
    fireEvent.keyDown(document.body, { key: "," });
    expect(screen.getByTestId("pathname")).toHaveTextContent("/sample/10");
  });

  it("steps only within the active sample's experiment (corpus order, no cross-experiment leak)", () => {
    useAppState.setState({ activeSampleId: 12 });
    render(<Harness samples={undefined} />, {
      wrapper: wrapperAt("/sample/12", clientWithCorpus(CORPUS)),
    });
    // sample 12 is the last of experiment 1; sample 90 (experiment 2) follows
    // it in corpus order but is NOT a sibling — the step must not wrap or leak.
    fireEvent.keyDown(document.body, { key: "." });
    expect(screen.getByTestId("pathname")).toHaveTextContent("/sample/12");
  });

  it("does not wrap past the first sample on the focus route", () => {
    useAppState.setState({ activeSampleId: 10 });
    render(<Harness samples={undefined} />, {
      wrapper: wrapperAt("/sample/10", clientWithCorpus(CORPUS)),
    });
    fireEvent.keyDown(document.body, { key: "," });
    expect(screen.getByTestId("pathname")).toHaveTextContent("/sample/10");
  });

  it("no-ops gracefully when the corpus cache is cold", () => {
    render(<Harness samples={undefined} />, {
      wrapper: wrapperAt("/sample/11", makeClient()),
    });
    expect(() => {
      fireEvent.keyDown(document.body, { key: "." });
      fireEvent.keyDown(document.body, { key: "," });
    }).not.toThrow();
    expect(screen.getByTestId("pathname")).toHaveTextContent("/sample/11");
  });

  it("no-ops gracefully when the active sample is unknown to the corpus", () => {
    useAppState.setState({ activeSampleId: 99999 });
    render(<Harness samples={undefined} />, {
      wrapper: wrapperAt("/sample/99999", clientWithCorpus(CORPUS)),
    });
    fireEvent.keyDown(document.body, { key: "." });
    expect(screen.getByTestId("pathname")).toHaveTextContent("/sample/99999");
  });

  it("F-STALEURL: bogus /sample/:id URL with a stale VALID store sample -> '.' no-ops", () => {
    // Mid-session navigation to /sample/99999 never seeds the store, so the
    // previous valid activeSampleId (11) survives. The topbar hides its
    // stepper via resolveRouteSampleStatus; the shortcut must gate on the
    // SAME predicate — otherwise it would invisibly step relative to a sample
    // the URL does not show, out of a "Sample not found" page.
    useAppState.setState({ activeSampleId: 11 });
    render(<Harness samples={undefined} />, {
      wrapper: wrapperAt("/sample/99999", clientWithCorpus(CORPUS)),
    });
    fireEvent.keyDown(document.body, { key: "." });
    expect(screen.getByTestId("pathname")).toHaveTextContent("/sample/99999");
    fireEvent.keyDown(document.body, { key: "," });
    expect(screen.getByTestId("pathname")).toHaveTextContent("/sample/99999");
  });

  it("F-STALEURL: non-numeric /sample/:id with a stale valid store sample -> '.' no-ops", () => {
    useAppState.setState({ activeSampleId: 11 });
    render(<Harness samples={undefined} />, {
      wrapper: wrapperAt("/sample/not-a-number", clientWithCorpus(CORPUS)),
    });
    fireEvent.keyDown(document.body, { key: "." });
    expect(screen.getByTestId("pathname")).toHaveTextContent("/sample/not-a-number");
  });

  it("ignores the samplesInExperiment parameter on the focus route (corpus wins)", () => {
    // A stale per-experiment list (e.g. from a previously-picked experiment)
    // must not drive the focus-route step — the corpus-derived siblings do.
    const staleOtherExperiment = [sample(90, 2), sample(91, 2)];
    render(<Harness samples={staleOtherExperiment} />, {
      wrapper: wrapperAt("/sample/11", clientWithCorpus(CORPUS)),
    });
    fireEvent.keyDown(document.body, { key: "." });
    expect(screen.getByTestId("pathname")).toHaveTextContent("/sample/12");
  });
});

// Corpus surfaces keep the store-driven, samplesInExperiment-parameterized step.
describe("useGlobalShortcuts — ,/. sample step (corpus surfaces, store-driven)", () => {
  const SAMPLES = [sample(10), sample(11), sample(12)];

  beforeEach(() => {
    useAppState.setState({ activeSampleId: 11 });
  });

  it("on a corpus route, '.' sets store state without navigating", () => {
    render(<Harness samples={SAMPLES} />, { wrapper: wrapperAt("/samples") });
    fireEvent.keyDown(document.body, { key: "." });
    // Store advanced, URL unchanged (corpus surfaces are store-driven).
    expect(useAppState.getState().activeSampleId).toBe(12);
    expect(screen.getByTestId("pathname")).toHaveTextContent("/samples");
  });

  it("on a corpus route, ',' steps the store back without navigating", () => {
    render(<Harness samples={SAMPLES} />, { wrapper: wrapperAt("/samples") });
    fireEvent.keyDown(document.body, { key: "," });
    expect(useAppState.getState().activeSampleId).toBe(10);
    expect(screen.getByTestId("pathname")).toHaveTextContent("/samples");
  });
});

// LO-KEYD: the shared suppressGlobalKeys guard must NOT swallow modifier
// chords — ⌘K is the one binding here that REQUIRES metaKey to pass through.
describe("useGlobalShortcuts — ⌘K passes the suppression guard", () => {
  beforeEach(() => {
    useAppState.setState({ navModalOpen: false, activeExperimentId: undefined });
  });

  it("opens the nav modal on ⌘K from a non-editing target", () => {
    renderHook(() => useGlobalShortcuts(undefined), {
      wrapper: wrapperAt("/samples"),
    });
    fireEvent.keyDown(document.body, { key: "k", metaKey: true });
    expect(useAppState.getState().navModalOpen).toBe(true);
  });
});
