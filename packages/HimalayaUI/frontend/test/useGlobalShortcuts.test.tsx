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

function corpusOf(...samples: Sample[]): CorpusSample[] {
  return samples as CorpusSample[];
}

function clientWithCorpus(corpus: CorpusSample[]): QueryClient {
  const client = makeClient();
  client.setQueryData(queryKeys.corpusSamples, corpus);
  return client;
}

// Renders the hook alongside a pathname probe so a test can observe whether a
// key navigated the URL (focus route) vs. only set store state.
function Harness(): JSX.Element {
  useGlobalShortcuts();
  const loc = useLocation();
  return <div data-testid="pathname">{loc.pathname}</div>;
}

// I5.1 (#182): the dual-nav `activePage` model + its ArrowLeft/Right page-tab
// step are deleted. The arrow keys are not bound by useGlobalShortcuts, so
// surfaces (Loupe frames, Focus exposures) own them outright via the shortcut
// library. These tests guard against a regression that re-binds them globally.
describe("useGlobalShortcuts — arrow keys are unbound (no page-tab step)", () => {
  beforeEach(() => {
    useAppState.setState({ activeSampleId: undefined });
  });

  it("does not mutate any nav state on a corpus route (the surface owns the arrows)", () => {
    const before = useAppState.getState();
    renderHook(() => useGlobalShortcuts(), {
      wrapper: wrapperAt("/samples/loupe/7"),
    });
    fireEvent.keyDown(document.body, { key: "ArrowRight" });
    fireEvent.keyDown(document.body, { key: "ArrowLeft" });
    const after = useAppState.getState();
    expect(after.activeSampleId).toBe(before.activeSampleId);
    expect(after.activeExperimentId).toBe(before.activeExperimentId);
  });

  it("is a no-op on a non-corpus route too (no global arrow binding remains)", () => {
    renderHook(() => useGlobalShortcuts(), {
      wrapper: wrapperAt("/totally/unknown"),
    });
    expect(() => {
      fireEvent.keyDown(document.body, { key: "ArrowRight" });
      fireEvent.keyDown(document.body, { key: "ArrowLeft" });
    }).not.toThrow();
  });

  // R0a (#221): the `T` theme-toggle shortcut is retired with the dark theme.
  it("does not bind T (theme toggle retired with the dark theme)", () => {
    renderHook(() => useGlobalShortcuts(), {
      wrapper: wrapperAt("/samples"),
    });
    const before = useAppState.getState();
    expect("theme" in (before as unknown as Record<string, unknown>)).toBe(false);
    expect(() => {
      fireEvent.keyDown(document.body, { key: "T" });
      fireEvent.keyDown(document.body, { key: "t" });
    }).not.toThrow();
    expect(useAppState.getState().activeSampleId).toBe(before.activeSampleId);
  });
});

// KEYS-LIB step 2: `,`/`.` is RETIRED as the sample stepper. `[`/`]` is now the
// one sample-step gesture, owned by the surfaces that have it (Loupe, Focus) via
// the shortcut library — NOT by this global hook. useGlobalShortcuts keeps only
// the find/nav-modal chord. These tests pin that `,`/`.` (and `[`/`]`) do
// nothing here: the focus-route case SEEDS the corpus so the retired code path
// WOULD have navigated — proving the binding is gone, not merely dormant.
describe("useGlobalShortcuts — ,/. sample step is retired (page owns [ ])", () => {
  const CORPUS = corpusOf(sample(10), sample(11), sample(12), sample(90, 2));

  beforeEach(() => {
    useAppState.setState({ activeSampleId: 11 });
  });

  it("does not navigate on '.' from /sample/:id even with a seeded corpus", () => {
    render(<Harness />, {
      wrapper: wrapperAt("/sample/11", clientWithCorpus(CORPUS)),
    });
    fireEvent.keyDown(document.body, { key: "." });
    expect(screen.getByTestId("pathname")).toHaveTextContent("/sample/11");
  });

  it("does not navigate on ',' from /sample/:id even with a seeded corpus", () => {
    render(<Harness />, {
      wrapper: wrapperAt("/sample/11", clientWithCorpus(CORPUS)),
    });
    fireEvent.keyDown(document.body, { key: "," });
    expect(screen.getByTestId("pathname")).toHaveTextContent("/sample/11");
  });

  it("does not change the active sample on a corpus route via '.' or ','", () => {
    render(<Harness />, { wrapper: wrapperAt("/samples", clientWithCorpus(CORPUS)) });
    fireEvent.keyDown(document.body, { key: "." });
    fireEvent.keyDown(document.body, { key: "," });
    expect(useAppState.getState().activeSampleId).toBe(11);
  });

  it("does not bind [ or ] globally (the sample step is page-owned, not global)", () => {
    render(<Harness />, {
      wrapper: wrapperAt("/sample/11", clientWithCorpus(CORPUS)),
    });
    fireEvent.keyDown(document.body, { key: "[" });
    fireEvent.keyDown(document.body, { key: "]" });
    // No navigation and no store step from the global hook.
    expect(screen.getByTestId("pathname")).toHaveTextContent("/sample/11");
    expect(useAppState.getState().activeSampleId).toBe(11);
  });
});

// The find/jump chord is the ONE binding useGlobalShortcuts keeps. LO-KEYD: the
// shared suppressGlobalKeys guard must NOT swallow modifier chords — ⌘K is the
// one binding here that REQUIRES metaKey to pass through.
describe("useGlobalShortcuts — find/nav-modal chord (⌘K and /)", () => {
  beforeEach(() => {
    useAppState.setState({ navModalOpen: false, activeExperimentId: undefined });
  });

  it("opens the nav modal on ⌘K from a non-editing target", () => {
    renderHook(() => useGlobalShortcuts(), { wrapper: wrapperAt("/samples") });
    fireEvent.keyDown(document.body, { key: "k", metaKey: true });
    expect(useAppState.getState().navModalOpen).toBe(true);
  });

  it("opens the nav modal on '/' from a non-editing target", () => {
    renderHook(() => useGlobalShortcuts(), { wrapper: wrapperAt("/samples") });
    fireEvent.keyDown(document.body, { key: "/" });
    expect(useAppState.getState().navModalOpen).toBe(true);
  });
});
