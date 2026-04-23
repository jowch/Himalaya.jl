# HimalayaUI Plan 3.5 — Frontend Hardening Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix four real structural weaknesses in the Plan 3 frontend before Plans 4–6 compound them: server-state handling, styling system, duplicated username state, and missing top-level error handling.

**Architecture:** Four independent, commit-per-task changes. Each task stands alone and leaves the app fully working and tested.

1. **Username single-source** — delete the `_username` module-global in `api.ts`; read from Zustand at call sites.
2. **Error boundary** — wrap `<App />` in `<ErrorBoundary>` so a render error shows a message instead of a blank page.
3. **TanStack Query** — replace the manual `useEffect(() => fetch())` pattern with `useQuery` so Plans 4–6 get caching, refetch, and mutation invalidation for free.
4. **Tailwind v4 + CSS-var theme** — migrate `styles.css` class-by-class to Tailwind utilities, keeping the existing dark palette as CSS variables on `:root`.

**Tech Stack:** React 18, TypeScript 5.6 (strict, `exactOptionalPropertyTypes: true`), Zustand 4.5, **TanStack Query 5** (new), **Tailwind v4** (new), Vitest + React Testing Library, Playwright.

**Assumptions about the starting state** (verified at plan-write time):

- `packages/HimalayaUI/frontend/package.json` deps: `react`, `react-dom`, `zustand`. Scripts: `dev`, `build` (`tsc --noEmit -p tsconfig.build.json && vite build`), `test`, `test:watch`, `e2e`.
- `src/` layout: `App.tsx`, `api.ts`, `state.ts`, `main.tsx`, `styles.css`, `components/{Navbar,Layout,SampleList,UserModal}.tsx`.
- `test/` at `frontend/test/`: `SampleList.test.tsx`, `api.test.ts`, `state.test.ts`, `smoke.test.tsx`, `setup.ts`.
- `e2e/smoke.spec.ts` mocks `/api/users`, `/api/experiments/1`, `/api/experiments/1/samples` via `page.route`.
- `App.tsx` imports `./styles.css` and uses class names like `nav`, `layout`, `col-left`, `sample-row`, `modal-dialog`.
- `main.tsx` renders `<StrictMode><App /></StrictMode>` into `#app`.
- All tests currently pass; build is clean.

**Run all commands from:** `packages/HimalayaUI/frontend/`.

---

## Task 1: Username single-source via Zustand

**Motivation:** `api.ts` keeps its own `_username` module-level variable, synced from `App.tsx` via a `useEffect`. Two copies of the same string drift easily. Read from Zustand at the one call site instead.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/api.ts` (remove `_username`/`setUsername`/`getUsername`; accept optional `username` on mutating calls)
- Modify: `packages/HimalayaUI/frontend/src/App.tsx` (drop the sync effect; no longer passes username into api)
- Modify: `packages/HimalayaUI/frontend/src/components/UserModal.tsx` (drop `api.setUsername`; `onSelect` is the only handoff)
- Modify: `packages/HimalayaUI/frontend/test/api.test.ts` (update expectations)

### Approach

Change `api.ts` mutating functions (`createUser`, `updateExperiment`, `updateSample`, `addSampleTag`, `removeSampleTag`) to accept an optional `opts: { username?: string }` parameter. The `X-Username` header is attached from `opts.username`, not from module state. Call sites read `useAppState.getState().username` (Zustand outside React) or `useAppState((s) => s.username)` (inside React).

### Steps

- [ ] **Step 1: Write failing test for api without module state**

Edit `packages/HimalayaUI/frontend/test/api.test.ts`. Remove any existing test that uses `setUsername`/`getUsername` and replace the "sends X-Username header" test with:

```ts
import { describe, it, expect, vi, beforeEach } from "vitest";
import * as api from "../src/api";

describe("api", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("sends X-Username header on mutating calls when opts.username provided", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({ id: 1, username: "alice" }), { status: 200 }),
    );
    await api.createUser("alice", { username: "alice" });
    const init = fetchSpy.mock.calls[0]![1] as RequestInit;
    expect((init.headers as Record<string, string>)["X-Username"]).toBe("alice");
  });

  it("omits X-Username on GET requests even when opts passed", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response("[]", { status: 200 }),
    );
    await api.listUsers();
    const init = fetchSpy.mock.calls[0]![1] as RequestInit;
    const headers = (init.headers as Record<string, string>) ?? {};
    expect(headers["X-Username"]).toBeUndefined();
  });

  it("does not export setUsername/getUsername", () => {
    expect((api as unknown as Record<string, unknown>).setUsername).toBeUndefined();
    expect((api as unknown as Record<string, unknown>).getUsername).toBeUndefined();
  });
});
```

- [ ] **Step 2: Run test to verify failure**

Run: `npm test -- api.test`

Expected: FAIL — `createUser` currently takes one arg; `setUsername` still exported.

- [ ] **Step 3: Rewrite `api.ts`**

Replace `src/api.ts` contents:

```ts
export interface User { id: number; username: string; }

export interface Experiment {
  id: number;
  name: string | null;
  path: string;
  data_dir: string;
  analysis_dir: string;
  manifest_path: string | null;
  created_at: string;
}

export interface SampleTag {
  id: number;
  key: string;
  value: string;
  source: string;
}

export interface Sample {
  id: number;
  experiment_id: number;
  label: string | null;
  name: string | null;
  notes: string | null;
  tags: SampleTag[];
}

export class ApiError extends Error {
  constructor(public status: number, message: string, public body: unknown) {
    super(message);
    this.name = "ApiError";
  }
}

export interface AuthOpts { username?: string }

async function request<T>(
  method: string,
  path: string,
  body?: unknown,
  opts?: AuthOpts,
): Promise<T> {
  const headers: Record<string, string> = {};
  if (body !== undefined) headers["Content-Type"] = "application/json";
  if (opts?.username && method !== "GET") headers["X-Username"] = opts.username;

  const init: RequestInit = { method, headers };
  if (body !== undefined) init.body = JSON.stringify(body);

  const res = await fetch(path, init);
  if (res.status === 204) return undefined as T;

  const text = await res.text();
  const parsed = text ? safeJson(text) : null;
  if (!res.ok) {
    const msg = parsed && typeof parsed === "object" && parsed !== null && "error" in parsed
      ? String((parsed as { error: unknown }).error)
      : `${method} ${path} failed with ${res.status}`;
    throw new ApiError(res.status, msg, parsed);
  }
  return parsed as T;
}

function safeJson(s: string): unknown {
  try { return JSON.parse(s); } catch { return s; }
}

// Users
export const listUsers  = () => request<User[]>("GET", "/api/users");
export const createUser = (username: string, opts?: AuthOpts) =>
  request<User>("POST", "/api/users", { username }, opts);

// Experiments
export const getExperiment = (id: number) =>
  request<Experiment>("GET", `/api/experiments/${id}`);
export const updateExperiment = (
  id: number,
  patch: Partial<Pick<Experiment, "name" | "data_dir" | "analysis_dir" | "manifest_path">>,
  opts?: AuthOpts,
) => request<Experiment>("PATCH", `/api/experiments/${id}`, patch, opts);

// Samples
export const listSamples    = (experiment_id: number) =>
  request<Sample[]>("GET", `/api/experiments/${experiment_id}/samples`);
export const updateSample   = (id: number, patch: { name?: string; notes?: string }, opts?: AuthOpts) =>
  request<Sample>("PATCH", `/api/samples/${id}`, patch, opts);
export const addSampleTag   = (id: number, key: string, value: string, opts?: AuthOpts) =>
  request<SampleTag>("POST", `/api/samples/${id}/tags`, { key, value }, opts);
export const removeSampleTag = (id: number, tag_id: number, opts?: AuthOpts) =>
  request<void>("DELETE", `/api/samples/${id}/tags/${tag_id}`, undefined, opts);
```

- [ ] **Step 4: Run api test to verify pass**

Run: `npm test -- api.test`
Expected: PASS (all three tests green).

- [ ] **Step 5: Update `App.tsx` — drop the sync effect**

Edit `packages/HimalayaUI/frontend/src/App.tsx`. Remove the line `useEffect(() => { api.setUsername(username); }, [username]);` and its now-unused imports. App continues to read `username` from Zustand for display in the navbar, and passes it to UserModal via `onSelect` as before. The full file after edit:

```tsx
import { useEffect, useState } from "react";
import "./styles.css";
import * as api from "./api";
import type { Experiment, Sample } from "./api";
import { useAppState } from "./state";
import { Navbar } from "./components/Navbar";
import { Layout } from "./components/Layout";
import { SampleList } from "./components/SampleList";
import { UserModal } from "./components/UserModal";

const EXPERIMENT_ID = 1;

export function App(): JSX.Element {
  const username        = useAppState((s) => s.username);
  const activeSampleId  = useAppState((s) => s.activeSampleId);
  const setUsername     = useAppState((s) => s.setUsername);
  const setActiveSample = useAppState((s) => s.setActiveSample);

  const [experiment, setExperiment] = useState<Experiment | null>(null);
  const [samples, setSamples]       = useState<Sample[]>([]);
  const [bootError, setBootError]   = useState<string | null>(null);
  const [modalOpen, setModalOpen]   = useState<boolean>(!username);

  useEffect(() => {
    void (async () => {
      try {
        const [exp, list] = await Promise.all([
          api.getExperiment(EXPERIMENT_ID),
          api.listSamples(EXPERIMENT_ID),
        ]);
        setExperiment(exp);
        setSamples(list);
      } catch (e) {
        setBootError((e as Error).message);
      }
    })();
  }, []);

  const activeSample = samples.find((s) => s.id === activeSampleId);
  const breadcrumb   = bootError
    ? `Error: ${bootError}`
    : (experiment?.name ?? "experiment")
        + (activeSample
          ? ` › ${activeSample.label ?? ""} ${activeSample.name ?? ""}`.trimEnd()
          : "");

  return (
    <>
      <Navbar
        breadcrumb={breadcrumb}
        username={username}
        onUserClick={() => setModalOpen(true)}
      />
      <Layout
        left={
          <SampleList
            samples={samples}
            activeId={activeSampleId}
            onSelect={setActiveSample}
          />
        }
      />
      <UserModal
        open={modalOpen}
        onSelect={(name) => {
          setUsername(name);
          setModalOpen(false);
        }}
        onClose={() => setModalOpen(false)}
      />
    </>
  );
}
```

- [ ] **Step 6: Update `UserModal.tsx` — pass username explicitly to createUser**

Edit `packages/HimalayaUI/frontend/src/components/UserModal.tsx`. Remove the `api.setUsername(username)` line inside `submit()`. Pass the new name to `createUser` via opts. Replacement `submit`:

```tsx
async function submit(): Promise<void> {
  setError(null);
  const username = selection === "__new__" ? newName.trim() : selection;
  if (!username) {
    setError("Username required");
    return;
  }
  try {
    await api.createUser(username, { username });
    onSelect(username);
  } catch (e) {
    setError(`Failed: ${(e as Error).message}`);
  }
}
```

- [ ] **Step 7: Run full test suite + build**

Run: `npm test && npm run build`
Expected: all unit tests pass; `tsc --noEmit` clean; `vite build` produces `dist/`.

- [ ] **Step 8: Commit**

```bash
cd packages/HimalayaUI/frontend
git add src/api.ts src/App.tsx src/components/UserModal.tsx test/api.test.ts
git commit -m "refactor(ui): drop _username module global; pass username per-call"
```

---

## Task 2: Top-level error boundary

**Motivation:** A render-time exception in any component currently blanks the page. A tiny class component keeps Plans 4–6 safe from regressions.

**Files:**
- Create: `packages/HimalayaUI/frontend/src/ErrorBoundary.tsx`
- Create: `packages/HimalayaUI/frontend/test/ErrorBoundary.test.tsx`
- Modify: `packages/HimalayaUI/frontend/src/main.tsx`

### Steps

- [ ] **Step 1: Write failing test**

Create `packages/HimalayaUI/frontend/test/ErrorBoundary.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { ErrorBoundary } from "../src/ErrorBoundary";

function Boom(): JSX.Element { throw new Error("kaboom"); }

describe("ErrorBoundary", () => {
  let errSpy: ReturnType<typeof vi.spyOn>;
  beforeEach(() => { errSpy = vi.spyOn(console, "error").mockImplementation(() => {}); });
  afterEach(() => { errSpy.mockRestore(); });

  it("renders children when they do not throw", () => {
    render(<ErrorBoundary><p>ok</p></ErrorBoundary>);
    expect(screen.getByText("ok")).toBeInTheDocument();
  });

  it("renders fallback when a child throws", () => {
    render(<ErrorBoundary><Boom /></ErrorBoundary>);
    expect(screen.getByRole("alert")).toHaveTextContent(/kaboom/i);
  });
});
```

- [ ] **Step 2: Run test to verify failure**

Run: `npm test -- ErrorBoundary`
Expected: FAIL — `ErrorBoundary` does not exist.

- [ ] **Step 3: Implement `ErrorBoundary.tsx`**

Create `packages/HimalayaUI/frontend/src/ErrorBoundary.tsx`:

```tsx
import { Component, type ErrorInfo, type ReactNode } from "react";

interface Props { children: ReactNode }
interface State { error: Error | null }

export class ErrorBoundary extends Component<Props, State> {
  state: State = { error: null };

  static getDerivedStateFromError(error: Error): State { return { error }; }

  componentDidCatch(error: Error, info: ErrorInfo): void {
    console.error("ErrorBoundary caught:", error, info);
  }

  render(): ReactNode {
    if (this.state.error) {
      return (
        <div role="alert" style={{ padding: 24, fontFamily: "monospace" }}>
          <h2>Something went wrong.</h2>
          <pre>{this.state.error.message}</pre>
          <button onClick={() => this.setState({ error: null })}>Try again</button>
        </div>
      );
    }
    return this.props.children;
  }
}
```

- [ ] **Step 4: Run test to verify pass**

Run: `npm test -- ErrorBoundary`
Expected: PASS (2/2).

- [ ] **Step 5: Wrap App in `main.tsx`**

Replace `packages/HimalayaUI/frontend/src/main.tsx`:

```tsx
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { App } from "./App";
import { ErrorBoundary } from "./ErrorBoundary";

const root = document.getElementById("app");
if (!root) throw new Error("#app root missing");
createRoot(root).render(
  <StrictMode>
    <ErrorBoundary>
      <App />
    </ErrorBoundary>
  </StrictMode>,
);
```

- [ ] **Step 6: Run full suite + build**

Run: `npm test && npm run build`
Expected: all green.

- [ ] **Step 7: Commit**

```bash
cd packages/HimalayaUI/frontend
git add src/ErrorBoundary.tsx src/main.tsx test/ErrorBoundary.test.tsx
git commit -m "feat(ui): add top-level ErrorBoundary"
```

---

## Task 3: TanStack Query for server state

**Motivation:** Plans 4–6 add exposures, peaks, and indices — all fetched, all mutated. The current `useEffect(() => fetch())` pattern doesn't cache, doesn't refetch on mutation, and doesn't dedupe. TanStack Query is the standard solution and plays nicely with Zustand (Zustand keeps client state like `activeSampleId`; Query keeps server state).

**Files:**
- Modify: `packages/HimalayaUI/frontend/package.json` (add `@tanstack/react-query`)
- Modify: `packages/HimalayaUI/frontend/src/main.tsx` (wrap in `QueryClientProvider`)
- Create: `packages/HimalayaUI/frontend/src/queries.ts` (query-key + hook definitions)
- Modify: `packages/HimalayaUI/frontend/src/App.tsx` (replace `useEffect` with `useQuery`)
- Modify: `packages/HimalayaUI/frontend/test/smoke.test.tsx` (wrap in provider)
- Create: `packages/HimalayaUI/frontend/test/test-utils.tsx` (shared `renderWithProviders`)

### Steps

- [ ] **Step 1: Install `@tanstack/react-query`**

```bash
cd packages/HimalayaUI/frontend
npm install @tanstack/react-query@^5.59.0
```

Expected: `package.json` and `package-lock.json` updated; no install errors.

- [ ] **Step 2: Create `queries.ts`**

Create `packages/HimalayaUI/frontend/src/queries.ts`:

```ts
import { useQuery } from "@tanstack/react-query";
import * as api from "./api";

export const queryKeys = {
  experiment: (id: number) => ["experiment", id] as const,
  samples:    (experimentId: number) => ["experiment", experimentId, "samples"] as const,
};

export function useExperiment(id: number) {
  return useQuery({
    queryKey: queryKeys.experiment(id),
    queryFn: () => api.getExperiment(id),
  });
}

export function useSamples(experimentId: number) {
  return useQuery({
    queryKey: queryKeys.samples(experimentId),
    queryFn: () => api.listSamples(experimentId),
  });
}
```

- [ ] **Step 3: Create shared test helper `test-utils.tsx`**

Create `packages/HimalayaUI/frontend/test/test-utils.tsx`:

```tsx
import type { ReactElement, ReactNode } from "react";
import { render, type RenderOptions } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";

export function makeClient(): QueryClient {
  return new QueryClient({
    defaultOptions: { queries: { retry: false, gcTime: 0, staleTime: 0 } },
  });
}

export function renderWithProviders(
  ui: ReactElement,
  options: RenderOptions & { client?: QueryClient } = {},
) {
  const { client = makeClient(), ...rest } = options;
  const Wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return { client, ...render(ui, { wrapper: Wrapper, ...rest }) };
}
```

- [ ] **Step 4: Update `App.tsx` to consume queries**

Replace `packages/HimalayaUI/frontend/src/App.tsx`:

```tsx
import { useState } from "react";
import "./styles.css";
import { useAppState } from "./state";
import { useExperiment, useSamples } from "./queries";
import { Navbar } from "./components/Navbar";
import { Layout } from "./components/Layout";
import { SampleList } from "./components/SampleList";
import { UserModal } from "./components/UserModal";

const EXPERIMENT_ID = 1;

export function App(): JSX.Element {
  const username        = useAppState((s) => s.username);
  const activeSampleId  = useAppState((s) => s.activeSampleId);
  const setUsername     = useAppState((s) => s.setUsername);
  const setActiveSample = useAppState((s) => s.setActiveSample);

  const [modalOpen, setModalOpen] = useState<boolean>(!username);

  const experimentQ = useExperiment(EXPERIMENT_ID);
  const samplesQ    = useSamples(EXPERIMENT_ID);

  const samples      = samplesQ.data ?? [];
  const activeSample = samples.find((s) => s.id === activeSampleId);
  const bootError    = experimentQ.error ?? samplesQ.error;

  const breadcrumb = bootError
    ? `Error: ${(bootError as Error).message}`
    : (experimentQ.data?.name ?? "experiment")
        + (activeSample
          ? ` › ${activeSample.label ?? ""} ${activeSample.name ?? ""}`.trimEnd()
          : "");

  return (
    <>
      <Navbar
        breadcrumb={breadcrumb}
        username={username}
        onUserClick={() => setModalOpen(true)}
      />
      <Layout
        left={
          <SampleList
            samples={samples}
            activeId={activeSampleId}
            onSelect={setActiveSample}
          />
        }
      />
      <UserModal
        open={modalOpen}
        onSelect={(name) => {
          setUsername(name);
          setModalOpen(false);
        }}
        onClose={() => setModalOpen(false)}
      />
    </>
  );
}
```

- [ ] **Step 5: Wrap app in `QueryClientProvider`**

Replace `packages/HimalayaUI/frontend/src/main.tsx`:

```tsx
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { App } from "./App";
import { ErrorBoundary } from "./ErrorBoundary";

const queryClient = new QueryClient({
  defaultOptions: {
    queries: { staleTime: 30_000, retry: 1, refetchOnWindowFocus: false },
  },
});

const root = document.getElementById("app");
if (!root) throw new Error("#app root missing");
createRoot(root).render(
  <StrictMode>
    <ErrorBoundary>
      <QueryClientProvider client={queryClient}>
        <App />
      </QueryClientProvider>
    </ErrorBoundary>
  </StrictMode>,
);
```

- [ ] **Step 6: Update `smoke.test.tsx` to use the provider**

Replace `packages/HimalayaUI/frontend/test/smoke.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, waitFor } from "@testing-library/react";
import { App } from "../src/App";
import { renderWithProviders } from "./test-utils";

function mockFetch(map: Record<string, unknown>): void {
  vi.spyOn(global, "fetch").mockImplementation(async (input) => {
    const url = typeof input === "string" ? input : (input as Request).url;
    const key = Object.keys(map).find((k) => url.endsWith(k));
    if (!key) return new Response("not found", { status: 404 });
    return new Response(JSON.stringify(map[key]), {
      status: 200, headers: { "Content-Type": "application/json" },
    });
  });
}

describe("App smoke", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    localStorage.clear();
    mockFetch({
      "/api/users": [],
      "/api/experiments/1": {
        id: 1, name: "demo", path: "/x", data_dir: "/x/data",
        analysis_dir: "/x/analysis", manifest_path: null,
        created_at: "2026-04-22T00:00:00Z",
      },
      "/api/experiments/1/samples": [
        { id: 10, experiment_id: 1, label: "A1", name: "s1", notes: null, tags: [] },
      ],
    });
  });

  it("renders breadcrumb and sample list after boot", async () => {
    renderWithProviders(<App />);
    await waitFor(() => expect(screen.getByText(/demo/)).toBeInTheDocument());
    expect(screen.getByText("A1")).toBeInTheDocument();
  });
});
```

- [ ] **Step 7: Run full suite**

Run: `npm test`
Expected: all tests pass (smoke test now exercises the provider).

- [ ] **Step 8: Run build**

Run: `npm run build`
Expected: clean `tsc --noEmit`, successful `vite build`.

- [ ] **Step 9: Run E2E to confirm live-app path still works**

Run: `npm run e2e`
Expected: PASS (the existing `/api/*` mocks in `e2e/smoke.spec.ts` still satisfy the Query-driven fetches — same URLs, same responses).

- [ ] **Step 10: Commit**

```bash
cd packages/HimalayaUI/frontend
git add package.json package-lock.json src/queries.ts src/App.tsx src/main.tsx \
        test/smoke.test.tsx test/test-utils.tsx
git commit -m "feat(ui): adopt TanStack Query for server state"
```

---

## Task 4: Tailwind v4 migration

**Motivation:** `styles.css` is a 150-line blob of `:root` vars plus bespoke class names. Tailwind v4 is the React-community default, integrates natively with Vite, and keeps CSS variables as the theme source via `@theme`. This is a mechanical class-for-class swap — no visual redesign.

**Files:**
- Modify: `packages/HimalayaUI/frontend/package.json` (add `tailwindcss`, `@tailwindcss/vite`)
- Modify: `packages/HimalayaUI/frontend/vite.config.ts` (register plugin)
- Modify: `packages/HimalayaUI/frontend/src/styles.css` (Tailwind directives + `@theme`)
- Modify: every component + `App.tsx` (swap `className` literals)
- Delete all bespoke class rules from `styles.css`

### Approach

- Install Tailwind v4 and the Vite plugin.
- Keep existing CSS vars but move them inside `@theme { ... }` so Tailwind exposes them as utilities (`bg-bg`, `text-fg`, `border-border`, etc.).
- Rewrite each component's `className` to use Tailwind utilities. Where a pattern repeats (e.g., pane layout), extract a tiny helper variable in the same file — don't introduce a `classnames` dependency.
- Delete the CSS rules from `styles.css` as each component migrates.
- Keep `styles.css` as the single stylesheet (now just Tailwind + `@theme` + minimal base).

### Steps

- [ ] **Step 1: Install Tailwind v4 + Vite plugin**

```bash
cd packages/HimalayaUI/frontend
npm install -D tailwindcss@^4.0.0 @tailwindcss/vite@^4.0.0
```

Expected: installs cleanly.

- [ ] **Step 2: Register Tailwind in `vite.config.ts`**

Replace `packages/HimalayaUI/frontend/vite.config.ts`:

```ts
import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";
import tailwindcss from "@tailwindcss/vite";

export default defineConfig({
  plugins: [react(), tailwindcss()],
  server: {
    port: 5173,
    proxy: {
      "/api": { target: "http://127.0.0.1:8080", changeOrigin: false },
    },
  },
  build: {
    outDir: "dist",
    emptyOutDir: true,
    sourcemap: true,
  },
});
```

- [ ] **Step 3: Replace `styles.css` with Tailwind + `@theme`**

Replace `packages/HimalayaUI/frontend/src/styles.css`:

```css
@import "tailwindcss";

@theme {
  --color-bg: #1a1a1a;
  --color-bg-elevated: #232323;
  --color-bg-hover: #2a2a2a;
  --color-fg: #e8e6e3;
  --color-fg-muted: #9a9894;
  --color-fg-dim: #666666;
  --color-border: #333333;
  --color-accent: #d97706;
  --color-success: #22c55e;
  --color-warning: #eab308;
  --color-error: #ef4444;
  --radius: 6px;
  --font-sans: -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif;
  --font-mono: ui-monospace, SFMono-Regular, Menlo, monospace;
}

@layer base {
  html, body, #app { height: 100%; }
  body {
    background: var(--color-bg);
    color: var(--color-fg);
    font-family: var(--font-sans);
    font-size: 14px;
    line-height: 1.5;
    -webkit-font-smoothing: antialiased;
  }
}
```

(Everything else is deleted. Utilities replace each former class.)

- [ ] **Step 4: Migrate `Navbar.tsx`**

Replace `packages/HimalayaUI/frontend/src/components/Navbar.tsx`:

```tsx
export interface NavbarProps {
  breadcrumb: string;
  username: string | undefined;
  onUserClick: () => void;
}

export function Navbar({ breadcrumb, username, onUserClick }: NavbarProps): JSX.Element {
  return (
    <header className="grid grid-cols-[auto_1fr_auto] items-center h-11 px-4 border-b border-border bg-bg">
      <div className="flex items-center gap-6">
        <span className="font-semibold tracking-wide">Himalaya</span>
        <nav className="flex gap-4">
          <a className="text-fg border-b border-accent px-1 py-0.5" href="#">Analysis</a>
        </nav>
      </div>
      <div className="text-center text-fg-muted">
        <span className="font-mono text-[13px]">{breadcrumb}</span>
      </div>
      <div className="flex">
        <button
          className="min-w-20 border border-border rounded-md px-2.5 py-1 hover:bg-bg-hover focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent"
          onClick={onUserClick}
        >
          {username ?? "Sign in"}
        </button>
      </div>
    </header>
  );
}
```

- [ ] **Step 5: Migrate `Layout.tsx`**

Replace `packages/HimalayaUI/frontend/src/components/Layout.tsx`:

```tsx
import type { ReactNode } from "react";

export interface LayoutProps {
  left: ReactNode;
  centerTop?: ReactNode;
  centerBottom?: ReactNode;
  rightTop?: ReactNode;
  rightBottom?: ReactNode;
}

const pane = "flex-1 overflow-auto p-4 flex flex-col min-h-0";
const paneBorder = "border-t border-border";
const placeholder = "flex items-center justify-center flex-1 text-fg-dim italic";

function Placeholder({ label }: { label: string }): JSX.Element {
  return <p className={placeholder}>{label}</p>;
}

export function Layout({
  left, centerTop, centerBottom, rightTop, rightBottom,
}: LayoutProps): JSX.Element {
  return (
    <main className="grid grid-cols-[280px_1fr_360px] h-[calc(100vh-44px)] overflow-hidden">
      <aside className="flex flex-col overflow-hidden border-r border-border">{left}</aside>
      <div className="flex flex-col overflow-hidden border-r border-border min-w-0">
        <section className={pane}>{centerTop ?? <Placeholder label="Trace viewer (Plan 4)" />}</section>
        <section className={`${pane} ${paneBorder}`}>{centerBottom ?? <Placeholder label="Properties panel (Plan 6)" />}</section>
      </div>
      <div className="flex flex-col overflow-hidden min-w-0">
        <section className={pane}>{rightTop ?? <Placeholder label="Miller-index plot (Plan 5)" />}</section>
        <section className={`${pane} ${paneBorder}`}>{rightBottom ?? <Placeholder label="Phase panel (Plan 5)" />}</section>
      </div>
    </main>
  );
}
```

- [ ] **Step 6: Migrate `SampleList.tsx`**

Replace `packages/HimalayaUI/frontend/src/components/SampleList.tsx`:

```tsx
import { useState, useMemo } from "react";
import type { Sample } from "../api";

export interface SampleListProps {
  samples: Sample[];
  activeId: number | undefined;
  onSelect: (id: number) => void;
}

type SampleStatus = "unanalyzed" | "candidates" | "confirmed";

const statusDot: Record<SampleStatus, string> = {
  unanalyzed: "bg-fg-dim",
  candidates: "bg-warning",
  confirmed:  "bg-success",
};

function matches(s: Sample, filter: string): boolean {
  if (!filter) return true;
  const n = filter.toLowerCase();
  if (s.label?.toLowerCase().includes(n)) return true;
  if (s.name?.toLowerCase().includes(n))  return true;
  return s.tags.some((t) =>
    t.key.toLowerCase().includes(n) || t.value.toLowerCase().includes(n),
  );
}

export function SampleList({ samples, activeId, onSelect }: SampleListProps): JSX.Element {
  const [filter, setFilter] = useState("");
  const filtered = useMemo(() => samples.filter((s) => matches(s, filter)), [samples, filter]);

  return (
    <div className="flex flex-col h-full">
      <div className="p-2 border-b border-border">
        <input
          className="w-full bg-bg-elevated border border-border rounded-md px-2 py-1 focus:outline focus:outline-1 focus:outline-accent focus:border-accent"
          type="text"
          placeholder="Filter samples…"
          value={filter}
          onChange={(e) => setFilter(e.target.value)}
        />
      </div>
      <ul className="list-none overflow-y-auto flex-1">
        {filtered.map((s) => {
          const status: SampleStatus = "unanalyzed";
          const active = s.id === activeId;
          return (
            <li
              key={s.id}
              className={
                "grid grid-cols-[12px_auto_1fr] items-center gap-2 py-1.5 px-4 cursor-pointer border-l-2 " +
                (active ? "bg-bg-elevated border-accent" : "border-transparent hover:bg-bg-hover")
              }
              data-sample-id={s.id}
              data-status={status}
              onClick={() => onSelect(s.id)}
            >
              <span className={`w-2 h-2 rounded-full ${statusDot[status]}`} />
              <span className="font-medium">{s.label ?? ""}</span>
              <span className="text-fg-muted">{s.name ?? ""}</span>
            </li>
          );
        })}
      </ul>
    </div>
  );
}
```

- [ ] **Step 7: Migrate `UserModal.tsx`**

Replace `packages/HimalayaUI/frontend/src/components/UserModal.tsx`:

```tsx
import { useEffect, useState } from "react";
import * as api from "../api";
import type { User } from "../api";

export interface UserModalProps {
  open: boolean;
  onSelect: (username: string) => void;
  onClose: () => void;
}

export function UserModal({ open, onSelect, onClose }: UserModalProps): JSX.Element | null {
  const [users, setUsers]     = useState<User[]>([]);
  const [selection, setSel]   = useState<string>("__new__");
  const [newName, setNewName] = useState("");
  const [error, setError]     = useState<string | null>(null);

  useEffect(() => {
    if (!open) return;
    setError(null);
    void (async () => {
      try {
        const list = await api.listUsers();
        setUsers(list);
        setSel(list.length === 0 ? "__new__" : list[0]!.username);
      } catch (e) {
        setError(`Failed to load users: ${(e as Error).message}`);
      }
    })();
  }, [open]);

  useEffect(() => {
    function onKey(e: KeyboardEvent): void {
      if (e.key === "Escape" && open) onClose();
    }
    window.addEventListener("keydown", onKey);
    return () => { window.removeEventListener("keydown", onKey); };
  }, [open, onClose]);

  if (!open) return null;

  async function submit(): Promise<void> {
    setError(null);
    const username = selection === "__new__" ? newName.trim() : selection;
    if (!username) {
      setError("Username required");
      return;
    }
    try {
      await api.createUser(username, { username });
      onSelect(username);
    } catch (e) {
      setError(`Failed: ${(e as Error).message}`);
    }
  }

  const inputClass =
    "w-full bg-bg-elevated border border-border rounded-md px-2 py-1 focus:outline focus:outline-1 focus:outline-accent focus:border-accent";

  return (
    <div
      className="fixed inset-0 bg-black/60 flex items-center justify-center z-50"
      role="presentation"
      onClick={(e) => { if (e.target === e.currentTarget) onClose(); }}
    >
      <div className="bg-bg-elevated border border-border rounded-md p-6 min-w-[360px] max-w-[480px] flex flex-col gap-4" role="dialog" aria-modal="true">
        <h2 className="text-base font-semibold">Who are you?</h2>
        <p className="text-fg-muted">
          Your name is stored with every change you make so others can see what you've done.
        </p>
        <select
          className={inputClass}
          value={selection}
          onChange={(e) => setSel(e.target.value)}
        >
          {users.map((u) => (
            <option key={u.id} value={u.username}>{u.username}</option>
          ))}
          <option value="__new__">+ New user…</option>
        </select>
        {selection === "__new__" && (
          <input
            className={inputClass}
            type="text"
            placeholder="Enter username"
            value={newName}
            onChange={(e) => setNewName(e.target.value)}
            autoFocus
          />
        )}
        {error && <p className="text-error text-[13px]">{error}</p>}
        <div className="flex justify-end gap-2">
          <button
            className="bg-accent border border-accent text-white rounded-md px-2.5 py-1 hover:brightness-110"
            onClick={submit}
          >
            Continue
          </button>
        </div>
      </div>
    </div>
  );
}
```

- [ ] **Step 8: Run tests**

Run: `npm test`
Expected: PASS — the tests assert on text content and `data-*` attributes, not bespoke class names, so they still hold. If any existing test asserted on an old class literal (check `SampleList.test.tsx`), update the assertion to reference stable attributes (`data-sample-id`, `data-status`) instead.

- [ ] **Step 9: Run build**

Run: `npm run build`
Expected: clean `tsc --noEmit`; Tailwind emits utilities into the final bundle; `vite build` succeeds.

- [ ] **Step 10: Visual smoke (dev server)**

Run in one terminal: `npm run dev`
Open `http://localhost:5173`.
Expected: app renders with the same dark theme, breadcrumb, left sample list, three placeholder panes, and working user modal. Pixel-perfect parity is not required; no console errors, no missing layout.

Stop the dev server (Ctrl-C).

- [ ] **Step 11: Run E2E to confirm selectors still match**

Run: `npm run e2e`
Expected: PASS — E2E uses stable selectors (text content, `data-sample-id`); class changes don't break it.

- [ ] **Step 12: Commit**

```bash
cd packages/HimalayaUI/frontend
git add package.json package-lock.json vite.config.ts src/styles.css \
        src/App.tsx src/components/Navbar.tsx src/components/Layout.tsx \
        src/components/SampleList.tsx src/components/UserModal.tsx
git commit -m "refactor(ui): migrate styles to Tailwind v4"
```

---

## Final verification

After all four tasks land, from `packages/HimalayaUI/frontend/`:

```bash
npm test && npm run build && npm run e2e
```

Expected: all unit tests pass, build is clean, E2E passes. From repo root:

```bash
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
```

Expected: backend still passes (untouched by this plan, but worth a sanity run).
