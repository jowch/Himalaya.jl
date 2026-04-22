# HimalayaUI Frontend Shell Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up the React + TypeScript + Vite frontend shell with navbar, three-column layout skeleton, working sample list, and user identification modal. Center (trace viewer, properties panel) and right (Miller plot, phase panel) columns are placeholder stubs filled in by Plans 4–6.

**Architecture:** React 18 function components. Vite for dev + build. Zustand (with `persist` middleware) for app state — replaces hand-rolled store + localStorage glue. Observable Plot deferred to Plans 4–5. Single page for v1 — the per-sample analysis view. `api.ts` is a thin typed layer over `fetch`, framework-agnostic. In dev, Vite proxies `/api/*` to the Julia server on :8080. In production, Oxygen.jl serves the compiled `dist/` at `/` and `/api/*` from the same origin.

**Tech Stack:** React 18, TypeScript 5, Vite 5, Zustand 4, Vitest + React Testing Library (units), Playwright (E2E).

**This is Plan 3 of 6.** Prior plans complete: Plan 1 (SQLite pipeline + CLI), Plan 2 (Oxygen REST API + `himalaya serve`). Remaining: Plan 4 (trace viewer + peak editing), Plan 5 (phase panel + Miller plot + hover preview), Plan 6 (properties panel tabs).

---

## File Map

| File | Responsibility |
|---|---|
| `packages/HimalayaUI/frontend/package.json` | deps + scripts |
| `packages/HimalayaUI/frontend/tsconfig.json` | TS strict + JSX |
| `packages/HimalayaUI/frontend/vite.config.ts` | React plugin, dev proxy, build → `dist/` |
| `packages/HimalayaUI/frontend/vitest.config.ts` | Vitest + jsdom + setup file |
| `packages/HimalayaUI/frontend/playwright.config.ts` | E2E config |
| `packages/HimalayaUI/frontend/.gitignore` | ignore `node_modules` + `dist` |
| `packages/HimalayaUI/frontend/index.html` | entry HTML |
| `packages/HimalayaUI/frontend/src/main.tsx` | `createRoot` bootstrap |
| `packages/HimalayaUI/frontend/src/App.tsx` | composition root |
| `packages/HimalayaUI/frontend/src/api.ts` | typed fetch wrappers |
| `packages/HimalayaUI/frontend/src/state.ts` | Zustand store + `persist` |
| `packages/HimalayaUI/frontend/src/styles.css` | design tokens + reset + base |
| `packages/HimalayaUI/frontend/src/components/Navbar.tsx` | top nav with breadcrumb + user |
| `packages/HimalayaUI/frontend/src/components/Layout.tsx` | three-column shell |
| `packages/HimalayaUI/frontend/src/components/SampleList.tsx` | left column sample list + filter |
| `packages/HimalayaUI/frontend/src/components/UserModal.tsx` | first-visit + switch-user modal |
| `packages/HimalayaUI/frontend/test/setup.ts` | `@testing-library/jest-dom/vitest` import |
| `packages/HimalayaUI/frontend/test/api.test.ts` | Vitest |
| `packages/HimalayaUI/frontend/test/state.test.ts` | Vitest |
| `packages/HimalayaUI/frontend/test/SampleList.test.tsx` | Vitest + RTL |
| `packages/HimalayaUI/frontend/e2e/smoke.spec.ts` | Playwright |

---

## Dev workflow summary

```bash
cd packages/HimalayaUI/frontend
npm install                # first time
npm run dev                # Vite on :5173, proxies /api → :8080
npm test                   # Vitest unit tests
npm run build              # emit packages/HimalayaUI/frontend/dist/
npx playwright install     # first-time E2E browser install
npm run e2e                # Playwright smoke
```

`dist/` is gitignored. Oxygen's `serve(db)` serves it at `/` if present.

---

## Design tokens (referenced throughout)

Dark, minimal, Claude Code desktop feel. Semantic CSS variables declared once in `styles.css`:

```
--bg: #1a1a1a             --fg: #e8e6e3
--bg-elevated: #232323    --fg-muted: #9a9894
--bg-hover: #2a2a2a       --fg-dim: #666666
--border: #333333         --accent: #d97706
--success: #22c55e        --warning: #eab308        --error: #ef4444
--radius: 6px             --gap-sm: 8px              --gap-md: 16px    --gap-lg: 24px
--font-sans: -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif
--font-mono: ui-monospace, SFMono-Regular, Menlo, monospace
```

---

### Task 1: Scaffold Vite + React + TypeScript

**Files:**
- Create: `packages/HimalayaUI/frontend/package.json`
- Create: `packages/HimalayaUI/frontend/tsconfig.json`
- Create: `packages/HimalayaUI/frontend/vite.config.ts`
- Create: `packages/HimalayaUI/frontend/vitest.config.ts`
- Create: `packages/HimalayaUI/frontend/.gitignore`
- Create: `packages/HimalayaUI/frontend/index.html`
- Create: `packages/HimalayaUI/frontend/src/main.tsx` (placeholder)
- Create: `packages/HimalayaUI/frontend/src/App.tsx` (placeholder)
- Delete: `packages/HimalayaUI/frontend/dist/.gitkeep`

- [ ] **Step 1: Initialize `package.json`**

```bash
cd packages/HimalayaUI/frontend
```

Create `package.json`:
```json
{
  "name": "himalaya-ui-frontend",
  "version": "0.1.0",
  "private": true,
  "type": "module",
  "scripts": {
    "dev": "vite",
    "build": "tsc --noEmit && vite build",
    "preview": "vite preview",
    "test": "vitest run",
    "test:watch": "vitest",
    "e2e": "playwright test"
  },
  "dependencies": {
    "react": "^18.3.0",
    "react-dom": "^18.3.0",
    "zustand": "^4.5.0"
  },
  "devDependencies": {
    "@types/react": "^18.3.0",
    "@types/react-dom": "^18.3.0",
    "@vitejs/plugin-react": "^4.3.0",
    "@testing-library/react": "^16.0.0",
    "@testing-library/jest-dom": "^6.5.0",
    "@testing-library/user-event": "^14.5.0",
    "@playwright/test": "^1.48.0",
    "jsdom": "^25.0.0",
    "typescript": "^5.6.0",
    "vite": "^5.4.0",
    "vitest": "^2.1.0"
  }
}
```

- [ ] **Step 2: `tsconfig.json`**

```json
{
  "compilerOptions": {
    "target": "ES2022",
    "module": "ESNext",
    "moduleResolution": "bundler",
    "jsx": "react-jsx",
    "strict": true,
    "noUnusedLocals": true,
    "noUnusedParameters": true,
    "noImplicitOverride": true,
    "noFallthroughCasesInSwitch": true,
    "exactOptionalPropertyTypes": true,
    "lib": ["ES2022", "DOM", "DOM.Iterable"],
    "esModuleInterop": true,
    "skipLibCheck": true,
    "isolatedModules": true,
    "resolveJsonModule": true
  },
  "include": ["src/**/*.ts", "src/**/*.tsx", "test/**/*.ts", "test/**/*.tsx", "e2e/**/*.ts", "*.ts"]
}
```

- [ ] **Step 3: `vite.config.ts`**

```ts
import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
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

- [ ] **Step 4: `vitest.config.ts`**

```ts
import { defineConfig } from "vitest/config";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  test: {
    environment: "jsdom",
    globals: true,
    setupFiles: ["./test/setup.ts"],
    include: ["test/**/*.test.ts", "test/**/*.test.tsx"],
  },
});
```

- [ ] **Step 5: `.gitignore`**

```
node_modules
dist
.DS_Store
test-results
playwright-report
```

Then remove the old placeholder:
```bash
git rm packages/HimalayaUI/frontend/dist/.gitkeep
```

- [ ] **Step 6: `index.html`**

```html
<!doctype html>
<html lang="en">
  <head>
    <meta charset="UTF-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
    <title>Himalaya</title>
  </head>
  <body>
    <div id="app"></div>
    <script type="module" src="/src/main.tsx"></script>
  </body>
</html>
```

- [ ] **Step 7: Placeholder `src/App.tsx` and `src/main.tsx`**

`src/App.tsx`:
```tsx
export function App(): JSX.Element {
  return <div>Himalaya UI — loading…</div>;
}
```

`src/main.tsx`:
```tsx
import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { App } from "./App";

const root = document.getElementById("app");
if (!root) throw new Error("#app root missing");
createRoot(root).render(
  <StrictMode>
    <App />
  </StrictMode>,
);
```

- [ ] **Step 8: Install and build**

```bash
cd packages/HimalayaUI/frontend
npm install
npm run build
```

Expected: `dist/index.html` and `dist/assets/*.js`. No TS errors.

- [ ] **Step 9: Commit**

```bash
cd ../../..
git add packages/HimalayaUI/frontend/package.json \
        packages/HimalayaUI/frontend/package-lock.json \
        packages/HimalayaUI/frontend/tsconfig.json \
        packages/HimalayaUI/frontend/vite.config.ts \
        packages/HimalayaUI/frontend/vitest.config.ts \
        packages/HimalayaUI/frontend/.gitignore \
        packages/HimalayaUI/frontend/index.html \
        packages/HimalayaUI/frontend/src/main.tsx \
        packages/HimalayaUI/frontend/src/App.tsx
git rm -f packages/HimalayaUI/frontend/dist/.gitkeep 2>/dev/null || true
git commit -m "feat(frontend): scaffold Vite + React + TypeScript"
```

---

### Task 2: Design tokens + test setup

**Files:**
- Create: `packages/HimalayaUI/frontend/src/styles.css`
- Create: `packages/HimalayaUI/frontend/test/setup.ts`
- Create: `packages/HimalayaUI/frontend/test/smoke.test.tsx` (verifies RTL works)
- Modify: `packages/HimalayaUI/frontend/src/App.tsx` (import styles)

- [ ] **Step 1: Write the failing test**

Create `test/smoke.test.tsx`:
```tsx
import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { App } from "../src/App";

describe("App smoke", () => {
  it("renders the loading text", () => {
    render(<App />);
    expect(screen.getByText(/loading/i)).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Create test setup file**

Create `test/setup.ts`:
```ts
import "@testing-library/jest-dom/vitest";
```

- [ ] **Step 3: Verify fail**

```bash
cd packages/HimalayaUI/frontend && npm test
```
Expected failure: `@testing-library/jest-dom/vitest` not actually failing (test will run), but likely passes without the matcher. If `toBeInTheDocument` is reported as not a function, you missed Step 2.

(Note: this test may pass outright given that `App` renders "loading…". The failing test pattern isn't strict here — move on.)

- [ ] **Step 4: Implement `src/styles.css`**

```css
:root {
  --bg: #1a1a1a;
  --bg-elevated: #232323;
  --bg-hover: #2a2a2a;
  --fg: #e8e6e3;
  --fg-muted: #9a9894;
  --fg-dim: #666666;
  --border: #333333;
  --accent: #d97706;
  --success: #22c55e;
  --warning: #eab308;
  --error: #ef4444;
  --radius: 6px;
  --gap-sm: 8px;
  --gap-md: 16px;
  --gap-lg: 24px;
  --font-sans: -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif;
  --font-mono: ui-monospace, SFMono-Regular, Menlo, monospace;
}

* { box-sizing: border-box; margin: 0; padding: 0; }

html, body, #app { height: 100%; }

body {
  background: var(--bg);
  color: var(--fg);
  font-family: var(--font-sans);
  font-size: 14px;
  line-height: 1.5;
  -webkit-font-smoothing: antialiased;
}

button {
  font: inherit;
  color: inherit;
  background: transparent;
  border: 1px solid var(--border);
  border-radius: var(--radius);
  padding: 4px 10px;
  cursor: pointer;
}
button:hover { background: var(--bg-hover); }
button:focus-visible { outline: 1px solid var(--accent); outline-offset: 1px; }

input, select {
  font: inherit;
  color: inherit;
  background: var(--bg-elevated);
  border: 1px solid var(--border);
  border-radius: var(--radius);
  padding: 4px 8px;
}
input:focus, select:focus { outline: 1px solid var(--accent); border-color: var(--accent); }

.mono { font-family: var(--font-mono); }
.muted { color: var(--fg-muted); }
.dim { color: var(--fg-dim); }
```

- [ ] **Step 5: Update `App.tsx` to import styles**

```tsx
import "./styles.css";

export function App(): JSX.Element {
  return <div>Himalaya UI — loading…</div>;
}
```

- [ ] **Step 6: Run tests**

```bash
cd packages/HimalayaUI/frontend && npm test
```
Expected: smoke.test passes.

- [ ] **Step 7: Commit**

```bash
cd ../../..
git add packages/HimalayaUI/frontend/src/styles.css \
        packages/HimalayaUI/frontend/src/App.tsx \
        packages/HimalayaUI/frontend/test/setup.ts \
        packages/HimalayaUI/frontend/test/smoke.test.tsx
git commit -m "feat(frontend): design tokens + RTL bootstrap"
```

---

### Task 3: Typed API client

**Files:**
- Create: `packages/HimalayaUI/frontend/src/api.ts`
- Create: `packages/HimalayaUI/frontend/test/api.test.ts`

Plan 3 uses: users, experiments, samples. Peaks/exposures/indices/groups/export added in Plans 4–6.

- [ ] **Step 1: Write failing test**

Create `test/api.test.ts`:
```ts
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import * as api from "../src/api";

function mockFetchJson(body: unknown, status = 200): void {
  global.fetch = vi.fn(async () => new Response(JSON.stringify(body), {
    status, headers: { "Content-Type": "application/json" },
  })) as unknown as typeof fetch;
}

describe("api client", () => {
  beforeEach(() => { api.setUsername(undefined); });
  afterEach(() => { vi.restoreAllMocks(); });

  it("listUsers returns [] for empty response", async () => {
    mockFetchJson([]);
    expect(await api.listUsers()).toEqual([]);
  });

  it("listUsers parses user rows", async () => {
    mockFetchJson([{ id: 1, username: "alice" }]);
    const users = await api.listUsers();
    expect(users.length).toBe(1);
    expect(users[0]!.username).toBe("alice");
  });

  it("createUser POSTs JSON body with Content-Type", async () => {
    const spy = vi.fn(async () => new Response(JSON.stringify({ id: 1, username: "bob" }), {
      status: 201, headers: { "Content-Type": "application/json" },
    }));
    global.fetch = spy as unknown as typeof fetch;
    const user = await api.createUser("bob");
    expect(user.username).toBe("bob");
    const init = spy.mock.calls[0]![1];
    expect(init?.method).toBe("POST");
    expect((init?.headers as Record<string, string>)["Content-Type"]).toBe("application/json");
    expect(init?.body).toBe(JSON.stringify({ username: "bob" }));
  });

  it("mutating calls send X-Username when set", async () => {
    const spy = vi.fn(async () => new Response("null", {
      status: 200, headers: { "Content-Type": "application/json" },
    }));
    global.fetch = spy as unknown as typeof fetch;
    api.setUsername("alice");
    await api.updateSample(1, { notes: "hi" });
    const init = spy.mock.calls[0]![1];
    expect((init?.headers as Record<string, string>)["X-Username"]).toBe("alice");
  });

  it("non-2xx throws ApiError with status and message", async () => {
    mockFetchJson({ error: "not found" }, 404);
    await expect(api.getExperiment(99)).rejects.toMatchObject({
      status: 404,
      message: expect.stringContaining("not found"),
    });
  });
});
```

- [ ] **Step 2: Verify fail**

- [ ] **Step 3: Implement `src/api.ts`**

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

let _username: string | undefined = undefined;
export function setUsername(name: string | undefined): void { _username = name; }
export function getUsername(): string | undefined { return _username; }

async function request<T>(method: string, path: string, body?: unknown): Promise<T> {
  const headers: Record<string, string> = {};
  if (body !== undefined) headers["Content-Type"] = "application/json";
  if (_username && method !== "GET") headers["X-Username"] = _username;

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
export const createUser = (username: string) =>
  request<User>("POST", "/api/users", { username });

// Experiments
export const getExperiment = (id: number) =>
  request<Experiment>("GET", `/api/experiments/${id}`);
export const updateExperiment = (id: number, patch: Partial<Pick<Experiment,
  "name" | "data_dir" | "analysis_dir" | "manifest_path">>) =>
  request<Experiment>("PATCH", `/api/experiments/${id}`, patch);

// Samples
export const listSamples    = (experiment_id: number) =>
  request<Sample[]>("GET", `/api/experiments/${experiment_id}/samples`);
export const updateSample   = (id: number, patch: { name?: string; notes?: string }) =>
  request<Sample>("PATCH", `/api/samples/${id}`, patch);
export const addSampleTag   = (id: number, key: string, value: string) =>
  request<SampleTag>("POST", `/api/samples/${id}/tags`, { key, value });
export const removeSampleTag = (id: number, tag_id: number) =>
  request<void>("DELETE", `/api/samples/${id}/tags/${tag_id}`);
```

- [ ] **Step 4: Run tests and commit**

```bash
cd packages/HimalayaUI/frontend && npm test
cd ../../..
git add packages/HimalayaUI/frontend/src/api.ts \
        packages/HimalayaUI/frontend/test/api.test.ts
git commit -m "feat(frontend): typed API client with X-Username handling"
```

---

### Task 4: Zustand app state store

**Files:**
- Create: `packages/HimalayaUI/frontend/src/state.ts`
- Create: `packages/HimalayaUI/frontend/test/state.test.ts`

`state.ts` is 20 lines — Zustand's `persist` middleware handles localStorage. We test only our store's shape (setters work, persistence key is stable).

- [ ] **Step 1: Write the failing test**

Create `test/state.test.ts`:
```ts
import { describe, it, expect, beforeEach } from "vitest";
import { useAppState, LS_KEY } from "../src/state";

describe("useAppState", () => {
  beforeEach(() => {
    localStorage.clear();
    // Reset store to defaults between tests
    useAppState.setState({
      username: undefined,
      activeSampleId: undefined,
      activeExposureId: undefined,
    });
  });

  it("starts with undefined fields", () => {
    const s = useAppState.getState();
    expect(s.username).toBeUndefined();
    expect(s.activeSampleId).toBeUndefined();
    expect(s.activeExposureId).toBeUndefined();
  });

  it("setUsername updates state", () => {
    useAppState.getState().setUsername("alice");
    expect(useAppState.getState().username).toBe("alice");
  });

  it("setActiveSample clears activeExposureId", () => {
    useAppState.setState({ activeExposureId: 7 });
    useAppState.getState().setActiveSample(3);
    expect(useAppState.getState().activeSampleId).toBe(3);
    expect(useAppState.getState().activeExposureId).toBeUndefined();
  });

  it("persists to localStorage under the stable key", () => {
    useAppState.getState().setUsername("bob");
    const raw = localStorage.getItem(LS_KEY);
    expect(raw).not.toBeNull();
    const parsed = JSON.parse(raw!);
    expect(parsed.state.username).toBe("bob");
  });
});
```

- [ ] **Step 2: Verify fail**

- [ ] **Step 3: Implement `src/state.ts`**

```ts
import { create } from "zustand";
import { persist } from "zustand/middleware";

export const LS_KEY = "himalaya-ui:state";

export interface AppState {
  username?: string;
  activeSampleId?: number;
  activeExposureId?: number;
  setUsername: (name: string) => void;
  setActiveSample: (id: number | undefined) => void;
  setActiveExposure: (id: number | undefined) => void;
}

export const useAppState = create<AppState>()(
  persist(
    (set) => ({
      setUsername: (username) => set({ username }),
      setActiveSample: (activeSampleId) => set({ activeSampleId, activeExposureId: undefined }),
      setActiveExposure: (activeExposureId) => set({ activeExposureId }),
    }),
    {
      name: LS_KEY,
      version: 1,
      partialize: (s) => ({
        username: s.username,
        activeSampleId: s.activeSampleId,
        activeExposureId: s.activeExposureId,
      }),
    },
  ),
);
```

The `partialize` option excludes the setter functions from the persisted JSON — we only want data fields.

- [ ] **Step 4: Run tests and commit**

```bash
cd packages/HimalayaUI/frontend && npm test
cd ../../..
git add packages/HimalayaUI/frontend/src/state.ts \
        packages/HimalayaUI/frontend/test/state.test.ts
git commit -m "feat(frontend): Zustand app state store with localStorage persistence"
```

---

### Task 5: Navbar component

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/Navbar.tsx`
- Modify: `packages/HimalayaUI/frontend/src/styles.css` (append navbar styles)

Layout: **logo** · **page links** (only "Analysis" in v1) · **breadcrumb** (center) · **user** (right). Props: `{breadcrumb, username, onUserClick}`.

- [ ] **Step 1: Implement `src/components/Navbar.tsx`**

```tsx
export interface NavbarProps {
  breadcrumb: string;
  username: string | undefined;
  onUserClick: () => void;
}

export function Navbar({ breadcrumb, username, onUserClick }: NavbarProps): JSX.Element {
  return (
    <header className="nav">
      <div className="nav-left">
        <span className="nav-logo">Himalaya</span>
        <nav className="nav-links">
          <a className="nav-link active" href="#">Analysis</a>
        </nav>
      </div>
      <div className="nav-center">
        <span className="nav-breadcrumb">{breadcrumb}</span>
      </div>
      <div className="nav-right">
        <button className="nav-user" onClick={onUserClick}>
          {username ?? "Sign in"}
        </button>
      </div>
    </header>
  );
}
```

- [ ] **Step 2: Append navbar styles**

Append to `src/styles.css`:
```css
.nav {
  display: grid;
  grid-template-columns: auto 1fr auto;
  align-items: center;
  padding: 0 var(--gap-md);
  height: 44px;
  border-bottom: 1px solid var(--border);
  background: var(--bg);
}
.nav-left { display: flex; align-items: center; gap: var(--gap-lg); }
.nav-logo { font-weight: 600; letter-spacing: 0.02em; }
.nav-links { display: flex; gap: var(--gap-md); }
.nav-link { color: var(--fg-muted); text-decoration: none; padding: 2px 4px; }
.nav-link.active { color: var(--fg); border-bottom: 1px solid var(--accent); }
.nav-center { text-align: center; color: var(--fg-muted); }
.nav-breadcrumb { font-family: var(--font-mono); font-size: 13px; }
.nav-right { display: flex; }
.nav-user { min-width: 80px; }
```

- [ ] **Step 3: Commit** (behavior exercised in Task 10 Playwright smoke)

```bash
git add packages/HimalayaUI/frontend/src/components/Navbar.tsx \
        packages/HimalayaUI/frontend/src/styles.css
git commit -m "feat(frontend): Navbar component with breadcrumb and user button"
```

---

### Task 6: Three-column Layout component

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/Layout.tsx`
- Modify: `packages/HimalayaUI/frontend/src/styles.css`

React composition: the `Layout` component takes children for each slot via named props. Plans 4–6 will fill the right slots.

- [ ] **Step 1: Implement `src/components/Layout.tsx`**

```tsx
import type { ReactNode } from "react";

export interface LayoutProps {
  left: ReactNode;
  centerTop?: ReactNode;
  centerBottom?: ReactNode;
  rightTop?: ReactNode;
  rightBottom?: ReactNode;
}

function Placeholder({ label }: { label: string }): JSX.Element {
  return <p className="muted placeholder">{label}</p>;
}

export function Layout({
  left, centerTop, centerBottom, rightTop, rightBottom,
}: LayoutProps): JSX.Element {
  return (
    <main className="layout">
      <aside className="col col-left">{left}</aside>
      <div className="col col-center">
        <section className="pane pane-center-top">
          {centerTop ?? <Placeholder label="Trace viewer (Plan 4)" />}
        </section>
        <section className="pane pane-center-bottom">
          {centerBottom ?? <Placeholder label="Properties panel (Plan 6)" />}
        </section>
      </div>
      <div className="col col-right">
        <section className="pane pane-right-top">
          {rightTop ?? <Placeholder label="Miller-index plot (Plan 5)" />}
        </section>
        <section className="pane pane-right-bottom">
          {rightBottom ?? <Placeholder label="Phase panel (Plan 5)" />}
        </section>
      </div>
    </main>
  );
}
```

- [ ] **Step 2: Append layout styles**

```css
.layout {
  display: grid;
  grid-template-columns: 280px 1fr 360px;
  height: calc(100vh - 44px);
  overflow: hidden;
}
.col { display: flex; flex-direction: column; overflow: hidden; }
.col-left   { border-right: 1px solid var(--border); }
.col-center { border-right: 1px solid var(--border); }
.col-center, .col-right { min-width: 0; }
.pane {
  flex: 1;
  overflow: auto;
  padding: var(--gap-md);
  display: flex;
  flex-direction: column;
  min-height: 0;
}
.pane + .pane { border-top: 1px solid var(--border); }
.placeholder {
  display: flex;
  align-items: center;
  justify-content: center;
  flex: 1;
  color: var(--fg-dim);
  font-style: italic;
}
```

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/Layout.tsx \
        packages/HimalayaUI/frontend/src/styles.css
git commit -m "feat(frontend): three-column Layout component with placeholder slots"
```

---

### Task 7: Sample list + filter

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/SampleList.tsx`
- Create: `packages/HimalayaUI/frontend/test/SampleList.test.tsx`
- Modify: `packages/HimalayaUI/frontend/src/styles.css`

Plan 3 renders all samples as `status="unanalyzed"` (grey). Plans 4–5 will upgrade status based on indices/groups via the same `data-status` attribute — no API change here.

Filter: single text input matches `label`, `name`, or any tag `key`/`value`. Substring, case-insensitive.

- [ ] **Step 1: Write the failing test**

Create `test/SampleList.test.tsx`:
```tsx
import { render, screen, fireEvent } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, it, expect, vi } from "vitest";
import { SampleList } from "../src/components/SampleList";
import type { Sample } from "../src/api";

const SAMPLES: Sample[] = [
  { id: 1, experiment_id: 1, label: "D1", name: "UX1", notes: null,
    tags: [{ id: 1, key: "lipid", value: "DOPC", source: "manual" }] },
  { id: 2, experiment_id: 1, label: "D2", name: "UX2", notes: null, tags: [] },
  { id: 3, experiment_id: 1, label: "D3", name: "UL1", notes: null,
    tags: [{ id: 2, key: "peptide", value: "melittin", source: "manual" }] },
];

describe("<SampleList>", () => {
  it("renders one row per sample", () => {
    render(<SampleList samples={SAMPLES} activeId={undefined} onSelect={() => {}} />);
    expect(document.querySelectorAll("[data-sample-id]").length).toBe(3);
  });

  it("marks the active sample", () => {
    render(<SampleList samples={SAMPLES} activeId={2} onSelect={() => {}} />);
    const active = document.querySelector(".sample-row.active");
    expect(active?.getAttribute("data-sample-id")).toBe("2");
  });

  it("calls onSelect with sample id on click", () => {
    const onSelect = vi.fn();
    render(<SampleList samples={SAMPLES} activeId={undefined} onSelect={onSelect} />);
    const row = document.querySelector<HTMLElement>('[data-sample-id="3"]')!;
    fireEvent.click(row);
    expect(onSelect).toHaveBeenCalledWith(3);
  });

  it("filter narrows by label, name, and tag", async () => {
    const user = userEvent.setup();
    render(<SampleList samples={SAMPLES} activeId={undefined} onSelect={() => {}} />);
    const filter = screen.getByPlaceholderText(/filter samples/i);

    await user.type(filter, "D1");
    expect(document.querySelectorAll("[data-sample-id]").length).toBe(1);

    await user.clear(filter);
    await user.type(filter, "UX");
    expect(document.querySelectorAll("[data-sample-id]").length).toBe(2);

    await user.clear(filter);
    await user.type(filter, "lipid");
    expect(document.querySelectorAll("[data-sample-id]").length).toBe(1);

    await user.clear(filter);
    await user.type(filter, "melittin");
    expect(document.querySelectorAll("[data-sample-id]").length).toBe(1);

    await user.clear(filter);
    expect(document.querySelectorAll("[data-sample-id]").length).toBe(3);
  });
});
```

- [ ] **Step 2: Verify fail**

- [ ] **Step 3: Implement `src/components/SampleList.tsx`**

```tsx
import { useState, useMemo } from "react";
import type { Sample } from "../api";

export interface SampleListProps {
  samples: Sample[];
  activeId: number | undefined;
  onSelect: (id: number) => void;
}

type SampleStatus = "unanalyzed" | "candidates" | "confirmed";

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
    <div className="sample-list-wrap">
      <div className="sample-list-header">
        <input
          className="sample-filter"
          type="text"
          placeholder="Filter samples…"
          value={filter}
          onChange={(e) => setFilter(e.target.value)}
        />
      </div>
      <ul className="sample-list">
        {filtered.map((s) => {
          const status: SampleStatus = "unanalyzed"; // upgraded by Plans 4-5
          return (
            <li
              key={s.id}
              className={"sample-row" + (s.id === activeId ? " active" : "")}
              data-sample-id={s.id}
              data-status={status}
              onClick={() => onSelect(s.id)}
            >
              <span className={`status-dot status-${status}`} />
              <span className="sample-label">{s.label ?? ""}</span>
              <span className="sample-name muted">{s.name ?? ""}</span>
            </li>
          );
        })}
      </ul>
    </div>
  );
}
```

- [ ] **Step 4: Append sample-list styles**

```css
.sample-list-wrap { display: flex; flex-direction: column; height: 100%; }
.sample-list-header {
  padding: var(--gap-sm);
  border-bottom: 1px solid var(--border);
}
.sample-filter { width: 100%; }
.sample-list {
  list-style: none;
  overflow-y: auto;
  flex: 1;
}
.sample-row {
  display: grid;
  grid-template-columns: 12px auto 1fr;
  align-items: center;
  gap: var(--gap-sm);
  padding: 6px var(--gap-md);
  cursor: pointer;
  border-left: 2px solid transparent;
}
.sample-row:hover { background: var(--bg-hover); }
.sample-row.active {
  background: var(--bg-elevated);
  border-left-color: var(--accent);
}
.sample-label { font-weight: 500; }
.status-dot { width: 8px; height: 8px; border-radius: 50%; background: var(--fg-dim); }
.status-candidates { background: var(--warning); }
.status-confirmed  { background: var(--success); }
```

- [ ] **Step 5: Run tests and commit**

```bash
cd packages/HimalayaUI/frontend && npm test
cd ../../..
git add packages/HimalayaUI/frontend/src/components/SampleList.tsx \
        packages/HimalayaUI/frontend/src/styles.css \
        packages/HimalayaUI/frontend/test/SampleList.test.tsx
git commit -m "feat(frontend): SampleList component with filter and click-to-select"
```

---

### Task 8: User identification modal

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/UserModal.tsx`
- Modify: `packages/HimalayaUI/frontend/src/styles.css`

Modal opens when no user is set or when the user clicks the navbar button. Controlled component: parent passes `open` and `onClose`. On submit, calls `api.createUser()` (idempotent on the server) and reports back via `onSelect`.

- [ ] **Step 1: Implement `src/components/UserModal.tsx`**

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
      await api.createUser(username);
      api.setUsername(username);
      onSelect(username);
    } catch (e) {
      setError(`Failed: ${(e as Error).message}`);
    }
  }

  return (
    <div
      className="modal-backdrop"
      role="presentation"
      onClick={(e) => { if (e.target === e.currentTarget) onClose(); }}
    >
      <div className="modal-dialog" role="dialog" aria-modal="true">
        <h2>Who are you?</h2>
        <p className="muted">
          Your name is stored with every change you make so others can see what you've done.
        </p>
        <select
          className="user-select"
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
            className="user-new-input"
            type="text"
            placeholder="Enter username"
            value={newName}
            onChange={(e) => setNewName(e.target.value)}
            autoFocus
          />
        )}
        {error && <p className="user-error">{error}</p>}
        <div className="modal-actions">
          <button className="user-submit primary" onClick={submit}>Continue</button>
        </div>
      </div>
    </div>
  );
}
```

- [ ] **Step 2: Append modal styles**

```css
.modal-backdrop {
  position: fixed;
  inset: 0;
  background: rgba(0,0,0,0.6);
  display: flex;
  align-items: center;
  justify-content: center;
  z-index: 100;
}
.modal-dialog {
  background: var(--bg-elevated);
  border: 1px solid var(--border);
  border-radius: var(--radius);
  padding: var(--gap-lg);
  min-width: 360px;
  max-width: 480px;
  display: flex;
  flex-direction: column;
  gap: var(--gap-md);
}
.modal-dialog h2 { font-size: 16px; font-weight: 600; }
.user-select, .user-new-input { width: 100%; }
.user-error { color: var(--error); font-size: 13px; }
.modal-actions { display: flex; justify-content: flex-end; gap: var(--gap-sm); }
button.primary { background: var(--accent); border-color: var(--accent); color: #fff; }
button.primary:hover { filter: brightness(1.1); }
```

- [ ] **Step 3: Commit** (behavior covered by Task 10 Playwright smoke)

```bash
git add packages/HimalayaUI/frontend/src/components/UserModal.tsx \
        packages/HimalayaUI/frontend/src/styles.css
git commit -m "feat(frontend): UserModal component"
```

---

### Task 9: App.tsx composition root

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/App.tsx`

Wires the data flow: load experiment + samples on mount, pass to components. Opens the user modal on first visit (no persisted username).

- [ ] **Step 1: Replace `src/App.tsx`**

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
  const username         = useAppState((s) => s.username);
  const activeSampleId   = useAppState((s) => s.activeSampleId);
  const setUsername      = useAppState((s) => s.setUsername);
  const setActiveSample  = useAppState((s) => s.setActiveSample);

  const [experiment, setExperiment] = useState<Experiment | null>(null);
  const [samples, setSamples]       = useState<Sample[]>([]);
  const [bootError, setBootError]   = useState<string | null>(null);
  const [modalOpen, setModalOpen]   = useState<boolean>(!username);

  // Sync api.ts's X-Username state with the store.
  useEffect(() => { api.setUsername(username); }, [username]);

  // Boot data
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
        + (activeSample ? ` › ${activeSample.label ?? ""} ${activeSample.name ?? ""}`.trimEnd() : "");

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

- [ ] **Step 2: Verify types**

```bash
cd packages/HimalayaUI/frontend && npx tsc --noEmit
```

- [ ] **Step 3: Quick dev smoke (optional)** — in one terminal run the Julia server on 8080 against the Plan 2 test experiment, in another `npm run dev`. Browse to http://localhost:5173.

- [ ] **Step 4: Run unit tests and commit**

```bash
cd packages/HimalayaUI/frontend && npm test
cd ../../..
git add packages/HimalayaUI/frontend/src/App.tsx
git commit -m "feat(frontend): App composition root with data boot and state wiring"
```

---

### Task 10: Playwright E2E smoke with mocked /api

**Files:**
- Create: `packages/HimalayaUI/frontend/playwright.config.ts`
- Create: `packages/HimalayaUI/frontend/e2e/smoke.spec.ts`

Playwright runs against the Vite dev server. All `/api/*` requests are mocked via `page.route`. Keeps Plan 3 self-contained — live-backend E2E comes in a later plan once the full UI exists.

- [ ] **Step 1: Install Playwright browsers** (first run only)

```bash
cd packages/HimalayaUI/frontend && npx playwright install chromium
```

- [ ] **Step 2: Write `playwright.config.ts`**

```ts
import { defineConfig } from "@playwright/test";

export default defineConfig({
  testDir: "./e2e",
  timeout: 15_000,
  use: {
    baseURL: "http://127.0.0.1:5173",
    browserName: "chromium",
    headless: true,
  },
  webServer: {
    command: "npm run dev",
    url: "http://127.0.0.1:5173",
    reuseExistingServer: !process.env["CI"],
    stdout: "pipe",
    stderr: "pipe",
  },
});
```

- [ ] **Step 3: Write E2E tests**

Create `e2e/smoke.spec.ts`:
```ts
import { test, expect, type Page } from "@playwright/test";

const EXP = {
  id: 1, name: "TestRun", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-01-01",
};
const SAMPLES = [
  { id: 1, experiment_id: 1, label: "D1", name: "UX1", notes: null,
    tags: [{ id: 1, key: "lipid", value: "DOPC", source: "manual" }] },
  { id: 2, experiment_id: 1, label: "D2", name: "UX2", notes: null, tags: [] },
];

async function mockApi(page: Page, users: { id: number; username: string }[] = []): Promise<void> {
  await page.route("**/api/users", async (route) => {
    const req = route.request();
    if (req.method() === "GET") {
      await route.fulfill({ status: 200, contentType: "application/json",
        body: JSON.stringify(users) });
    } else if (req.method() === "POST") {
      const data = JSON.parse(req.postData() ?? "{}");
      await route.fulfill({ status: 201, contentType: "application/json",
        body: JSON.stringify({ id: users.length + 1, username: data.username }) });
    } else {
      await route.continue();
    }
  });
  await page.route("**/api/experiments/1", async (route) =>
    route.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(EXP) }));
  await page.route("**/api/experiments/1/samples", async (route) =>
    route.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(SAMPLES) }));
}

test.beforeEach(async ({ page }) => {
  await page.addInitScript(() => { localStorage.clear(); });
});

test("shell loads with navbar, layout placeholders, and sample list", async ({ page }) => {
  await mockApi(page);
  await page.goto("/");

  await expect(page.locator(".nav-logo")).toHaveText("Himalaya");
  await expect(page.locator(".placeholder")).toHaveCount(4);
  await expect(page.locator("[data-sample-id]")).toHaveCount(2);
  await expect(page.locator('[data-sample-id="1"] .sample-label')).toHaveText("D1");
});

test("first-visit prompts user modal and submits new user", async ({ page }) => {
  await mockApi(page);
  await page.goto("/");

  await expect(page.locator(".modal-dialog")).toBeVisible();
  await expect(page.locator(".modal-dialog h2")).toHaveText("Who are you?");

  await page.locator(".user-new-input").fill("alice");
  await page.locator(".user-submit").click();

  await expect(page.locator(".modal-dialog")).not.toBeVisible();
  await expect(page.locator(".nav-user")).toHaveText("alice");
});

test("clicking a sample updates the active state and breadcrumb", async ({ page }) => {
  await mockApi(page, [{ id: 1, username: "alice" }]);
  await page.addInitScript(() => {
    // Zustand `persist` stores { state, version }
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({ state: { username: "alice" }, version: 1 }),
    );
  });
  await page.goto("/");

  await expect(page.locator(".modal-dialog")).not.toBeVisible();
  await page.locator('[data-sample-id="2"]').click();
  await expect(page.locator('[data-sample-id="2"]')).toHaveClass(/active/);
  await expect(page.locator(".nav-breadcrumb")).toContainText("D2");
});

test("filter narrows the sample list", async ({ page }) => {
  await mockApi(page, [{ id: 1, username: "alice" }]);
  await page.addInitScript(() => {
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({ state: { username: "alice" }, version: 1 }),
    );
  });
  await page.goto("/");

  await page.locator(".sample-filter").fill("DOPC");
  await expect(page.locator("[data-sample-id]")).toHaveCount(1);
  await expect(page.locator("[data-sample-id]")).toHaveAttribute("data-sample-id", "1");
});
```

Note the Zustand localStorage shape: `{ state, version }`, not the bare object. That's Zustand's `persist` format — important to get right or the hydration silently fails.

- [ ] **Step 4: Run E2E**

```bash
cd packages/HimalayaUI/frontend && npm run e2e
```
Expected: all 4 tests pass. If any fail, inspect `playwright-report/`.

- [ ] **Step 5: Final checks**

```bash
cd packages/HimalayaUI/frontend
npx tsc --noEmit
npm test
npm run build
```
Expected: clean types, unit tests pass, `dist/` emitted.

- [ ] **Step 6: Commit**

```bash
cd ../../..
git add packages/HimalayaUI/frontend/playwright.config.ts \
        packages/HimalayaUI/frontend/e2e/smoke.spec.ts
git commit -m "test(frontend): Playwright E2E smoke with mocked /api"
```

---

## Self-Review

**Spec coverage (§5 Frontend):**

| Requirement | Covered by |
|---|---|
| React 18 + TypeScript + Vite | Task 1 |
| Zustand with persist | Task 4 |
| Dark Claude Code aesthetic via CSS tokens | Task 2 |
| localStorage state persistence | Task 4 (via Zustand persist) |
| Single-page analysis layout | Task 6 |
| Navbar (logo, page links, breadcrumb, user) | Task 5 |
| Sample list w/ filter and status dot | Task 7 |
| User identification modal | Task 8 |
| Served from `frontend/dist/` by Oxygen | Plan 2 Task 10 static mount; `npm run build` emits here |
| Vitest + RTL units + Playwright E2E | Tasks 2–4, 7, 10 |

**Intentionally deferred to Plans 4–6 (placeholder panes exist):**
- Trace viewer (Observable Plot log-log, peak markers, click-to-add, hover preview) — Plan 4
- Miller-index plot + phase panel with +/- buttons + hover preview on alternatives — Plan 5
- Properties panel (Exposures / Peaks / Tags / Notes tabs) — Plan 6

**Status-dot caveat:** Plan 3 renders all samples grey. Once Plans 4–5 wire in indices/groups, status upgrades to yellow/green via the same `data-status` attribute — no API change needed.

**Naming consistency:**
- Components: `Navbar`, `Layout`, `SampleList`, `UserModal` — PascalCase `.tsx`.
- Store hook: `useAppState` — selector-style access (`useAppState((s) => s.username)`) to minimize re-renders.
- API functions: `listUsers / createUser / getExperiment / listSamples / updateSample / addSampleTag / removeSampleTag`.
- Plan 4 will add `listExposures / listPeaks / addPeak / removePeak`, etc.

**Zustand persist gotcha:** the localStorage shape is `{ state: {...}, version: N }`, not a bare object. The Playwright tests and Task 4 unit test are explicit about this.

**No placeholders in step bodies** — every step has the complete code or exact command.

**Ready for execution.**
