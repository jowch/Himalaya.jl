# HimalayaUI Frontend Shell Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up the TypeScript + Vite frontend shell with the navigation bar, three-column layout skeleton, working sample list, and user identification modal. Center (trace viewer, properties panel) and right (Miller plot, phase panel) columns are placeholder stubs filled in by Plans 4–6.

**Architecture:** Vanilla TypeScript (no framework). Vite for dev server + production build. Observable Plot deferred to Plans 4–5. Single page for v1 — the per-sample analysis view. State lives in module-level closures with a thin `state.ts` that persists app state (active sample id, active exposure id) to `localStorage`. The API client (`api.ts`) is a thin typed layer over `fetch`. In dev, Vite proxies `/api/*` to the Julia server on :8080. In production, Oxygen.jl serves the compiled `dist/` at `/` and handles `/api/*` itself (single-origin).

**Tech Stack:** TypeScript 5, Vite 5, Vitest (unit tests), Playwright (E2E), no UI framework.

**This is Plan 3 of 6.** Prior plans (complete): Plan 1 (SQLite pipeline + CLI), Plan 2 (Oxygen REST API + `himalaya serve`). Remaining: Plan 4 (trace viewer + peak editing), Plan 5 (phase panel + Miller plot + hover preview), Plan 6 (properties panel tabs).

---

## File Map

| File | Responsibility |
|---|---|
| `packages/HimalayaUI/frontend/package.json` | npm deps + scripts |
| `packages/HimalayaUI/frontend/tsconfig.json` | TS strict config |
| `packages/HimalayaUI/frontend/vite.config.ts` | dev proxy, build → `dist/` |
| `packages/HimalayaUI/frontend/playwright.config.ts` | E2E config |
| `packages/HimalayaUI/frontend/.gitignore` | ignore `node_modules` + `dist` |
| `packages/HimalayaUI/frontend/index.html` | entry HTML |
| `packages/HimalayaUI/frontend/src/main.ts` | composition root |
| `packages/HimalayaUI/frontend/src/dom.ts` | `h()` element builder |
| `packages/HimalayaUI/frontend/src/api.ts` | typed fetch wrappers |
| `packages/HimalayaUI/frontend/src/state.ts` | app state + localStorage |
| `packages/HimalayaUI/frontend/src/styles.css` | design tokens + reset + base |
| `packages/HimalayaUI/frontend/src/components/navbar.ts` | top nav with breadcrumb + user |
| `packages/HimalayaUI/frontend/src/components/layout.ts` | three-column shell |
| `packages/HimalayaUI/frontend/src/components/sample_list.ts` | left column sample list + filter |
| `packages/HimalayaUI/frontend/src/components/user_modal.ts` | first-visit + switch-user modal |
| `packages/HimalayaUI/frontend/test/api.test.ts` | Vitest |
| `packages/HimalayaUI/frontend/test/state.test.ts` | Vitest |
| `packages/HimalayaUI/frontend/e2e/smoke.spec.ts` | Playwright |

---

## Dev workflow summary

```bash
cd packages/HimalayaUI/frontend
npm install               # first time
npm run dev               # Vite on :5173, proxies /api → :8080
npm test                  # Vitest unit tests
npm run build             # emit packages/HimalayaUI/frontend/dist/
npx playwright install    # first-time E2E browser install
npm run e2e               # Playwright smoke
```

`dist/` is gitignored. Oxygen's `serve(db)` will pick it up on disk if it exists and serve it at `/`.

---

## Design tokens (referenced throughout)

Dark, minimal, Claude Code desktop feel. Semantic CSS variables declared once in `styles.css`:

```
--bg: #1a1a1a              --fg: #e8e6e3
--bg-elevated: #232323     --fg-muted: #9a9894
--bg-hover: #2a2a2a        --fg-dim: #666666
--border: #333333          --accent: #d97706
--success: #22c55e         --warning: #eab308        --error: #ef4444
--radius: 6px              --gap-sm: 8px              --gap-md: 16px    --gap-lg: 24px
--font-sans: -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif
--font-mono: ui-monospace, SFMono-Regular, Menlo, monospace
```

---

### Task 1: Scaffold Vite + TypeScript

**Files:**
- Create: `packages/HimalayaUI/frontend/package.json`
- Create: `packages/HimalayaUI/frontend/tsconfig.json`
- Create: `packages/HimalayaUI/frontend/vite.config.ts`
- Create: `packages/HimalayaUI/frontend/.gitignore`
- Create: `packages/HimalayaUI/frontend/index.html`
- Delete: `packages/HimalayaUI/frontend/dist/.gitkeep`
- Create: `packages/HimalayaUI/frontend/src/main.ts` (minimal — just inserts "Hello" to verify build)

- [ ] **Step 1: Initialize package.json**

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
  "devDependencies": {
    "typescript": "^5.6.0",
    "vite": "^5.4.0",
    "vitest": "^2.1.0",
    "jsdom": "^25.0.0",
    "@playwright/test": "^1.48.0"
  }
}
```

- [ ] **Step 2: tsconfig.json**

```json
{
  "compilerOptions": {
    "target": "ES2022",
    "module": "ESNext",
    "moduleResolution": "bundler",
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
    "resolveJsonModule": true,
    "types": ["vitest/globals"]
  },
  "include": ["src/**/*.ts", "test/**/*.ts", "e2e/**/*.ts", "*.ts"]
}
```

- [ ] **Step 3: vite.config.ts**

```ts
import { defineConfig } from "vite";

export default defineConfig({
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
  test: {
    environment: "jsdom",
    globals: true,
    include: ["test/**/*.test.ts"],
  },
});
```

- [ ] **Step 4: .gitignore**

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

- [ ] **Step 5: index.html**

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
    <script type="module" src="/src/main.ts"></script>
  </body>
</html>
```

- [ ] **Step 6: minimal src/main.ts**

```ts
const app = document.getElementById("app");
if (app) app.textContent = "Himalaya UI — loading…";
```

- [ ] **Step 7: Install and build**

```bash
cd packages/HimalayaUI/frontend
npm install
npm run build
```

Expected: `dist/index.html` and `dist/assets/*.js` produced. No TS errors.

- [ ] **Step 8: Commit**

```bash
cd ../../..
git add packages/HimalayaUI/frontend/package.json \
        packages/HimalayaUI/frontend/package-lock.json \
        packages/HimalayaUI/frontend/tsconfig.json \
        packages/HimalayaUI/frontend/vite.config.ts \
        packages/HimalayaUI/frontend/.gitignore \
        packages/HimalayaUI/frontend/index.html \
        packages/HimalayaUI/frontend/src/main.ts
git rm -f packages/HimalayaUI/frontend/dist/.gitkeep 2>/dev/null || true
git commit -m "feat(frontend): scaffold Vite + TypeScript frontend"
```

---

### Task 2: Design tokens CSS + DOM helper + Vitest bootstrap

**Files:**
- Create: `packages/HimalayaUI/frontend/src/styles.css`
- Create: `packages/HimalayaUI/frontend/src/dom.ts`
- Create: `packages/HimalayaUI/frontend/test/dom.test.ts`
- Modify: `packages/HimalayaUI/frontend/src/main.ts` (import styles)

`h()` is a terse `createElement` wrapper — the only DOM abstraction we need for vanilla TS.

- [ ] **Step 1: Write failing test**

Create `test/dom.test.ts`:
```ts
import { describe, it, expect } from "vitest";
import { h } from "../src/dom";

describe("h()", () => {
  it("creates an element with a tag", () => {
    const el = h("div");
    expect(el.tagName).toBe("DIV");
  });

  it("applies props (className, textContent, dataset, on*)", () => {
    const clicks: string[] = [];
    const el = h("button", {
      className: "primary",
      textContent: "click",
      dataset: { id: "42" },
      onclick: () => clicks.push("x"),
    });
    expect(el.className).toBe("primary");
    expect(el.textContent).toBe("click");
    expect(el.dataset["id"]).toBe("42");
    el.click();
    expect(clicks).toEqual(["x"]);
  });

  it("appends children", () => {
    const el = h("ul", {}, [h("li", { textContent: "a" }), h("li", { textContent: "b" })]);
    expect(el.children.length).toBe(2);
    expect(el.children[0]!.textContent).toBe("a");
  });

  it("accepts string children", () => {
    const el = h("p", {}, ["hello ", h("b", { textContent: "world" })]);
    expect(el.textContent).toBe("hello world");
  });
});
```

- [ ] **Step 2: Verify failure**

```bash
cd packages/HimalayaUI/frontend && npm test
```
Expected: cannot resolve `../src/dom`.

- [ ] **Step 3: Implement `src/dom.ts`**

```ts
type Props<T extends HTMLElement> = Partial<
  Omit<T, "style" | "dataset" | "children">
> & {
  style?: Partial<CSSStyleDeclaration>;
  dataset?: Record<string, string>;
};

type Child = Node | string;

export function h<K extends keyof HTMLElementTagNameMap>(
  tag: K,
  props: Props<HTMLElementTagNameMap[K]> = {},
  children: Child[] = [],
): HTMLElementTagNameMap[K] {
  const el = document.createElement(tag);
  for (const [key, value] of Object.entries(props)) {
    if (value === undefined) continue;
    if (key === "style" && value && typeof value === "object") {
      Object.assign(el.style, value);
    } else if (key === "dataset" && value && typeof value === "object") {
      for (const [dk, dv] of Object.entries(value)) {
        el.dataset[dk] = String(dv);
      }
    } else {
      // Event handlers (onclick etc.) and plain DOM props
      (el as unknown as Record<string, unknown>)[key] = value;
    }
  }
  for (const c of children) {
    el.append(typeof c === "string" ? document.createTextNode(c) : c);
  }
  return el;
}

/** Replace all children of `parent` with `nodes`. */
export function mount(parent: HTMLElement, ...nodes: Child[]): void {
  parent.replaceChildren(
    ...nodes.map((n) => (typeof n === "string" ? document.createTextNode(n) : n)),
  );
}
```

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

- [ ] **Step 5: Update main.ts to import styles**

```ts
import "./styles.css";

const app = document.getElementById("app");
if (app) app.textContent = "Himalaya UI — loading…";
```

- [ ] **Step 6: Run tests**

```bash
cd packages/HimalayaUI/frontend && npm test
```
Expected: dom.test passes.

- [ ] **Step 7: Commit**

```bash
cd ../../..
git add packages/HimalayaUI/frontend/src/styles.css \
        packages/HimalayaUI/frontend/src/dom.ts \
        packages/HimalayaUI/frontend/src/main.ts \
        packages/HimalayaUI/frontend/test/dom.test.ts
git commit -m "feat(frontend): design tokens, h() DOM helper, Vitest bootstrap"
```

---

### Task 3: Typed API client

**Files:**
- Create: `packages/HimalayaUI/frontend/src/api.ts`
- Create: `packages/HimalayaUI/frontend/test/api.test.ts`

Plan 3 only uses: users (list/create), experiments (get), samples (list/patch/tag). Stubs for peaks/exposures/indices/groups/export are added in Plans 4–6.

- [ ] **Step 1: Write failing test**

Create `test/api.test.ts`:
```ts
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import * as api from "../src/api";

function mockFetch(response: { status?: number; body: unknown }): void {
  const status = response.status ?? 200;
  global.fetch = vi.fn(async (_input: RequestInfo | URL, init?: RequestInit) => {
    return new Response(JSON.stringify(response.body), {
      status,
      headers: { "Content-Type": "application/json" },
    }) as unknown as Response;
  }) as unknown as typeof fetch;
  return;
}

describe("api client", () => {
  beforeEach(() => {
    api.setUsername(undefined);
  });
  afterEach(() => {
    vi.restoreAllMocks();
  });

  it("listUsers returns [] for empty response", async () => {
    mockFetch({ body: [] });
    expect(await api.listUsers()).toEqual([]);
  });

  it("listUsers parses user rows", async () => {
    mockFetch({ body: [{ id: 1, username: "alice" }] });
    const users = await api.listUsers();
    expect(users.length).toBe(1);
    expect(users[0]!.username).toBe("alice");
  });

  it("createUser POSTs JSON with Content-Type", async () => {
    const spy = vi.fn(async () => new Response(JSON.stringify({ id: 1, username: "bob" }), {
      status: 201, headers: { "Content-Type": "application/json" },
    }));
    global.fetch = spy as unknown as typeof fetch;
    const user = await api.createUser("bob");
    expect(user.username).toBe("bob");
    const [, init] = spy.mock.calls[0]!;
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
    const [, init] = spy.mock.calls[0]!;
    expect((init?.headers as Record<string, string>)["X-Username"]).toBe("alice");
  });

  it("non-2xx throws ApiError with status and body", async () => {
    mockFetch({ status: 404, body: { error: "not found" } });
    await expect(api.getExperiment(99)).rejects.toMatchObject({
      status: 404,
      message: expect.stringContaining("not found"),
    });
  });

  it("getExperiment returns typed row", async () => {
    mockFetch({
      body: { id: 1, name: "E1", path: "/p", data_dir: "/p/data",
              analysis_dir: "/p/analysis", manifest_path: null, created_at: "..." },
    });
    const exp = await api.getExperiment(1);
    expect(exp.name).toBe("E1");
  });
});
```

- [ ] **Step 2: Verify fail**

```bash
cd packages/HimalayaUI/frontend && npm test
```

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

// ——— Users ———
export const listUsers   = () => request<User[]>("GET", "/api/users");
export const createUser  = (username: string) =>
  request<User>("POST", "/api/users", { username });

// ——— Experiments ———
export const getExperiment  = (id: number) =>
  request<Experiment>("GET", `/api/experiments/${id}`);
export const updateExperiment = (id: number, patch: Partial<Pick<Experiment,
  "name" | "data_dir" | "analysis_dir" | "manifest_path">>) =>
  request<Experiment>("PATCH", `/api/experiments/${id}`, patch);

// ——— Samples ———
export const listSamples    = (experiment_id: number) =>
  request<Sample[]>("GET", `/api/experiments/${experiment_id}/samples`);
export const updateSample   = (id: number, patch: { name?: string; notes?: string }) =>
  request<Sample>("PATCH", `/api/samples/${id}`, patch);
export const addSampleTag   = (id: number, key: string, value: string) =>
  request<SampleTag>("POST", `/api/samples/${id}/tags`, { key, value });
export const removeSampleTag = (id: number, tag_id: number) =>
  request<void>("DELETE", `/api/samples/${id}/tags/${tag_id}`);
```

- [ ] **Step 4: Run tests**

```bash
cd packages/HimalayaUI/frontend && npm test
```
Expected: all api tests pass.

- [ ] **Step 5: Commit**

```bash
cd ../../..
git add packages/HimalayaUI/frontend/src/api.ts \
        packages/HimalayaUI/frontend/test/api.test.ts
git commit -m "feat(frontend): typed API client with X-Username handling"
```

---

### Task 4: State module (localStorage)

**Files:**
- Create: `packages/HimalayaUI/frontend/src/state.ts`
- Create: `packages/HimalayaUI/frontend/test/state.test.ts`

`state.ts` is a tiny pub/sub store with localStorage persistence. Listeners are invoked on any state change.

- [ ] **Step 1: Write failing test**

Create `test/state.test.ts`:
```ts
import { describe, it, expect, beforeEach, vi } from "vitest";
import { createState, LS_KEY } from "../src/state";

describe("state", () => {
  beforeEach(() => { localStorage.clear(); });

  it("has sensible defaults", () => {
    const s = createState();
    expect(s.get().username).toBeUndefined();
    expect(s.get().activeSampleId).toBeUndefined();
    expect(s.get().activeExposureId).toBeUndefined();
  });

  it("hydrates from localStorage", () => {
    localStorage.setItem(LS_KEY, JSON.stringify({
      username: "alice", activeSampleId: 3, activeExposureId: 7,
    }));
    const s = createState();
    expect(s.get().username).toBe("alice");
    expect(s.get().activeSampleId).toBe(3);
    expect(s.get().activeExposureId).toBe(7);
  });

  it("persists on set", () => {
    const s = createState();
    s.set({ username: "bob" });
    const raw = JSON.parse(localStorage.getItem(LS_KEY)!);
    expect(raw.username).toBe("bob");
  });

  it("merges partial updates", () => {
    const s = createState();
    s.set({ username: "a" });
    s.set({ activeSampleId: 2 });
    expect(s.get().username).toBe("a");
    expect(s.get().activeSampleId).toBe(2);
  });

  it("notifies subscribers", () => {
    const s = createState();
    const spy = vi.fn();
    const unsub = s.subscribe(spy);
    s.set({ username: "c" });
    expect(spy).toHaveBeenCalledTimes(1);
    expect(spy.mock.calls[0]![0].username).toBe("c");
    unsub();
    s.set({ username: "d" });
    expect(spy).toHaveBeenCalledTimes(1);
  });

  it("ignores corrupt localStorage", () => {
    localStorage.setItem(LS_KEY, "{not json");
    const s = createState();
    expect(s.get().username).toBeUndefined();
  });
});
```

- [ ] **Step 2: Verify fail**

- [ ] **Step 3: Implement `src/state.ts`**

```ts
export interface AppState {
  username?: string;
  activeSampleId?: number;
  activeExposureId?: number;
}

export const LS_KEY = "himalaya-ui:state";

export interface Store {
  get(): AppState;
  set(patch: Partial<AppState>): void;
  subscribe(fn: (s: AppState) => void): () => void;
}

export function createState(): Store {
  let state: AppState = loadFromStorage();
  const listeners = new Set<(s: AppState) => void>();

  return {
    get: () => state,
    set(patch) {
      state = { ...state, ...patch };
      persist(state);
      for (const fn of listeners) fn(state);
    },
    subscribe(fn) {
      listeners.add(fn);
      return () => { listeners.delete(fn); };
    },
  };
}

function loadFromStorage(): AppState {
  const raw = localStorage.getItem(LS_KEY);
  if (!raw) return {};
  try {
    const parsed = JSON.parse(raw);
    return (parsed && typeof parsed === "object") ? parsed as AppState : {};
  } catch {
    return {};
  }
}

function persist(s: AppState): void {
  localStorage.setItem(LS_KEY, JSON.stringify(s));
}
```

- [ ] **Step 4: Run tests**

- [ ] **Step 5: Commit**

```bash
cd ../../..
git add packages/HimalayaUI/frontend/src/state.ts \
        packages/HimalayaUI/frontend/test/state.test.ts
git commit -m "feat(frontend): app state store with localStorage persistence"
```

---

### Task 5: Navbar component

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/navbar.ts`

Nav bar layout (left → right): **logo** · **page links** (just "Analysis" active in v1) · **breadcrumb** (center) · **user** (right, clickable). Emits callbacks via options — no global coupling.

- [ ] **Step 1: Implement `src/components/navbar.ts`**

```ts
import { h, mount } from "../dom";

export interface NavbarOptions {
  onUserClick: () => void;
}

export interface NavbarApi {
  el: HTMLElement;
  setBreadcrumb(crumb: string): void;
  setUsername(name: string | undefined): void;
}

export function renderNavbar(opts: NavbarOptions): NavbarApi {
  const breadcrumbEl = h("span", { className: "nav-breadcrumb" });
  const userBtn = h("button", {
    className: "nav-user",
    textContent: "Sign in",
    onclick: opts.onUserClick,
  });

  const el = h("header", { className: "nav" }, [
    h("div", { className: "nav-left" }, [
      h("span", { className: "nav-logo", textContent: "Himalaya" }),
      h("nav", { className: "nav-links" }, [
        h("a", { className: "nav-link active", textContent: "Analysis", href: "#" }),
      ]),
    ]),
    h("div", { className: "nav-center" }, [breadcrumbEl]),
    h("div", { className: "nav-right" }, [userBtn]),
  ]);

  return {
    el,
    setBreadcrumb(crumb) { mount(breadcrumbEl, crumb); },
    setUsername(name) { userBtn.textContent = name ?? "Sign in"; },
  };
}
```

- [ ] **Step 2: Append navbar styles to `src/styles.css`**

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

- [ ] **Step 3: Commit** (no tests yet — rendered shape is exercised in Task 10's Playwright smoke)

```bash
git add packages/HimalayaUI/frontend/src/components/navbar.ts \
        packages/HimalayaUI/frontend/src/styles.css
git commit -m "feat(frontend): navbar component with breadcrumb and user button"
```

---

### Task 6: Three-column layout shell

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/layout.ts`

Skeleton with named slots that Plans 4–6 will populate. Right now, center and right columns show placeholder text.

- [ ] **Step 1: Implement `src/components/layout.ts`**

```ts
import { h } from "../dom";

export interface LayoutApi {
  el: HTMLElement;
  leftSlot: HTMLElement;
  centerTopSlot: HTMLElement;
  centerBottomSlot: HTMLElement;
  rightTopSlot: HTMLElement;
  rightBottomSlot: HTMLElement;
}

export function renderLayout(): LayoutApi {
  const leftSlot         = h("aside", { className: "col col-left" });
  const centerTopSlot    = h("section", { className: "pane pane-center-top" },
    [h("p", { className: "muted placeholder", textContent: "Trace viewer (Plan 4)" })]);
  const centerBottomSlot = h("section", { className: "pane pane-center-bottom" },
    [h("p", { className: "muted placeholder", textContent: "Properties panel (Plan 6)" })]);
  const rightTopSlot     = h("section", { className: "pane pane-right-top" },
    [h("p", { className: "muted placeholder", textContent: "Miller-index plot (Plan 5)" })]);
  const rightBottomSlot  = h("section", { className: "pane pane-right-bottom" },
    [h("p", { className: "muted placeholder", textContent: "Phase panel (Plan 5)" })]);

  const el = h("main", { className: "layout" }, [
    leftSlot,
    h("div", { className: "col col-center" }, [centerTopSlot, centerBottomSlot]),
    h("div", { className: "col col-right" }, [rightTopSlot, rightBottomSlot]),
  ]);

  return { el, leftSlot, centerTopSlot, centerBottomSlot, rightTopSlot, rightBottomSlot };
}
```

- [ ] **Step 2: Append layout styles**

Append to `src/styles.css`:
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
git add packages/HimalayaUI/frontend/src/components/layout.ts \
        packages/HimalayaUI/frontend/src/styles.css
git commit -m "feat(frontend): three-column layout shell with named slots"
```

---

### Task 7: Sample list view + filter

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/sample_list.ts`
- Create: `packages/HimalayaUI/frontend/test/sample_list.test.ts`

Renders samples in the left column. Status dot logic for v1: all samples render with `status="unanalyzed"` (grey) unless they have at least one index (follow-up work when indices are wired). Plan 3 renders grey unconditionally — the data for analyzed/confirmed status comes from Plans 4-5. A `data-status` attribute is set so Plans 4-5 can upgrade it without modifying this file's API.

Filter: single text input matches against `label`, `name`, or any tag's `key`/`value`.

- [ ] **Step 1: Write failing test**

Create `test/sample_list.test.ts`:
```ts
import { describe, it, expect, vi } from "vitest";
import { renderSampleList } from "../src/components/sample_list";
import type { Sample } from "../src/api";

const SAMPLES: Sample[] = [
  { id: 1, experiment_id: 1, label: "D1", name: "UX1", notes: null,
    tags: [{ id: 1, key: "lipid", value: "DOPC", source: "manual" }] },
  { id: 2, experiment_id: 1, label: "D2", name: "UX2", notes: null, tags: [] },
  { id: 3, experiment_id: 1, label: "D3", name: "UL1", notes: null,
    tags: [{ id: 2, key: "peptide", value: "melittin", source: "manual" }] },
];

describe("renderSampleList", () => {
  it("renders one row per sample", () => {
    const onSelect = vi.fn();
    const list = renderSampleList({ onSelect });
    list.setSamples(SAMPLES);
    const rows = list.el.querySelectorAll("[data-sample-id]");
    expect(rows.length).toBe(3);
  });

  it("marks the active sample", () => {
    const list = renderSampleList({ onSelect: () => {} });
    list.setSamples(SAMPLES);
    list.setActive(2);
    const active = list.el.querySelector(".active");
    expect(active?.getAttribute("data-sample-id")).toBe("2");
  });

  it("calls onSelect with sample id on click", () => {
    const onSelect = vi.fn();
    const list = renderSampleList({ onSelect });
    list.setSamples(SAMPLES);
    const row = list.el.querySelector<HTMLElement>('[data-sample-id="3"]')!;
    row.click();
    expect(onSelect).toHaveBeenCalledWith(3);
  });

  it("filter narrows by label, name, and tag", () => {
    const list = renderSampleList({ onSelect: () => {} });
    list.setSamples(SAMPLES);

    list.setFilter("D1");
    expect(list.el.querySelectorAll("[data-sample-id]").length).toBe(1);

    list.setFilter("UX");
    expect(list.el.querySelectorAll("[data-sample-id]").length).toBe(2);

    list.setFilter("lipid");
    expect(list.el.querySelectorAll("[data-sample-id]").length).toBe(1);

    list.setFilter("melittin");
    expect(list.el.querySelectorAll("[data-sample-id]").length).toBe(1);

    list.setFilter("");
    expect(list.el.querySelectorAll("[data-sample-id]").length).toBe(3);
  });
});
```

- [ ] **Step 2: Verify fail**

- [ ] **Step 3: Implement `src/components/sample_list.ts`**

```ts
import { h, mount } from "../dom";
import type { Sample } from "../api";

export interface SampleListOptions {
  onSelect: (sampleId: number) => void;
}

export interface SampleListApi {
  el: HTMLElement;
  setSamples(samples: Sample[]): void;
  setActive(id: number | undefined): void;
  setFilter(text: string): void;
}

type SampleStatus = "unanalyzed" | "candidates" | "confirmed";

export function renderSampleList(opts: SampleListOptions): SampleListApi {
  let samples: Sample[] = [];
  let activeId: number | undefined = undefined;
  let filter = "";

  const listEl = h("ul", { className: "sample-list" });
  const filterInput = h("input", {
    className: "sample-filter",
    placeholder: "Filter samples…",
    type: "text",
    oninput: () => { filter = filterInput.value; render(); },
  });

  const el = h("div", { className: "sample-list-wrap" }, [
    h("div", { className: "sample-list-header" }, [filterInput]),
    listEl,
  ]);

  function matches(s: Sample): boolean {
    if (!filter) return true;
    const needle = filter.toLowerCase();
    if (s.label?.toLowerCase().includes(needle)) return true;
    if (s.name?.toLowerCase().includes(needle))  return true;
    return s.tags.some(
      (t) => t.key.toLowerCase().includes(needle) || t.value.toLowerCase().includes(needle),
    );
  }

  function render(): void {
    const rows = samples.filter(matches).map((s) => {
      const status: SampleStatus = "unanalyzed";   // upgraded by later plans
      const row = h("li", {
        className: "sample-row" + (s.id === activeId ? " active" : ""),
        dataset: { sampleId: String(s.id), status },
        onclick: () => opts.onSelect(s.id),
      }, [
        h("span", { className: `status-dot status-${status}` }),
        h("span", { className: "sample-label", textContent: s.label ?? "" }),
        h("span", { className: "sample-name muted", textContent: s.name ?? "" }),
      ]);
      return row;
    });
    mount(listEl, ...rows);
  }

  return {
    el,
    setSamples(s)  { samples = s; render(); },
    setActive(id)  { activeId = id; render(); },
    setFilter(t)   { filter = t; filterInput.value = t; render(); },
  };
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

- [ ] **Step 5: Run tests**

```bash
cd packages/HimalayaUI/frontend && npm test
```

- [ ] **Step 6: Commit**

```bash
cd ../../..
git add packages/HimalayaUI/frontend/src/components/sample_list.ts \
        packages/HimalayaUI/frontend/src/styles.css \
        packages/HimalayaUI/frontend/test/sample_list.test.ts
git commit -m "feat(frontend): sample list with filter and click-to-select"
```

---

### Task 8: User identification modal

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/user_modal.ts`

Modal opens when no user is set or when the user clicks the navbar user button. Behavior:
- Dropdown of known usernames from `listUsers()`.
- An option "New user…" reveals a text input.
- Selection calls `api.createUser()` (idempotent on the server) and reports back via `onSelect`.
- Escape key or backdrop click closes (but does not pick a user).

- [ ] **Step 1: Implement `src/components/user_modal.ts`**

```ts
import { h, mount } from "../dom";
import * as api from "../api";

export interface UserModalApi {
  el: HTMLElement;
  open(): Promise<void>;
  close(): void;
}

export interface UserModalOptions {
  onSelect: (username: string) => void;
}

export function renderUserModal(opts: UserModalOptions): UserModalApi {
  const select   = h("select", { className: "user-select" });
  const newInput = h("input", {
    className: "user-new-input",
    type: "text",
    placeholder: "Enter username",
    style: { display: "none" },
  });
  const errorEl = h("p", { className: "user-error", style: { display: "none" } });
  const submitBtn = h("button", { className: "user-submit primary", textContent: "Continue" });

  const dialog = h("div", { className: "modal-dialog", role: "dialog", ariaModal: "true" } as Partial<HTMLDivElement>, [
    h("h2", { textContent: "Who are you?" }),
    h("p", { className: "muted", textContent: "Your name is stored with every change you make so others can see what you've done." }),
    select,
    newInput,
    errorEl,
    h("div", { className: "modal-actions" }, [submitBtn]),
  ]);

  const backdrop = h("div", { className: "modal-backdrop", style: { display: "none" } }, [dialog]);

  function showError(msg: string): void {
    errorEl.textContent = msg;
    errorEl.style.display = "block";
  }

  function hideError(): void { errorEl.style.display = "none"; }

  async function refreshOptions(): Promise<void> {
    const users = await api.listUsers();
    mount(select,
      ...users.map((u) => h("option", { value: u.username, textContent: u.username })),
      h("option", { value: "__new__", textContent: "+ New user…" }),
    );
    if (users.length === 0) {
      select.value = "__new__";
      select.dispatchEvent(new Event("change"));
    }
  }

  select.onchange = () => {
    if (select.value === "__new__") {
      newInput.style.display = "block";
      newInput.focus();
    } else {
      newInput.style.display = "none";
    }
    hideError();
  };

  let resolvePromise: ((v: void) => void) | undefined;

  submitBtn.onclick = async () => {
    hideError();
    let username = select.value;
    if (username === "__new__") {
      username = newInput.value.trim();
      if (!username) { showError("Username required"); return; }
    }
    try {
      await api.createUser(username);
      api.setUsername(username);
      opts.onSelect(username);
      close();
    } catch (err) {
      showError(`Failed: ${(err as Error).message}`);
    }
  };

  backdrop.onclick = (ev) => { if (ev.target === backdrop) close(); };
  window.addEventListener("keydown", (ev) => {
    if (ev.key === "Escape" && backdrop.style.display !== "none") close();
  });

  function close(): void {
    backdrop.style.display = "none";
    resolvePromise?.();
    resolvePromise = undefined;
  }

  return {
    el: backdrop,
    async open() {
      await refreshOptions();
      backdrop.style.display = "flex";
      return new Promise<void>((resolve) => { resolvePromise = resolve; });
    },
    close,
  };
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

- [ ] **Step 3: Commit** (behavior is covered by Task 10's Playwright smoke — unit testing the modal against mocked `api.*` without integration isn't meaningful)

```bash
git add packages/HimalayaUI/frontend/src/components/user_modal.ts \
        packages/HimalayaUI/frontend/src/styles.css
git commit -m "feat(frontend): user identification modal"
```

---

### Task 9: Wire everything in `main.ts`

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/main.ts`

- [ ] **Step 1: Replace `src/main.ts`**

```ts
import "./styles.css";
import { mount } from "./dom";
import * as api from "./api";
import { createState } from "./state";
import { renderNavbar } from "./components/navbar";
import { renderLayout } from "./components/layout";
import { renderSampleList } from "./components/sample_list";
import { renderUserModal } from "./components/user_modal";

const EXPERIMENT_ID = 1;

async function main(): Promise<void> {
  const appEl = document.getElementById("app");
  if (!appEl) throw new Error("#app root missing");

  const state = createState();
  if (state.get().username) api.setUsername(state.get().username);

  const userModal = renderUserModal({
    onSelect: (username) => {
      state.set({ username });
      navbar.setUsername(username);
    },
  });

  const navbar = renderNavbar({
    onUserClick: () => { void userModal.open(); },
  });
  navbar.setUsername(state.get().username);

  const layout = renderLayout();

  const sampleList = renderSampleList({
    onSelect: (sampleId) => {
      state.set({ activeSampleId: sampleId, activeExposureId: undefined });
      renderBreadcrumb();
    },
  });
  mount(layout.leftSlot, sampleList.el);

  mount(appEl, navbar.el, layout.el, userModal.el);

  // Data boot
  try {
    const exp = await api.getExperiment(EXPERIMENT_ID);
    const samples = await api.listSamples(EXPERIMENT_ID);
    sampleList.setSamples(samples);
    if (state.get().activeSampleId) sampleList.setActive(state.get().activeSampleId);

    function renderBreadcrumb(): void {
      const sid = state.get().activeSampleId;
      const s = samples.find((x) => x.id === sid);
      const sampleBit = s ? ` › ${s.label ?? ""} ${s.name ?? ""}`.trimEnd() : "";
      navbar.setBreadcrumb((exp.name ?? "experiment") + sampleBit);
    }
    renderBreadcrumb();
    // Make available for the sample-list onSelect closure
    (window as unknown as { renderBreadcrumb: () => void }).renderBreadcrumb = renderBreadcrumb;
  } catch (err) {
    navbar.setBreadcrumb(`Error: ${(err as Error).message}`);
  }

  if (!state.get().username) void userModal.open();
}

void main();
```

Note: the small `window.renderBreadcrumb` trick is a deliberate shortcut — Plan 4 will refactor to a cleaner observer once more views need the breadcrumb. Mark it with a `// TODO: Plan 4 refactor` comment if you like, but it's not required.

Actually, there's a cleaner fix in place: replace the closure hack by having the sample-list `onSelect` call `state.set({...})` and subscribing a `renderBreadcrumb` listener directly.

Revise `main.ts` to use state subscription (preferred):

```ts
import "./styles.css";
import { mount } from "./dom";
import * as api from "./api";
import { createState } from "./state";
import { renderNavbar } from "./components/navbar";
import { renderLayout } from "./components/layout";
import { renderSampleList } from "./components/sample_list";
import { renderUserModal } from "./components/user_modal";

const EXPERIMENT_ID = 1;

async function main(): Promise<void> {
  const appEl = document.getElementById("app");
  if (!appEl) throw new Error("#app root missing");

  const state = createState();
  if (state.get().username) api.setUsername(state.get().username);

  const userModal = renderUserModal({
    onSelect: (username) => state.set({ username }),
  });

  const navbar = renderNavbar({
    onUserClick: () => { void userModal.open(); },
  });

  const layout = renderLayout();

  const sampleList = renderSampleList({
    onSelect: (sampleId) => state.set({ activeSampleId: sampleId, activeExposureId: undefined }),
  });
  mount(layout.leftSlot, sampleList.el);

  mount(appEl, navbar.el, layout.el, userModal.el);

  let experimentName = "experiment";
  let samples: api.Sample[] = [];

  function render(): void {
    const s = state.get();
    navbar.setUsername(s.username);
    const active = samples.find((x) => x.id === s.activeSampleId);
    const bit = active ? ` › ${active.label ?? ""} ${active.name ?? ""}`.trimEnd() : "";
    navbar.setBreadcrumb(experimentName + bit);
    sampleList.setActive(s.activeSampleId);
  }
  state.subscribe(render);

  try {
    const exp = await api.getExperiment(EXPERIMENT_ID);
    experimentName = exp.name ?? "experiment";
    samples = await api.listSamples(EXPERIMENT_ID);
    sampleList.setSamples(samples);
    render();
  } catch (err) {
    navbar.setBreadcrumb(`Error: ${(err as Error).message}`);
  }

  if (!state.get().username) void userModal.open();
}

void main();
```

Use **the second** (revised) `main.ts`. Discard the first example.

- [ ] **Step 2: Verify types**

```bash
cd packages/HimalayaUI/frontend && npx tsc --noEmit
```
Expected: no errors.

- [ ] **Step 3: Quick manual dev-server sanity check** (optional, skip if blocked)

In one terminal:
```bash
cd packages/HimalayaUI && julia --project=. -e '
using HimalayaUI
# Reuse /tmp/himalaya_test from Plan 2 smoke, or reseed
main(["serve", "/tmp/himalaya_test", "--port", "8080"])'
```

In another:
```bash
cd packages/HimalayaUI/frontend && npm run dev
```
Navigate to http://localhost:5173. Expect: navbar, sample list populated (D1 UX1), user modal prompting for identification, layout placeholders visible.

- [ ] **Step 4: Commit**

```bash
cd ../../..
git add packages/HimalayaUI/frontend/src/main.ts
git commit -m "feat(frontend): wire navbar, layout, sample list, user modal in main.ts"
```

---

### Task 10: Playwright E2E smoke with mocked API

**Files:**
- Create: `packages/HimalayaUI/frontend/playwright.config.ts`
- Create: `packages/HimalayaUI/frontend/e2e/smoke.spec.ts`

Playwright runs against the Vite dev server. All `/api/*` requests are mocked via `page.route`. This keeps Plan 3 self-contained — live-backend E2E will be added in a later plan once the full UI exists.

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

- [ ] **Step 3: Write failing E2E test**

Create `e2e/smoke.spec.ts`:
```ts
import { test, expect } from "@playwright/test";

const EXP = {
  id: 1, name: "TestRun", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-01-01",
};
const SAMPLES = [
  { id: 1, experiment_id: 1, label: "D1", name: "UX1", notes: null,
    tags: [{ id: 1, key: "lipid", value: "DOPC", source: "manual" }] },
  { id: 2, experiment_id: 1, label: "D2", name: "UX2", notes: null, tags: [] },
];

async function mockApi(page: import("@playwright/test").Page, users: { id: number; username: string }[] = []): Promise<void> {
  await page.route("**/api/users", async (route) => {
    if (route.request().method() === "GET") {
      await route.fulfill({ status: 200, contentType: "application/json", body: JSON.stringify(users) });
    } else if (route.request().method() === "POST") {
      const data = JSON.parse(route.request().postData() ?? "{}");
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

  // Modal should be visible on first visit
  await expect(page.locator(".modal-dialog")).toBeVisible();
  await expect(page.locator(".modal-dialog h2")).toHaveText("Who are you?");

  // Select "+ New user…" (auto-selected if empty list) and type
  await page.locator(".user-new-input").fill("alice");
  await page.locator(".user-submit").click();

  // Modal hides, navbar shows username
  await expect(page.locator(".modal-dialog")).not.toBeVisible();
  await expect(page.locator(".nav-user")).toHaveText("alice");
});

test("clicking a sample updates the active state and breadcrumb", async ({ page }) => {
  await mockApi(page, [{ id: 1, username: "alice" }]);
  await page.addInitScript(() => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({ username: "alice" }));
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
    localStorage.setItem("himalaya-ui:state", JSON.stringify({ username: "alice" }));
  });
  await page.goto("/");

  await page.locator(".sample-filter").fill("DOPC");
  await expect(page.locator("[data-sample-id]")).toHaveCount(1);
  await expect(page.locator("[data-sample-id]")).toHaveAttribute("data-sample-id", "1");
});
```

- [ ] **Step 4: Verify failing, then run**

```bash
cd packages/HimalayaUI/frontend && npm run e2e
```
If Vite isn't already running, Playwright's `webServer` will start it. Expected: all 4 tests pass.

If any fail:
- Inspect `playwright-report/` (auto-generated on failure)
- Check that selectors match the actual component class names from Tasks 5–8

- [ ] **Step 5: Final type check + unit test pass**

```bash
cd packages/HimalayaUI/frontend
npx tsc --noEmit
npm test
```

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
| TypeScript + Vite | Task 1 |
| Dark Claude Code aesthetic via CSS tokens | Task 2 |
| localStorage state persistence | Task 4 |
| Single-page (analysis) layout | Task 6 |
| Navbar with logo, page links, breadcrumb, user | Task 5 |
| Sample list w/ filter and status dot | Task 7 |
| User identification modal | Task 8 |
| Served from `frontend/dist/` by Oxygen | Plan 2 Task 10 static mount; `npm run build` emits here |
| Vitest units + Playwright E2E | Tasks 2–4, 7, 10 |

**Intentionally deferred to Plans 4–6 (placeholder panes exist):**
- Trace viewer (Observable Plot log-log, peak markers, click-to-add, hover preview) — Plan 4
- Miller-index plot + phase panel with +/- buttons + hover preview on alternatives — Plan 5
- Properties panel (Exposures / Peaks / Tags / Notes tabs) — Plan 6

**Status-dot caveat:** Plan 3 renders all samples as unanalyzed (grey). Once Plans 4–5 wire in indices/groups, status will upgrade to candidates (yellow) or confirmed (green). The `data-status` attribute on sample rows is the extension point — no API change needed in `sample_list.ts`.

**Naming consistency:**
- `createState()`, `Store.get/set/subscribe` — consistent.
- `renderNavbar / renderLayout / renderSampleList / renderUserModal` — one factory per component.
- All components expose `{ el, ...api }` objects; slots are assigned via `mount(layout.leftSlot, component.el)`.
- API functions: `listUsers / createUser / getExperiment / listSamples / updateSample / addSampleTag / removeSampleTag`. Plan 4 will add `listExposures / listPeaks / addPeak / removePeak`, etc.

**No placeholders in step bodies** — every step has the complete code or exact command.

**Ready for execution.**
