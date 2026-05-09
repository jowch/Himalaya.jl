# Permalink URLs — Design Spec

**Issue:** #89 (depends on #88, which has landed at `8ac2bf6`)
**Status:** Amended 2026-05-09 to integrate with the existing `react-router-dom` setup discovered during plan review. The original design said "no router library"; in fact `BrowserRouter` is already mounted at `main.tsx` and Compare pages use `<Route>` + `useParams` already. The amendment below preserves the original intent (URL ↔ Zustand sync) while routing through react-router primitives so the two systems coexist coherently.
**Worktree:** `.claude/worktrees/permalinks` (branch `permalinks`)

## 1. Problem

HimalayaUI is single-page (`/`) with all selection state in Zustand + localStorage. There are no addressable URLs. A user can't share "look at this sample" / "look at this exposure" / "look at this comparison" by copying the address bar. Bookmarks are useless. Refreshing the tab keeps you where you were — but only because Zustand persists; nothing about the URL communicates intent.

This spec adds slug-based permalinks for every page in the app. A user copies the address bar; the recipient pastes; the recipient lands at the same screen.

#88 is a hard prerequisite: pre-#88, `samples.name` was UI-mutable and not unique within an experiment, so any URL keyed on it would break on rename and collide on duplicate names. Post-#88, `samples.name` is stable, convention-enforced (`^[A-Za-z0-9._-]+$`), and unique within experiment. Slug-safety for samples is a property of the data, not a runtime escape-hatch.

`experiments.name` and `exposures.filename` do **not** yet have the same guarantees (no UNIQUE constraints on either; no convention validation on `experiments.name`). §3.1 below addresses the resulting ambiguity with a deterministic tiebreaker plus one targeted schema tightening for exposures. Full experiment-name stabilization (a #88-style migration for experiments) is explicitly out of scope (§10).

## 2. URL grammar

```
/                                                    redirect via persisted Zustand state, fallback /index
/index                                               empty (NavModal at experiment step)
/index/<experiment>                                  experiment chosen, NavModal at sample step
/index/<experiment>/<sample>                         full Index workspace
/inspect                                             empty
/inspect/<experiment>                                experiment chosen
/inspect/<experiment>/<sample>                       full Inspect page
/inspect/<experiment>/<sample>?exposure=<filename>   thumbnail deep link
/compare                                             comparisons list
/compare/new                                         new draft
/compare/<id>                                        review mode
/compare/<id>/edit                                   edit mode
/<anything-else>                                     stale URL → 404 page (§6)
```

**Slugs.**
- `<experiment>` — `experiments.name`. Free-form today; legacy DBs may contain non-conformant names, including duplicates within the table. Conformant data emits verbatim; non-conformant percent-encodes at the boundary. Duplicate-name resolution: see §3.1.
- `<sample>` — `samples.name`, `^[A-Za-z0-9._-]+$` post-#88, unique within experiment.
- `<filename>` — `exposures.filename`. The column stores the **stem** (the `{name}` substitution), not the full on-disk filename — `resolve_files` strips the integration_pattern suffix at ingest (`config.jl:219`). So `exposures.filename = "JC001-007"`, not `"JC001-007_int.dat"`. This is what manifests speak in and what operators expect to see in URLs. This spec adds a `UNIQUE INDEX exposures_unique_filename ON exposures(sample_id, filename)` to make the procedural upsert contract (`cli.jl:191`) declarative and remove `/api/resolve` ambiguity (§3.1).
- `<id>` — `comparisons.id`. Numeric. Comparisons are addressed by id; a slug column is a future enhancement (§10).

**Encoding.** Conformant slugs (matching `[A-Za-z0-9._-]+`) emit and parse verbatim — no percent-encoding either way. The `parseLocation` parser uses `decodeURIComponent` on every captured slug to handle the legacy non-conformant case; the URL-emit path uses `encodeURIComponent` symmetrically.

## 3. Backend

### 3.1 `GET /api/resolve`

New file `packages/HimalayaUI/src/routes_resolve.jl`, registered in `register_routes!` between `register_picker_routes!()` and the SPA catch-all (`/api/*` precedence is unambiguous because the catch-all explicitly returns 404 for `api/`-prefixed paths — see §3.2).

```
GET /api/resolve?experiment=<name>[&sample=<name>][&exposure=<filename>]
GET /api/resolve?experiment_id=<n>[&sample_id=<m>][&exposure_id=<k>]
```

Mutually exclusive: a request that mixes name- and id-form params for the same entity (`experiment` + `experiment_id`) returns **400** with `{ error: "ambiguous_params" }`. Mixing across entities (`experiment` + `sample_id`) is allowed — name-form for some, id-form for others.

```
→ 200 {
    experiment_id, experiment_name,
    sample_id?,    sample_name?,
    exposure_id?,  exposure_filename?
  }
→ 400 { error: "ambiguous_params" }            // name+id for same entity
→ 404 {
    error: "not_found",
    missing: "experiment" | "sample" | "exposure",
    missing_value: "<value>",
    experiment_resolved?: { id, name },        // present when chain partially resolved
    sample_resolved?:     { id, name }         // present when missing == "exposure"
  }
```

Names are always echoed in the 200 response so the cold-mount `/` redirect (§5) can build the slug URL without a follow-up fetch.

**Resolution order.** Walk the chain left-to-right. First miss is the 404. Earlier-resolved entities ride along under `experiment_resolved` / `sample_resolved` (only the relevant prefix is included — e.g. missing `"sample"` carries `experiment_resolved` only; missing `"exposure"` carries both).

**Tiebreakers for non-unique names.** `experiments.name` has no UNIQUE constraint and no convention enforcement; legacy data may contain duplicates. The resolve query selects deterministically: `WHERE name = ? ORDER BY id ASC LIMIT 1`. This is the same tiebreaker rule samples used pre-#88. The behavior is documented and pinned in `test_route_response_shapes.jl`. A future #88-style experiment-name migration would supplant this (§10).

`exposures.filename` gets a new `UNIQUE INDEX exposures_unique_filename ON exposures(sample_id, filename)` to make the procedural upsert contract declarative and remove resolve ambiguity at the source. The migration follows the **dedupe-then-enforce** pattern established by `migrate_samples_naming!` (`db.jl:910–977`) — duplicates (none expected; the upsert contract has been in place since v1, but legacy data is checked anyway) are renamed deterministically before the index is created, so the `CREATE UNIQUE INDEX` always succeeds against clean data. There is no permanent safety-net fallback in the resolve query.

Implemented as a named helper `migrate_exposures_unique_filename!(db)` in `db.jl`, called from `migrate_schema!` **after `_fix_fk_references_after_autoincrement_migration!`** (`db.jl:317`). This ordering is load-bearing: `migrate_pk_to_autoincrement!` (`db.jl:313`) renames `exposures` → `_migrate_old_exposures` and rebuilds via `create_schema!` on legacy DBs; placing the index-creation any earlier would attach the index to the soon-to-be-dropped table. Running after the FK-heal step ensures the AUTOINCREMENT rebuild has settled and the index attaches to the rebuilt `exposures` table that survives.

Direct invocation from tests follows the FK-heal regression-test pattern in CLAUDE.md (helper called directly, bypassing `open_db`'s full migration chain).

```sql
-- Step 1: Find duplicates via GROUP BY HAVING COUNT(*) > 1
-- Step 2: For each (sample_id, filename) duplicate group, ordered by id ASC,
--         the first row keeps the bare filename; second-and-later get
--         <filename>-2, <filename>-3, etc., with collision avoidance against
--         user-supplied filenames already shaped that way. Emit @warn per rename.
-- Step 3: CREATE UNIQUE INDEX IF NOT EXISTS exposures_unique_filename
--           ON exposures(sample_id, filename);
```

All three steps run inside `SQLite.transaction(db) do … end` so a Ctrl-C between steps 2 and 3 cannot leave a half-renamed DB without uniqueness enforcement. Idempotent on re-run: the duplicate-suffix pass becomes a no-op once names are unique, and `CREATE UNIQUE INDEX IF NOT EXISTS` is idempotent.

**Upsert invariant.** The new index is now an additional invariant the upsert in `cli.jl:191` relies on. The existing SELECT-then-INSERT pattern remains correct (the SELECT runs first). Concurrency surface: Oxygen.jl handlers run concurrently (`server.jl:140` `async = true`), and `routes_experiments.jl:87` does NOT wrap `reingest!` in `with_idempotency` — so two simultaneous `POST /api/experiments/:id/reingest` calls for the same experiment can in principle race on the SELECT-then-INSERT and trigger a UNIQUE-constraint failure. The race is rare (operators typically do not double-fire reingest) and pre-exists this spec — `samples_unique_name` from #88 has the same shape. Wrapping reingest in `with_idempotency` is out of scope for #89, but is the right fix and should be filed as a follow-up. If reingest is ever moved to a worker thread or made parallel, the upsert must be wrapped in `INSERT OR IGNORE` + post-check or moved into a transaction with retry-on-constraint-violation.

**Read-only.** No writes; no `with_idempotency` wrapping; no `apply_event!`; no SSE emission; no `client_op_id`; no `pendingDeferreds` participation. Three SELECTs and a response. The route is queue-orthogonal by design — confirmed against `docs/mutation-queue.md`. Future maintainers tempted to extend it with a write path: please don't; add a sibling endpoint instead.

**Why a single endpoint instead of three.** A single round trip for the common `/index/<exp>/<sample>` case. Three calls would race or require client-side chaining. The 404 body shape gives the frontend everything it needs to position the recovery UX in one shot.

### 3.2 SPA catch-all

The HTTP.jl router that Oxygen.jl wraps (`HTTP/src/Handlers.jl`) splits the request path on `/` and matches each segment independently. Conditional regex captures like `{rest:.*}` apply per-segment, so `.*` cannot span slashes — `/inspect/exp/sample` would never match `/{rest:.*}`. The multi-segment wildcard is `/**` (`doublestar`); it matches any number of trailing segments and must be the final segment in the route.

`Oxygen.dynamicfiles(dist_dir, "/")` registers each on-disk file under `dist_dir` as its own exact-match route. There is no implicit fallthrough to `index.html` for unknown paths. So we need an explicit catch-all that always returns the SPA shell for non-asset paths.

Add to `server.jl::register_routes!` immediately after the existing `Oxygen.dynamicfiles(dist_dir, "/")` block (only when `isdir(dist_dir)`):

```julia
@get "/**" function(req::HTTP.Request)
    path = HTTP.URI(req.target).path        # e.g. "/inspect/exp1/sampleA"
    rest = lstrip(path, '/')
    startswith(rest, "api/") && return HTTP.Response(404, "Not found")
    return HTTP.Response(200,
        ["Content-Type" => "text/html", "Cache-Control" => "no-store"],
        read(joinpath(dist_dir, "index.html")))
end
```

- Static assets (`/foo.png`, `/assets/main.abc123.js`) are served by the exact-match routes registered by `dynamicfiles`, which precede the catch-all. The catch-all never sees them. There is no `isfile` branch in the catch-all because it would be dead code.
- The `api/` guard returns 404 for any unregistered API path. Without this guard, `/api/typo` would fall through to `index.html`, masking server bugs as 200s and breaking client error handling.
- `Cache-Control: no-store` on the HTML shell ensures users get the latest after deploy. Vite-bundled JS/CSS keep their content-hashed long-cache headers via `dynamicfiles`.

**Test-environment guard.** The catch-all only mounts when `dist_dir` exists. The Julia test harness has no frontend dist, so `routes_resolve` is exercised directly via HTTP and the SPA fallback never runs in those tests. A separate test file in §8.4 covers the fallback paths by setting `ENV["HIMALAYA_FRONTEND_DIST"]` to a synthetic dir (containing only `index.html` and one asset) before `start_test_server!`, then deleting it on teardown — `server.jl:31–35` reads this env var inside `register_routes!` so the synthetic dir is picked up cleanly.

## 4. Frontend URL-sync layer

**Router integration.** `react-router-dom` v6 is already in use — `BrowserRouter` mounts at `main.tsx`; `<Routes>` lives in `AppShell.tsx` for Compare URLs (`/experiments/:eid/compare/...` and `/compare/all/...`); `TabRocker` uses `useNavigate`. The new permalink hooks integrate with this rather than fighting it:

- `useStateFromUrl` reads the URL via `useLocation()` (subscribes to react-router's location, so `popstate` and `useNavigate` both flow through) and dispatches Zustand setters.
- `useUrlFromState` writes the URL via `useNavigate()` with `{ replace: true | false }`. No raw `history.pushState`.
- New `<Route>` declarations cover the index/inspect URL shapes; the existing `<Route path="*" element={<ZustandShellPage />} />` fallback becomes the explicit `kind: "stale"` handler.

Total surface: ~150 lines across the parser + the two hooks + the page component + ~10 new `<Route>` lines in AppShell.

### 4.1 `parseLocation(pathname, search) → ParsedUrl`

Pure function in `packages/HimalayaUI/frontend/src/lib/url/parseLocation.ts`. No side effects.

The repo's TypeScript config has `exactOptionalPropertyTypes: true`, which makes `field?: T` and `field: T | undefined` non-equivalent. The codebase consistently uses the latter (e.g. `state.ts:32` `activeExperimentId: number | undefined`). The discriminated union below follows that convention so values produced from `string | undefined` upstream values can populate the parser output without compilation friction.

```ts
export type ParsedUrl =
  | { kind: "root" }
  | { kind: "index";   experiment: string | undefined; sample: string | undefined }
  | { kind: "inspect"; experiment: string | undefined; sample: string | undefined; exposure: string | undefined }
  | { kind: "compare"; view: "list" }
  | { kind: "compare"; view: "new" }
  | { kind: "compare"; view: "review";  id: number }
  | { kind: "compare"; view: "edit";    id: number }
  | { kind: "stale";   raw: string };          // unrecognized path

export function parseLocation(pathname: string, search: string): ParsedUrl;
```

Rules:
- Empty path or `/` → `{ kind: "root" }`.
- `/<page>` where page ∈ {index, inspect, compare}: split on `/`, `decodeURIComponent` each segment, populate fields per the grammar. Missing slug positions populate as `undefined`.
- `/compare/<n>` requires `n` to be a positive integer; otherwise `{ kind: "stale", raw }`.
- Anything else → `{ kind: "stale", raw: pathname + search }`.
- For `inspect`, the `?exposure=` query param is read only when both experiment and sample are present.

The 200 / 404 response shapes in §3.1 use `field?: T` syntax in the on-the-wire JSON for compactness; the TypeScript interfaces in `api.ts` declare every optional field as `T | undefined` to keep the strict-mode contract consistent.

### 4.2 `useStateFromUrl()`

Mounted in `App.tsx` (or AppShell — anywhere under `BrowserRouter`). Reads `useLocation()` from react-router so it reacts to both `popstate` and `useNavigate`:

1. `parseLocation(location.pathname, location.search)`.
2. If `kind: "root"` → run the §5 redirect.
3. If `kind: "stale"` → set Zustand `staleUrlContext = { kind: "unknown_path", raw }`. Active-* slots (`activeExperimentId`, `activeSampleId`, `activeExposureId`) are **not** cleared — the user just pasted an unparseable URL; preserving prior context lets the recovery CTA route them back to where they were if they want.
4. Otherwise → fetch `/api/resolve?…` with whichever slugs are present.
   - 200 → dispatch Zustand setters: `setActivePage`, `setActiveExperiment`, `setActiveSample` (or `undefined`, depending on what resolved), `setActiveExposure`. Setters clear `staleUrlContext` (see §6).
   - 404 → `setStaleUrlContext(body)`. Active-* slots are **not** cleared on 404 either — same rationale. The recovery action (`recoverFromStaleUrl` in §6) decides which slots to overwrite based on what the 404 payload knows resolved, so blanket clearing here would lose useful context.

**Origin-tagged fetches.** Each fetch is tagged with the `pathname + search` it was launched against. Before applying the response (200 or 404), the hook compares the tag to the current `location.pathname + search`; if it changed (the user navigated mid-flight via TabRocker, NavModal, or another `popstate`), the response is discarded. This closes the popstate-vs-click race that a bare `AbortController` doesn't catch — Zustand-driven mutations don't trigger `popstate` so the abort signal alone is insufficient.

```ts
const origin = location.pathname + location.search;
const ctl = new AbortController();
const data = await fetch(url, { signal: ctl.signal }).then(r => r.json());
if (location.pathname + location.search !== origin) return;   // stale
applyResolveResult(data);
```

**No automatic redirect on partial 404.** The earlier "redirect to fallback" model (per the original issue text) is replaced by the 404 page; the URL stays at the broken value until the user picks something.

**Render-state during fetch.** While a resolve is in flight, AppShell renders a `<ResolvingFallback />` placeholder instead of the previous page contents. This avoids the "flash of wrong content" that would otherwise occur when a deep-link cold mount paints the previously-active page (from persisted Zustand) for ~30–80ms before the resolve completes and swaps in the right one. The placeholder is intentionally minimal: an empty container matching the existing layout chrome's flex bones (so the AppHeader and TabRocker stay rendered above it without a layout shift) marked with `data-testid="resolving"` for the §8.2 test. No `<Skeleton>` from boneyard — this is a sub-paint flicker mask, not a load-state animation; a skeleton would itself flash on every navigation.

`resolving` lives in Zustand alongside `staleUrlContext` (one slot, one source of truth, AppShell reads via a single selector). Set to `true` at the start of `useStateFromUrl`'s recognized-kind fetch; cleared on success/failure/discard. The named-action discipline of `state.ts` covers it: `setResolving(value: boolean)` is the only mutator, and the various `setActive*` setters do not touch it.

**Pre-fetch `staleUrlContext` clear.** The first thing `useStateFromUrl` does on a recognized-kind URL (anything other than `kind: "stale"`) is `setStaleUrlContext(null)`. This prevents a stale 404 page from leaking into the next render cycle when the user navigates away from a stale URL — without this, the AppShell ladder below would briefly show the OLD `<StaleUrlPage>` over the new `<ResolvingFallback />` between popstate-fire and the named-action clear that lands when the resolve 200s.

### 4.3 `useUrlFromState()`

Mounted in `App.tsx`. Subscribes via TanStack `useQuery`s to `experiments` and to `samples` for the active experiment, so name-by-id lookups happen against live cache state — and so SSE-driven cache rewrites trigger a re-render of this hook (§7 invalidation). Also subscribes to the relevant Zustand selectors via `useAppState`.

On any change, computes the target URL and emits `navigate(target, { replace })` per the policy:

| Trigger | Push or replace |
|---|---|
| TabRocker click → `activePage` change | push, target `/<page>/<lastExp>/<lastSample>` (continuity) |
| NavModal commit (closes with full selection) | push |
| Intermediate NavModal state | no URL change (modal isn't committed) |
| Inspect thumbnail click → `activeExposureId` change | replace |
| App-init URL sync (initial render after `useStateFromUrl`) | replace |
| SSE-driven cache update that invalidates the current URL | replace, then re-resolve |
| Stale-URL recovery (`recoverFromStaleUrl` or `setActivePage` from StaleUrlPage CTA) | replace (drops the broken URL from history) |

The hook **does not emit** when the resulting URL equals the current one (string compare on `pathname + search`). Prevents spurious history entries when an SSE refetch hydrates names without changing them.

**Replay-as-rerun.** During `replayCoordinator.ts`'s rollback → applyRemoteToCache → re-run-onMutate cycle, two cases are possible for the URL:
1. **Foreign event does not affect the active entity** (the common case — e.g. another user edits a different sample). Optimistic ids re-create the same name mapping; the URL recompute sees the same slug after every step; the equality guard suppresses the trivial echo. The §8.1 replay-as-rerun no-spurious-emit test asserts exactly this case.
2. **Foreign event removes the active entity** (e.g. another user deletes the sample you're viewing). The URL legitimately changes — the active sample's slug becomes `undefined`. The equality check does not suppress this, nor should it: this is the §7 SSE-driven URL invalidation flow, which is the correct UX (replace to a still-valid URL, set `staleUrlContext`). The replay test does not cover this case; the §7 live-integration test does.

**Continuity rule.** When the user TabRockers from `/compare/123` to Index, the hook reads `activeExperimentId` / `activeSampleId` from Zustand and emits `/index/<exp>/<sample>` if both are present. If either is missing, emits the bare page (`/index`). Names come from the TanStack cache via `useQuery` subscriptions, so a name change (cache hydration or refetch) triggers a re-emit through normal React reactivity.

**Source-of-truth conflict.** Zustand and the URL are two views of the same state. `useStateFromUrl` writes to Zustand on URL change; `useUrlFromState` writes to the URL on Zustand change. The string-equality emit-guard suppresses the trivial echo case (URL → state → identical URL). The non-trivial case (a Zustand mutation arrives mid-resolve) is handled by the origin-tag in §4.2 — the in-flight resolve is discarded if the URL changed, so it cannot overwrite the user's mutation.

### 4.4 Mounting order and state-driven `staleUrlContext` clearing

```tsx
// App.tsx
useStateFromUrl();    // reads URL → writes Zustand
useUrlFromState();    // reads Zustand → writes URL
```

`useStateFromUrl` is called first so its synchronous setter dispatches inside the effect populate Zustand before `useUrlFromState`'s effect runs on cold mount. After first render, both hooks are reactive and the URL-equality guard makes order irrelevant.

**Clearing `staleUrlContext`.** The setters that own `activePage` / `activeExperimentId` / `activeSampleId` / `activeExposureId` in `state.ts` also clear `staleUrlContext` — implemented inside `state.ts` so it stays a single source of truth and doesn't leak into every call site. Example:

```ts
setActiveSample: (activeSampleId) =>
  set({ activeSampleId, activeExposureId: undefined, staleUrlContext: null }),
```

This cannot live as an externally-attached effect because the CLAUDE.md gotcha forbids reaching into the store via `setState` from outside `state.ts`. It lives at the named-action level so every commit path (NavModal, picker, SSE-driven invalidation) clears the slot uniformly.

## 5. `/` redirect

When `parseLocation` returns `{ kind: "root" }`:

1. Read `activePage`, `activeExperimentId`, `activeSampleId` from persisted Zustand.
2. If `activePage` is missing or unrecognized → `replaceState("/index")`.
3. Otherwise, build the slug URL:
   - Read names from the TanStack cache **synchronously** via `qc.getQueryData(queryKeys.experiments)` and `qc.getQueryData(queryKeys.sample(id))`. Do **not** use `useQuery` — we don't want to wait for the experiments query to finish loading before the URL is set; the redirect must complete before paint, and the entities exist almost certainly.
   - If both names are present in the cache → `replaceState` to the slug URL.
   - If either is missing → fetch `/api/resolve?experiment_id=N&sample_id=M` once. 200 → use the echoed names. 404 → `replaceState("/index")` (the persisted entity no longer exists).

**Why resolve-by-id over `getQueryData(queryKeys.sample(id))`.** The cached sample query exists but isn't guaranteed to have run on cold mount. Falling through to a single resolve-by-id request is one round trip vs. waiting for the experiments+sample queries to fire and fulfill, which would block paint behind two cache misses. The redundancy is intentional: `getQueryData` is the fast path; resolve-by-id is the cold-mount fallback.

## 6. StaleUrlPage (404 component)

**AppShell rendering ladder.** The chrome (AppHeader, TabRocker) renders unconditionally; only the page region swaps based on the ladder. This matches §4.2's "AppHeader and TabRocker stay rendered above [the resolving fallback] without a layout shift."

```tsx
// components/AppShell.tsx (sketch)
const resolving       = useAppState(s => s.resolving);
const staleUrlContext = useAppState(s => s.staleUrlContext);

return (
  <>
    <AppHeader />
    <TabRocker />
    <main>
      {resolving       ? <ResolvingFallback />
       : staleUrlContext ? <StaleUrlPage staleUrlContext={staleUrlContext} />
       : <PageRouter />}
    </main>
  </>
);
```

`resolving` outranks `staleUrlContext` so a navigation away from a stale URL doesn't briefly re-show the old 404 between popstate-fire and the §4.2 pre-fetch clear (defense-in-depth — the pre-fetch clear should already cover it). `staleUrlContext` outranks the page-router because a `kind: "stale"` URL has no `activePage`, so per-page mounting would leave the 404 unreachable.

URL stays at the stale path — the 404 *is* the response for that URL, so the link remains shareable as broken.

**`StaleUrlContext` type:**

```ts
export type StaleUrlContext =
  | { kind: "not_found"; missing: "experiment" | "sample" | "exposure";
      missing_value: string;
      experiment_resolved: { id: number; name: string } | undefined;
      sample_resolved:     { id: number; name: string } | undefined }
  | { kind: "unknown_path"; raw: string };
```

The `not_found` variant is populated from the `/api/resolve` 404 body. The `unknown_path` variant is set by `useStateFromUrl` when `parseLocation` returns `kind: "stale"` (i.e. the URL doesn't match any grammar rule).

```tsx
// components/StaleUrlPage.tsx (sketch)
<div role="alert"
     data-testid="stale-url-page"
     data-missing={dataMissing}>
  <h2>{header}</h2>
  <p>It may have been renamed or removed.</p>
  <button onClick={onPick} data-testid="stale-url-cta">
    {ctaLabel}<kbd>/</kbd>
  </button>
</div>
```

`data-testid="stale-url-page"` and `data-missing` follow the project's E2E selector convention (CLAUDE.md gotcha: never assert on Tailwind class strings).

**Per-variant copy (computed from `staleUrlContext`):**

| variant | `data-missing` | `header` | CTA |
|---|---|---|---|
| `kind: "not_found"`, `missing: "experiment"` | `"experiment"` | `Experiment '{missing_value}' not found.` | "Select an experiment" `/` → `recoverFromStaleUrl({ step: "experiment" })` |
| `kind: "not_found"`, `missing: "sample"` | `"sample"` | `Sample '{missing_value}' not found in '{experiment_resolved.name}'.` | "Select another sample" `/` → `recoverFromStaleUrl({ step: "sample", experimentId: experiment_resolved.id })` |
| `kind: "not_found"`, `missing: "exposure"` | `"exposure"` | `Exposure '{missing_value}' not found in '{sample_resolved.name}'.` | "Back to sample" `/` → `recoverFromStaleUrl({ step: "sample", experimentId: experiment_resolved.id, sampleId: sample_resolved.id, openModal: false })` — the prior sample is still valid, so we just snap back to it; Inspect filmstrip then handles per-exposure picking |
| `kind: "unknown_path"` | `"path"` | `Page not found.` | "Go to Index" → `setActivePage("index")` (clears `staleUrlContext` per §4.4; `useUrlFromState` emits `/index`) |

**`recoverFromStaleUrl(opts)` named action** in `state.ts` — single atomic Zustand transition. Clears `staleUrlContext`, sets active-* ids when provided, and optionally opens NavModal at the requested step. One transition means `useUrlFromState` recomputes once, no half-state where the slot is cleared but the modal hasn't opened (or vice versa):

```ts
type RecoverOpts = {
  step: NavModalStep;
  experimentId: number | undefined;
  sampleId: number | undefined;
  openModal?: boolean;         // default true; row "exposure" passes false
};

recoverFromStaleUrl: (opts: RecoverOpts) =>
  set((s) => ({
    staleUrlContext: null,
    activeExperimentId: opts.experimentId ?? s.activeExperimentId,
    activeSampleId:     opts.sampleId     ?? undefined,
    activeExposureId:   undefined,
    navModalOpen:       opts.openModal ?? true,
    navModalStep:       opts.step,
  }))
```

Notes on field semantics:
- **`experimentId` `?? s.activeExperimentId`** — row 1 omits `experimentId` so the previous active experiment stays set; the URL reads `/index/<previous-exp>` while NavModal is open at experiment step. Intentional: gives the user a visual anchor to where they were before pasting the bad URL.
- **`sampleId ?? undefined`** — row 1 and row 2 omit `sampleId` to clear it (fresh sample selection); row 3 passes the still-valid `sample_resolved.id` to preserve it.
- **`openModal ?? true`** — rows 1, 2 default to opening NavModal; row 3 explicitly passes `false` (user is snapped to the still-valid sample, no modal needed).
- **`navModalStep`** is set even when `openModal: false` (row 3) so any later `/` keypress opens NavModal at the right starting step.

This consolidates rows 1–3 into one action. Without it, row 1 (`openNavModal("experiment")` alone) wouldn't clear `staleUrlContext` — closing the NavModal without picking would leave the user back at the StaleUrlPage. Routing all three through `recoverFromStaleUrl` makes that bug structurally impossible. The `unknown_path` row uses `setActivePage("index")` because there's no NavModal to open — the user is taking themselves to the Index empty state, where NavModal will auto-open at the experiment step via the existing `IndexPage` mount logic.

Notes: the `unknown_path` variant has no `missing_value` and no scope; the simpler "Page not found." copy is intentional. The trailing-period punctuation is uniform across variants. The `unknown_path` CTA dispatches a Zustand action rather than a raw `replaceState` because `replaceState` does not fire `popstate`, so `useStateFromUrl`'s pre-fetch clear would not run and `staleUrlContext` would remain set.

**Keyboard.** `/` triggers the same handler as the button (matches `useGlobalShortcuts` line 35 binding for opening NavModal). The visible `<kbd>/</kbd>` chip teaches the shortcut.

**Slot lifecycle.** `staleUrlContext` is cleared by the `setActive*` actions in `state.ts` (§4.4). Successful navigation (NavModal commit, picker selection, programmatic redirect) automatically clears the slot. Not persisted to localStorage — a stale URL across reloads should re-resolve from scratch, since the underlying 404 might no longer apply.

## 7. SSE-driven URL invalidation

When an SSE event causes the current URL's slug to disappear (another user — or another tab of the same user — deletes the sample you're viewing), `useUrlFromState`'s subscriptions to the `samples` cache see the row removed; the recompute reads `undefined` for the slug; the URL-equality check fails; the hook emits `replaceState` to the closest still-valid URL and `setStaleUrlContext({ ... })` so the user sees the 404 page.

**No `client_id` self-echo filter on this path.** Unlike the `applyRemoteToCache` path (which filters self-echoes per `lib/clientId.ts`), URL invalidation runs on *every* SSE-driven cache change, including the same-user-different-tab case. This is the explicit intent: edits in one tab should refresh the URL state in another tab of the same user, not just other users' edits. CLAUDE.md gotcha: "Two tabs of the same user are treated as distinct subscribers — edits in one tab refresh the other."

In practice the URL-invalidation path fires only on entity deletion (names are stable post-#88 for samples, and unchanged through current flows for experiments and exposures).

## 8. Tests

### 8.1 Vitest (frontend unit)

- `parseLocation`: round-trip every URL shape in §2; edge cases (trailing slash, missing slug, invalid `id`, encoded slugs).
- `useStateFromUrl`: mocked `/api/resolve` for happy path each shape; 404 each missing-entity variant; popstate cancellation; **origin-tag race test** — start a fetch, mutate `location.pathname` mid-flight, assert the in-flight response is discarded.
- `useUrlFromState`: push/replace correctness per the §4.3 table; no emit when URL unchanged; continuity emit on TabRocker; **replay-as-rerun no-spurious-emit** — simulate a rollback → applyRemoteToCache → re-onMutate cycle and assert no extra `pushState` calls.
- `/` redirect: full state, missing pieces, missing entities, `getQueryData` cache hit vs. resolve-by-id fallback.
- `StaleUrlPage`: render variant per `(missing, scope)`; CTA dispatches the right NavModal step; `/` keypress triggers CTA. Cover all four `data-missing` values.
- `state.ts`: each `setActive*` clears `staleUrlContext`.

### 8.2 Playwright (mocked)

- Paste deep URL → land at right page. Assert no flash of wrong content via `data-testid="resolving"` appearing exactly once and the prior page's entity-keyed elements never painting under the new URL.
- Paste stale URL → 404 page; click CTA → NavModal opens at right step. Selectors use `data-testid="stale-url-page"` and `data-missing="<missing>"`.
- Back/forward through pushed states (TabRocker → NavModal commit → thumbnail click) restores each state.
- TabRocker continuity: navigate to a sample, switch tab, switch back → URL is `/<page>/<exp>/<sample>` not bare.

### 8.3 Playwright (live integration)

- Open `/inspect/<exp>/<sample>` in tab A. Delete the sample via reingest from tab B. Tab A receives the SSE; the URL replaces and `<StaleUrlPage>` renders.
- Same-user-different-tab: open same URL in two tabs (same `username`), delete from one, assert the other shows the 404 page (no `client_id` self-echo filter).
- Paste a URL referencing a sample deleted in another session. Land on the 404 page.

### 8.4 Julia (backend)

- `/api/resolve` happy path — name and id forms, plus mixed (`experiment` + `sample_id`).
- `/api/resolve` 404 — each missing variant emits the right `missing` + `missing_value` + partial-resolved fields.
- `/api/resolve` 400 — same-entity name+id collision (`experiment` + `experiment_id`).
- **Tiebreaker pinning** — insert two experiments with the same name, two ids; assert the resolve returns the lower id deterministically. (After future #88-style experiment migration this test should fail and be deleted.)
- **Response-shape contract** — add a row to `test_route_response_shapes.jl` pinning the exact key sets on the 200, 404, and 400 bodies. Mirror with a TS interface in `api.ts`.
- **Stale-name regression** — resolve an experiment by name, rename it via raw SQL, assert old name 404s and new name resolves.
- **Cold-mount Zustand staleness** — resolve-by-id with a `sample_id` whose row was deleted; assert 404 with `missing: "sample"`.
- **`exposures_unique_filename` migration** — direct invocation of `migrate_exposures_unique_filename!(db)` (FK-heal regression-test pattern):
  - Clean DB → index exists, no warnings.
  - Idempotent re-run → no-op.
  - Synthetic-duplicate fixture (insert two rows with identical `(sample_id, filename)` directly via SQL, bypassing the upsert) → second-and-later renamed to `<filename>-2`, …; `@warn` emitted per rename; index created.
  - Transactionality (whole-helper wrapped in `SQLite.transaction(db) do … end`) is verified by code review, not at test time — same as the `migrate_samples_naming!` precedent, where there is no Julia-test idiom for forcing a mid-transaction interruption short of monkey-patching SQLite.jl.
- SPA fallback (uses synthetic `dist_dir` with one HTML + one asset):
  - `/foo` → 200 `index.html`, `Cache-Control: no-store`.
  - `/inspect/exp/sample` → 200 `index.html` (multi-segment must reach the catch-all — pins the `/**` syntax fix).
  - `/api/foo` → 404 (does not fall through).
  - `/asset.png` → served by `dynamicfiles`, normal cache headers (not the no-store HTML path).

## 9. File-level changes

**Backend:**
- New: `packages/HimalayaUI/src/routes_resolve.jl`, `packages/HimalayaUI/test/test_routes_resolve.jl`, `packages/HimalayaUI/test/test_spa_fallback.jl`.
- Edit: `packages/HimalayaUI/src/server.jl` — `register_resolve_routes!()` call + `/**` SPA catch-all.
- Edit: `packages/HimalayaUI/src/HimalayaUI.jl` — `include("routes_resolve.jl")`.
- Edit: `packages/HimalayaUI/src/db.jl` — add named helper `migrate_exposures_unique_filename!(db)` following the `migrate_samples_naming!` precedent (dedupe-then-enforce inside `SQLite.transaction`); call from `migrate_schema!` **after `_fix_fk_references_after_autoincrement_migration!`** (so the AUTOINCREMENT rebuild has settled and the index attaches to the rebuilt `exposures` table — see §3.1 for the rationale).
- Edit: `packages/HimalayaUI/test/test_db.jl` — directly invoke `migrate_exposures_unique_filename!` on synthetic-duplicate fixtures (FK-heal regression-test pattern; bypasses `open_db`'s full migration chain).
- Edit: `packages/HimalayaUI/test/runtests.jl` — register the new test files.
- Edit: `packages/HimalayaUI/test/test_route_response_shapes.jl` — add resolve 200/400/404 shape rows.

**Frontend:**
- New: `lib/url/parseLocation.ts`, `lib/url/parseLocation.test.ts`.
- New: `hooks/useStateFromUrl.ts`, `hooks/useUrlFromState.ts`, plus `.test.tsx` siblings.
- New: `components/StaleUrlPage.tsx`, `components/StaleUrlPage.test.tsx`.
- New: `components/ResolvingFallback.tsx` — the near-empty placeholder rendered while a resolve is in flight (§4.2).
- Edit: `App.tsx` — mount the two hooks.
- Edit: `components/AppShell.tsx` — read `staleUrlContext` and `resolving`; render `<StaleUrlPage>` / `<ResolvingFallback>` / page-router accordingly.
- Edit: `state.ts` — add `staleUrlContext: StaleUrlContext | null` (type defined in §6), `setStaleUrlContext`; add `resolving: boolean`, `setResolving`; add `recoverFromStaleUrl(opts)` and export the `RecoverOpts` type alongside `StaleUrlContext` (signatures in §6) for the not_found:* CTAs; have `setActivePage` / `setActiveExperiment` / `setActiveSample` / `setActiveExposure` clear `staleUrlContext` (do NOT clear `resolving` — it's controlled by `useStateFromUrl` lifecycle). All new slots are **not** persisted (omit from `partialize`). No localStorage version bump (ephemeral slots).
- Edit: `api.ts` — `ResolveSuccess`, `ResolveError404`, `ResolveError400` types using `T | undefined` for optional fields.
- New: `e2e/permalinks.spec.ts` (mocked); `e2e/live/permalinks.spec.ts` (live).

## 10. Out of scope (preserved from issue + reviewer suggestions)

- Slug column on `comparisons` (`/compare/<id>/<title-slug>` Stack Overflow style).
- Mention chip click → navigate. `@sample` / `@exposure` chips currently have no click handler; once URLs exist, they could navigate. Separate UX decision.
- **Experiment-name stabilization (a #88-style migration for experiments)** — UNIQUE constraint, convention enforcement, stable/editable split. Deferred to a future issue. Until then, `/api/resolve` uses the deterministic tiebreaker described in §3.1.
- Per-comparison-list query state (filters, sort, page) in URL.
- Exposure-step picker in NavModal (filmstrip handles intra-sample exposure picking; cross-sample exposure picking would need a new picker).

## 11. Risks and migration

- **Additive schema migration.** The `exposures_unique_filename` UNIQUE INDEX is the only schema change, applied via `migrate_exposures_unique_filename!` (dedupe-then-enforce inside a single `SQLite.transaction`, mirroring the `migrate_samples_naming!` precedent from #88). Idempotent on re-run.
- **Reingest concurrency surface.** `routes_experiments.jl:87` does not wrap `reingest!` in `with_idempotency`; under fast double-fire of `POST /api/experiments/:id/reingest`, a SELECT-then-INSERT race against the new `exposures_unique_filename` (or the existing `samples_unique_name` from #88) can surface as HTTP 500 on a UNIQUE-constraint violation. Operationally rare; pre-exists this spec; flagged as a follow-up to wrap reingest in `with_idempotency`. Operators should expect this as a known transient until the follow-up lands.
- **Browser back-button surprises.** A user with a long edit session might walk back through dozens of pushed states. The push-only-on-commit policy keeps this manageable, but worth noting.
- **Persisted-URL on cold-mount mismatch.** If localStorage Zustand has a stale `activeSampleId` (the sample was deleted while the tab was closed), the cold-mount `/` redirect resolves through `/api/resolve?sample_id=M` and 404s. We handle it by falling through to `/index` — the user lands in the empty state without explicit feedback. Acceptable: reload-after-deletion is rare, and the user reaching the empty state and re-picking a sample is what they'd do anyway.
- **Non-conformant experiment names.** Backfill data may have spaces in `experiments.name`. Round-trip via `encode/decodeURIComponent` handles encoding; the `ORDER BY id ASC LIMIT 1` tiebreaker handles duplicates. The §8.4 stale-name regression test covers the rename case. A future #88-style experiment migration would supplant both these workarounds.
