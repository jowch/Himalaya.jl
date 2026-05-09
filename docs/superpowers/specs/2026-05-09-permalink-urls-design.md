# Permalink URLs — Design Spec

**Issue:** #89 (depends on #88, which has landed at `8ac2bf6`)
**Status:** Drafted 2026-05-09 (pre-review)
**Worktree:** `.claude/worktrees/permalinks` (branch `permalinks`)

## 1. Problem

HimalayaUI is single-page (`/`) with all selection state in Zustand + localStorage. There are no addressable URLs. A user can't share "look at this sample" / "look at this exposure" / "look at this comparison" by copying the address bar. Bookmarks are useless. Refreshing the tab keeps you where you were — but only because Zustand persists; nothing about the URL communicates intent.

This spec adds slug-based permalinks for every page in the app. A user copies the address bar; the recipient pastes; the recipient lands at the same screen.

#88 is a hard prerequisite: pre-#88, `samples.name` was UI-mutable and not unique within an experiment, so any URL keyed on it would break on rename and collide on duplicate names. Post-#88, `samples.name` is stable, convention-enforced (`^[A-Za-z0-9._-]+$`), and unique within experiment. Slug-safety is therefore a property of the data, not a runtime escape-hatch.

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
- `<experiment>` — `experiments.name`. Already free-form today; legacy DBs may contain non-conformant names. Conformant data emits verbatim; non-conformant percent-encodes at the boundary.
- `<sample>` — `samples.name`, `^[A-Za-z0-9._-]+$` post-#88, unique within experiment.
- `<filename>` — `exposures.filename`. Note: the column stores the **stem** (the `{name}` substitution), not the full on-disk filename — `resolve_files` strips the integration_pattern suffix at ingest (`config.jl:219`). So `exposures.filename = "JC001-007"`, not `"JC001-007_int.dat"`. This is what manifests speak in and what operators expect to see in URLs. Unique per sample by upsert contract (`cli.jl:191`).
- `<id>` — `comparisons.id`. Numeric. A future enhancement could add a slug column on comparisons (`/compare/<id>/<title-slug>` Stack Overflow style); out of scope here.

**Encoding.** Conformant slugs (matching `[A-Za-z0-9._-]+`) emit and parse verbatim — no percent-encoding either way. The `parseLocation` parser uses `decodeURIComponent` on every captured slug to handle the legacy non-conformant case.

## 3. Backend

### 3.1 `GET /api/resolve`

New file `packages/HimalayaUI/src/routes_resolve.jl`. Registered in `register_routes!` between `register_picker_routes!()` and the SPA catch-all (so `/api/*` precedence is unambiguous).

```
GET /api/resolve?experiment=<name>[&sample=<name>][&exposure=<filename>]
→ 200 { experiment_id, sample_id?, exposure_id? }
→ 404 {
    error: "not_found",
    missing: "experiment" | "sample" | "exposure",
    missing_value: "<value>",
    experiment_resolved?: { id, name },     // present when chain partially resolved
    sample_resolved?:     { id, name }      // present when missing == "exposure"
  }
```

**Resolution order.** Walk the chain left-to-right. First miss is the 404. Earlier-resolved entities ride along in `experiment_resolved` / `sample_resolved` (only the relevant prefix is included — e.g. missing `"sample"` includes `experiment_resolved` only; missing `"exposure"` includes both).

**Read-only.** No writes; no `with_idempotency` wrapping; no events emitted. Just three SELECTs.

**Why a single endpoint instead of three.** A single round trip for the common `/index/<exp>/<sample>` case. Three calls would race or require client-side chaining. The 404 body shape gives the frontend everything it needs to position the recovery UX in one shot.

### 3.2 SPA catch-all

Add to `server.jl::register_routes!` immediately after `Oxygen.dynamicfiles(dist_dir, "/")` (only when `isdir(dist_dir)`):

```julia
@get "/{rest:.*}" function(req::HTTP.Request, rest::String)
    startswith(rest, "api/") && return HTTP.Response(404, "Not found")
    file = joinpath(dist_dir, rest)
    return isfile(file) ?
        Oxygen.file(file) :
        HTTP.Response(200,
            ["Content-Type" => "text/html", "Cache-Control" => "no-store"],
            read(joinpath(dist_dir, "index.html")))
end
```

- The `api/` guard is critical — without it, an unregistered API path falls through to `index.html`, masking 404s as 200s and breaking client error handling.
- The `isfile(file)` branch lets static asset URLs (`/foo.png`) keep their normal cache headers via `Oxygen.file`.
- `Cache-Control: no-store` on the HTML shell ensures users get the latest after deploy. Vite-bundled JS/CSS keep their content-hashed long-cache headers.

**Test-environment guard.** The catch-all only mounts when `dist_dir` exists. The Julia test harness has no frontend dist; routes_resolve is exercised directly via HTTP and the SPA fallback never runs in tests.

## 4. Frontend URL-sync layer

No router library. Two hooks plus a pure parser, mounted at `App.tsx` under the `QueryClient`. Total surface: ~150 lines across the parser + the two hooks + the page component.

### 4.1 `parseLocation(pathname, search) → ParsedUrl`

Pure function in `packages/HimalayaUI/frontend/src/lib/url/parseLocation.ts`. No side effects.

```ts
export type ParsedUrl =
  | { kind: "root" }
  | { kind: "index";   experiment?: string; sample?: string }
  | { kind: "inspect"; experiment?: string; sample?: string; exposure?: string }
  | { kind: "compare"; view: "list" }
  | { kind: "compare"; view: "new" }
  | { kind: "compare"; view: "review";  id: number }
  | { kind: "compare"; view: "edit";    id: number }
  | { kind: "stale";   raw: string };          // unrecognized path

export function parseLocation(pathname: string, search: string): ParsedUrl;
```

Rules:
- Empty path or `/` → `{ kind: "root" }`.
- `/<page>` where page ∈ {index, inspect, compare}: split on `/`, `decodeURIComponent` each segment, populate fields per the grammar.
- `/compare/<n>` requires `n` to be a positive integer; otherwise `{ kind: "stale" }`.
- Anything else → `{ kind: "stale", raw: pathname + search }`.
- For `inspect`, the `?exposure=` query param is read only when both experiment and sample are present.

### 4.2 `useStateFromUrl()`

Mounted in `App.tsx`. On mount and on every `popstate`:

1. `parseLocation(location.pathname, location.search)`.
2. If `kind: "root"` → run the §5 redirect.
3. If `kind: "stale"` → set Zustand `staleUrlContext = { kind: "unknown_path", raw }`.
4. Otherwise → `fetch("/api/resolve?…")` with whichever slugs are present.
   - 200 → dispatch Zustand setters: `setActivePage`, `setActiveExperiment`, `setActiveSample` (or undefined, depending on what resolved), `setActiveExposure`. Clear `staleUrlContext`.
   - 404 → set `staleUrlContext = body` (the full 404 payload). The page renders `<StaleUrlPage>` per §6.

**Race avoidance.** A `popstate` arriving mid-resolve (e.g. user holds the back button) cancels the in-flight fetch via `AbortController`. The latest `popstate` always wins.

**No automatic redirect on partial 404.** The earlier "redirect to fallback" model (per the original issue text) is replaced by the 404 page; the URL stays at the broken value until the user picks something.

### 4.3 `useUrlFromState()`

Mounted in `App.tsx`. Subscribes to `(activePage, activeExperimentId, activeSampleId, activeExposureId)` and to the relevant TanStack queries (`experiments`, `samples` for the active experiment) so it can look slugs up by id without a fetch.

On any change, computes the target URL and emits `pushState` or `replaceState` per the policy:

| Trigger | Push or replace |
|---|---|
| TabRocker click → `activePage` change | push, target `/<page>/<lastExp>/<lastSample>` (continuity) |
| NavModal commit (closes with full selection) | push |
| Intermediate NavModal state | no URL change (modal isn't committed) |
| Inspect thumbnail click → `activeExposureId` change | replace |
| App-init URL sync (initial render after `useStateFromUrl`) | replace |
| SSE-driven cache update that invalidates the current URL | replace, then re-resolve |

The hook **does not emit** when the resulting URL equals the current one (checked via string compare on `pathname + search`). Prevents spurious history entries when, for example, an SSE refetch hydrates names without changing them.

**Continuity rule.** When the user TabRockers from `/compare/123` to Index, we read `activeExperimentId` / `activeSampleId` from Zustand and emit `/index/<exp>/<sample>` if both are present. If either is missing, emit the bare page (`/index`).

**Source of truth conflict.** Zustand and the URL are two views of the same state. `useStateFromUrl` writes to Zustand on URL change; `useUrlFromState` writes to the URL on Zustand change. To prevent feedback loops, `useUrlFromState` no-ops when the computed URL equals the current URL. There is no "lock" or sequence number — the equality check is sufficient because the URL→state and state→URL paths each canonicalize before comparing.

### 4.4 Mounting order

```tsx
// App.tsx
useStateFromUrl();    // reads URL → writes Zustand
useUrlFromState();    // reads Zustand → writes URL
```

`useStateFromUrl` runs first so the initial Zustand read in `useUrlFromState` reflects the URL on cold mount. After that, both hooks are reactive and the order doesn't matter.

## 5. `/` redirect

When `parseLocation` returns `{ kind: "root" }`:

1. Read `activePage`, `activeExperimentId`, `activeSampleId` from persisted Zustand.
2. If `activePage` is missing or unrecognized → `replaceState("/index")`.
3. Otherwise, fetch the names: prefer the cached `experiments` query, fall back to a single round-trip lookup via `/api/resolve?experiment_id=N&sample_id=M` (extension to the resolve endpoint — see below). On 404, fall through to `/index`.
4. `replaceState` to the slug URL.

**Resolve-by-id extension.** Add `experiment_id` / `sample_id` query params to `/api/resolve` for the cold-mount path. Mutually exclusive with name params; the 200 response shape is identical (`{ experiment_id, sample_id?, exposure_id? }`) plus `experiment_name` / `sample_name` / `exposure_filename` so the redirect can build the URL without a follow-up fetch. Keeps the cold-mount path to a single round trip.

## 6. StaleUrlPage (404 component)

Mounted in the page slot under `AppShell` when Zustand `staleUrlContext` is non-null. URL stays at the stale path — the 404 *is* the response for that URL, so the link remains shareable as broken.

```tsx
// components/StaleUrlPage.tsx (sketch)
<div role="alert" className="…">
  <h2><EntityType> '<value>' not found<scope>.</h2>
  <p>It may have been renamed or removed.</p>
  <button onClick={onPick}>
    Select another <entity>
    <kbd>/</kbd>
  </button>
</div>
```

**Variant by `staleUrlContext.missing`:**

| missing | scope | CTA |
|---|---|---|
| `"experiment"` | `""` | "Select an experiment" → `openNavModal("experiment")` |
| `"sample"` | `" in '<experiment>'"` (from `experiment_resolved.name`) | "Select another sample" → preselect experiment via `setActiveExperiment(experiment_resolved.id)`, then `openNavModal("sample")` |
| `"exposure"` | `" in '<sample>'"` (from `sample_resolved.name`) | "Select another sample" → same as above; the Inspect filmstrip handles per-exposure picking once the user is on a valid sample |
| `"unknown_path"` | _(no scope)_ | "Go to Index" → `replaceState("/index")` |

**Keyboard.** `/` triggers the same handler as the button (matches the existing `useGlobalShortcuts` binding for opening NavModal). The visible `<kbd>/</kbd>` chip teaches the shortcut.

**Clearing the slot.** `staleUrlContext` is cleared by:
- A successful `useStateFromUrl` round-trip on a subsequent navigation.
- Any direct dispatch of `setActivePage` / `setActiveExperiment` / `setActiveSample` (NavModal commit, picker selection).

The slot is **not** persisted to localStorage — a stale URL across reloads should re-resolve from scratch (the 404 might no longer apply).

## 7. SSE-driven URL invalidation

When an SSE event causes the current URL's slug to disappear (e.g. another user deletes the sample you're viewing), `useUrlFromState` emits `replaceState` to the closest still-valid URL plus sets `staleUrlContext` so the user sees the 404 page. In practice this fires only on entity deletion — names are immutable post-#88 (samples) or unchanged through the relevant flows (experiments, exposures).

The detection is implicit: when `activeSampleId` is non-null but the corresponding row is missing from the `samples` cache after an SSE-driven invalidation, the recompute in `useUrlFromState` reads `undefined` for the slug, sees a URL mismatch with the current location, and falls through to the stale path.

Live-integration test (§8) covers the rename / delete cases.

## 8. Tests

### 8.1 Vitest (frontend unit)

- `parseLocation`: round-trip every URL shape in §2; edge cases (trailing slash, missing slug, invalid `id`, encoded slugs).
- `useStateFromUrl`: mocked `/api/resolve` for happy path each shape; 404 each missing-entity variant; popstate cancellation.
- `useUrlFromState`: push/replace correctness per the §4.3 table; no emit when URL unchanged; continuity emit on TabRocker.
- `/` redirect: full state, missing pieces, missing entities, cache hits vs. fallback.
- `StaleUrlPage`: render variant per `(missing, scope)`; CTA dispatches the right NavModal step; `/` keypress triggers CTA.

### 8.2 Playwright (mocked)

- Paste deep URL → land at right page, no flash of wrong content.
- Paste stale URL → 404 page; click CTA → NavModal opens at right step.
- Back/forward through pushed states (TabRocker → NavModal commit → thumbnail click) restores each state.
- TabRocker continuity: navigate to a sample, switch tab, switch back → URL is `/<page>/<exp>/<sample>` not bare.

### 8.3 Playwright (live integration)

- Open `/inspect/<exp>/<sample>` in tab A. Delete the sample via reingest from tab B (or directly via DB). Tab A receives the SSE; the URL replaces to `/inspect/<exp>` and `<StaleUrlPage>` renders.
- Paste a URL referencing a sample deleted in another session. Land on the 404 page. (Requires the test harness's no-mock backend.)

### 8.4 Julia (backend)

- `/api/resolve`:
  - Happy path: experiment-only, experiment+sample, experiment+sample+exposure.
  - 404: each missing variant emits the right `missing` + `missing_value` + partial-resolved fields.
  - Resolve-by-id: `experiment_id=N&sample_id=M` returns names alongside ids.
- SPA fallback:
  - `/foo` → 200 `index.html`.
  - `/api/foo` → 404 (does not fall through).
  - `/dist-asset.png` (when present) → 200 with normal cache headers (not the no-store HTML path).

## 9. File-level changes

**Backend:**
- New: `packages/HimalayaUI/src/routes_resolve.jl`, `packages/HimalayaUI/test/test_routes_resolve.jl`.
- Edit: `packages/HimalayaUI/src/server.jl` — `register_resolve_routes!()` call + SPA catch-all.
- Edit: `packages/HimalayaUI/src/HimalayaUI.jl` — `include("routes_resolve.jl")`.
- Edit: `packages/HimalayaUI/test/runtests.jl` — register the new test file.

**Frontend:**
- New: `lib/url/parseLocation.ts`, `lib/url/parseLocation.test.ts`.
- New: `hooks/useStateFromUrl.ts`, `hooks/useUrlFromState.ts`, plus their `.test.tsx` siblings.
- New: `components/StaleUrlPage.tsx`, `components/StaleUrlPage.test.tsx`.
- Edit: `App.tsx` — mount the two hooks and the page.
- Edit: `state.ts` — add `staleUrlContext` slot + `setStaleUrlContext` action; add `setActiveExposure` if missing (already exists per §state.ts:236).
- Edit: `pages/IndexPage.tsx`, `pages/InspectPage.tsx`, `pages/ComparePage.tsx` — render `<StaleUrlPage>` when `staleUrlContext` is non-null (or do this at `AppShell` level — decided at write time).
- New: `e2e/permalinks.spec.ts` (mocked); `e2e/live/permalinks.spec.ts` (live).

## 10. Out of scope (preserved from issue)

- Slug column on `comparisons` (`/compare/<id>/<title-slug>` Stack Overflow style).
- Mention chip click → navigate. `@sample` / `@exposure` chips currently have no click handler; once URLs exist, they could navigate. Separate UX decision.
- Experiment-name editable/stable split (parallel to #88's sample model). Experiments are renamed less often; deferred.
- Per-comparison-list query state (filters, sort, page) in URL.
- Exposure-step picker in NavModal (filmstrip handles intra-sample exposure picking; cross-sample exposure picking would need a new picker).

## 11. Risks and migration

- **No DB migration.** All required columns exist post-#88. The resolve endpoint is read-only.
- **Browser back-button surprises.** A user with a long edit session might walk back through dozens of pushed states. The push-only-on-commit policy keeps this manageable, but worth noting.
- **Persisted-URL on cold-mount mismatch.** If localStorage Zustand has a stale `activeSampleId` (the sample was deleted while the tab was closed), the cold-mount `/` redirect resolves through `/api/resolve?sample_id=M` and 404s. We handle it by falling through to `/index` — but this is the one case where the user lands in the empty state without explicit feedback. Acceptable: reload-after-deletion is rare, and the user reaching the empty state and re-picking a sample is what they'd do anyway.
- **Test gap on legacy non-conformant names.** Backfill data may have spaces in `experiments.name`. The percent-encoding boundary should be tested with at least one such case to prove round-trip safety.
