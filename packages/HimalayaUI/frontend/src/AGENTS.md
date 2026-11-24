# packages/HimalayaUI/frontend/src — Frontend

React 18 + Vite + TypeScript strict + TailwindCSS 4. TanStack Query for server state, Zustand for client state.

## Where to look

| Task | Location | Notes |
|------|----------|-------|
| App entry | `main.tsx` | `index.html → main.tsx`; `StrictMode > ErrorBoundary > QueryClientProvider > BrowserRouter useTransitions={false} > App`; mounts `#app`. Keep `useTransitions={false}`: React Router 7's default transitions let an arrow key pressed right after a click act on the previous sample (`e2e/workspace-focus.spec.ts`, intermittent) |
| App shell | `App.tsx` (`App`) | Composition root: `AppRoutes` + SSE + mutation-queue effects + shell siblings |
| Server state | `queries.ts` | TanStack Query hooks; `authOpts(username, clientId, clientOpId)` helper lives in `lib/authOpts.ts` |
| API layer | `api.ts` | Fetch wrappers; AuthOpts per-call for mutations |
| Client state | `state.ts` | Zustand store — **use named actions only** |
| Phase palette | `phases.ts` | phase → color mapping |
| Tailwind theme | `styles.css` | `@theme { --color-* … }` |
| App shell + routing | `shell/` | `AppShell`, `TopNav`, `AppRoutes`, `ExperimentEntry`, `LinkGone`, `ResolvingFallback`, `OnboardingFlow`, `ShellPalette`, `InfrastructureBanner`. See [shell/AGENTS.md](shell/AGENTS.md) |
| Composites | `components/` | Page-composing components (rails, plates, panels, rows, modals) built from the `ui/` primitives |
| Triage workspace | `workspace/` | The three-pane workspace (`WorkspacePage` at `/experiments/:id/workspace/:sampleId?`, outside `AppShell`; pure `model.ts`). Its components live in `ui/` with every other appearance-owning component; `WorkspacePage` stays placement-only. The plot (`ui/WorkspacePlot.tsx`) draws the approved mock-up directly on `plot/PlotFrame` + `projection.ts`; there is no shared trace-plot engine any more. |
| UI primitives | `ui/` | Closed-look design-system primitives (Button, Card, PhaseChip, ModalShell, Kicker, IconButton, Dot, ToastContainer, …). Appearance lives here; consumer `className` is placement-only. See "Design system" below. |
| Render layers | `{plot,detector}/` | Appearance-authoring render layers (plot frame and projection, detector image). The `lint:design` appearance guard excludes them |
| Pages | `pages/` | `ExperimentsHomePage`, `NewExperimentPage` (under the single `AppShell`) |
| Hooks | `hooks/` | `useFocusTrap`, … |
| Library | `lib/` | URL helpers, plot helpers, peak refs, custom index, labels, bridges |
| Mutation queue | `lib/queue/` | See [lib/queue/AGENTS.md](lib/queue/AGENTS.md) |

## TypeScript strict + `exactOptionalPropertyTypes`

`set({ username: undefined })` fails — optional fields declared as `string | undefined` (rather than `username?: string`) keep this ergonomic. For passing optional values through (e.g. `AuthOpts`), use the `authOpts(username, clientId, clientOpId)` helper in `lib/authOpts.ts`, which omits any undefined key (never `{ username: undefined }`).

## Zustand — named actions

Use the store's named actions (`clearUsername`, `setTutorialSeen`, `openPalette`, …). Avoid `useAppState.setState({ ... })` — direct setState bypasses encapsulation and triggers lint warnings. New state transitions go in `state.ts` as named actions.

## Zustand persist — version-bump wipe-guard

Two invisible invariants in `state.ts` (`persist` config, `queries.ts`'s sibling — see `state.ts:494-530`). Adding a persisted field + naively bumping `version` silently wipes **every existing user's** state on first load:

- **Bumping `version` WITHOUT a `migrate` discards the entire persisted blob.** Zustand v4 drops a blob whose version it can't migrate — all prefs (username/tutorialSeen/active*/…) gone on upgrade. Every version bump must ship a `migrate`. The current migrate's only job across all prior versions is unconditional dead-key strips (`activePage` I5.2, `theme` R0a), so there's no `switch (version)`.
- **`migrate` must return the ORIGINAL object on malformed/non-object input — never `{}`.** Handing `{}` is itself a partial wipe; returning the original lets `merge` fold whatever survived (or fall back to defaults). The wipe-guard is `if (persisted && typeof persisted === "object") { …strip… } return persisted` — the `else` returns untouched.

## State split (load-bearing)

- **Zustand owns *client* state**: the user, tutorial, palette and ingest-progress fields. The URL owns the current experiment, sample and frame.
- **TanStack Query owns *server* state**: experiments, samples, exposures, peaks, indices, assignment.
- Mutations invalidate scoped query keys (`queryKeys.peaks(id)`, `queryKeys.indices(id)`). Don't mix the two concerns in the same hook.

## Per-tab SSE identity

SSE self-echo filtering uses a per-tab `client_id` minted into `sessionStorage` on first load (see `lib/clientId.ts`). Audit identity (`actor` / `X-Username`) is unchanged. **Two tabs of the same user are distinct subscribers** — edits in one refresh the other. `client_id` lives for the tab session: survives reload, scoped to one tab. See `docs/event-log.md` §"Client side".

## Imperative renderers in effects

Wrap any function that is both defined inside a component AND used as a `useEffect` dependency in `useCallback` with its true deps. The effect then depends on `[theCallback]` alone — no redundant dep list, no eslint-disable.

**Async re-entry guards need a ref, not state.** Block double-invocation with a synchronous `inFlight` ref (`if (inFlight.current) return` before `setPending(true)`), *not* the `pending` state alone — state flips asynchronously, so two clicks in one tick both pass a state-only guard. The ref blocks re-entry; the state drives the UI. Same shape for any async-action hook.

## Tailwind v4 theming

"The Print" palette (light warm paper, terracotta accent hue 38, Source Sans 3 throughout) is defined once in `styles.css` via `@theme { --color-* ... }` (R0a, #221). The app has a dark theme selected by `data-theme`, set on `<html>` by `useDocumentTheme` (`lib/theme.ts`, mounted in `App`) so every page follows the pick, and on the workspace root from the same store; `useDocumentTheme` sets `data-theme` on `<html>`: the dark blocks at the end of `styles.css` restate the chrome, status and `--color-phase-*` tokens, and `phaseColor()` returns `var(--color-phase-…)` so phase colour follows the theme. Any token built from another token must be restated there too. The legacy dark-era neutral-ramp shim (the duplicated `bg-`/`text-`/`border-` neutral utilities that mirrored the canonical Print names) was excised in R3-F (#259): use the canonical Print utilities directly — `bg-paper`/`bg-paper-sunk`/`bg-plate`, `text-ink`/`text-ink-soft`/`text-ink-faint`, `border-hair`/`border-hair-strong`, `text-print-accent` (or `accent` for terracotta). Reintroducing an old name is self-revealing: Tailwind won't generate the utility and its `--color-*` custom property no longer resolves. Titles use the `text-display`/`text-headline` roles. To add a new color, add it to `@theme` first.

## Design system — closed-look primitives + the design guard (ENFORCED)

The Print's recurring patterns live as **closed-look** primitives in `ui/` (Button, Card, PhaseChip, ModalShell, Kicker, IconButton, Dot, ToastContainer). They own their appearance via semantic props (`variant` / `size` / `tone` / domain props); a consumer's `className` is **placement-only** (margins, position, grid). To change how a primitive looks, build a variant *into the primitive* — the idiom is a `Record<Variant,string>` map + a tiny local `cx()` join helper (no cva/clsx/tailwind-merge). Don't restyle from the outside.

**Ratchets (2026-11-19).** `npm run lint:ratchet` (`scripts/check-ratchets.mjs`, also in `npm run build`) counts two things per file against `scripts/ratchet-baseline.json` and fails if any count grows:

- **Appearance tokens** (counted per token, not per `className`) in every string or template literal inside a consumer `className=` expression, including `cx(...)` and ternaries. Named tokens count too (`text-ink-soft`, `text-caption`, `font-*`, `bg-*`, `border*`, `rounded*`, …); text alignment and wrapping are placement. Exempt: `{ui,plot,detector}/`, `ui/` and stories. `lint:design` catches bracket values and raw colours; the ratchet stops the named-token count from growing. It does not see styling built outside a `className=` expression.
- **Comment lines** in every file `git ls-files` lists with a Julia, TypeScript, JavaScript or CSS extension, minus the short `EXCLUDED` list in the script. That covers every package, tests, e2e, scripts, stories, config files and anything added later; a test fails if a tracked source file is neither measured nor excluded. TypeScript is read through its AST, so trailing and JSX comments count and strings do not. Julia counts `#`, `#= =#` and docstrings in every form (`"""…"""`, a one-line `"…"` above a definition, `@doc` and `raw"""`); SQL strings, a string a block returns, and `#` inside interpolation or command literals do not count. CSS counts `/* */`. Agents don't write comments: put the reason in a name, a small function, a type, a test or an assertion, because comments go stale and code doesn't (see CLAUDE.md "Conventions"). Exempt: comments that *start* with `physics:`, and one-line `eslint-disable… -- <reason>`, `@ts-expect-error` and `/// <reference>` directives of at most 120 characters. A reason-less directive, a longer one and `@ts-ignore` count.

A count that falls must be locked in with `npm run lint:ratchet -- --update` (commit the baseline). `--update` refuses any per-file increase, so allowance can't move between files: a new file or a split starts at zero. The one exception is a pure rename, where a new file has exactly the count of a baseline file that no longer exists.

`npm run lint` runs ESLint (`react-hooks/rules-of-hooks`, `react-hooks/exhaustive-deps`, `no-console`) on all of `src/`, stories included, plus `copy/ui-copy` outside stories (see "UI copy").

This is **mechanically enforced** (2026-05-29 extraction). `scripts/check-design.mjs` runs as a pure-absolute `lint:design` step prepended to `npm run build` (plus a warn-only PostToolUse hook), and **fails the build** on any inline appearance utility *outside* `ui/**`, `plot/**` and `detector/**`:

- arbitrary type size `text-[…]` → use a named scale role (`text-xs/sm/base/lg/xl/headline-lg/display`)
- arbitrary radius `rounded-[…]` → `rounded-sm` / `rounded-md` (both 5px) / `rounded-full`
- raw colour literal (`oklch(` / `rgba(` / quoted `#hex`) → a `--color-*` token utility
- side-stripe `border-l/r` > 1px → a full border + a leading icon/word instead

Only the colour-AUTHORING files are exempt (rules #3/#5 share an allowlist: `phases.ts`, `main.tsx`, plus the excluded `{ui,plot,detector}/` prefixes).

## UI copy (ENFORCED)

Copy follows the Microsoft, GitLab Pajamas and Polaris guides on the points below. `copy/ui-copy` (`scripts/eslint-ui-copy.mjs`, run by `npm run lint`) checks the mechanical rules, and its diagnostics carry the fix. `test/copy-no-emdash.test.ts` bans em dashes.

- **Hints are full imperative sentences.** Say what to do, not what the control does: "Click to hide.", "⌘-click to subtract from the current trace.", "Press Enter to add it or Esc to close." Never the telegraphic "Click hides." or "Enter removes it" (lint).
- **Sentences end with a period; names don't.** A tooltip that only names a control ("Keyboard shortcuts", "All experiments") is a short phrase with no period. A string that holds a sentence, or starts with an instruction ("Click to …", "Press …"), ends with one (lint). Placeholders are prompts with no period ("Type to filter").
- **Sentence case everywhere**: buttons, headings, menu items, tooltips. Capitalise only the first word and proper names (lint flags Title Case of three or more words).
- **House vocabulary** (lint, user-facing strings only): "image", not "frame" or "exposure" ("exposure time" stays); "label", not "tag"; "Drop"/"Restore", not "mute". Prefer the plain literal verb over a metaphor.
- **Tooltips add, never hold the only copy.** Don't repeat the visible label, and keep anything important in the visible UI. Tooltips never hold links or buttons; popovers do.
- **No em dashes.** Use a period, a comma or a colon. The standalone "—" missing-value placeholder is the one exception.
- **Name things the way the person sees them**: the sample's name, the setup's folder name, the keycap glyph (⌘ on a Mac, Ctrl elsewhere, via `MOD_KEY` from `lib/modKey`).

The rule reads JSX text, the copy props (`title`, `aria-label`, `placeholder`, `label`, `content`, `hint`, …), the copy keys of plain objects (`what`, `hint`, `label`, `message`, …) and any other string shaped like prose (a capitalised word followed by a lowercase one). Identifiers, paths, query keys and `className` are never read. A general prose linter (write-good, alex) was tried on every user-facing string on 2026-11-22: 17 flags, none a real problem, and it missed the telegraphic form, so it is not used.

## Mutation queue — load-bearing one-liners

Full architecture in `docs/mutation-queue.md`; queue internals in `lib/queue/AGENTS.md`. The invariants that bite UI code outside the queue:

- **Optimistic placeholder ids are NEGATIVE.** `Peak.id < 0` means "not yet confirmed by server"; SSE confirmation overwrites with the positive server id. UI code that filters or compares peak ids must handle negatives.
- **A peak is `(source, id)`.** Auto and manual peak ids come from separate sequences and collide on real data. Compare peaks with `samePeak` and key them with `peakKey` (`lib/peakRef.ts`); plot callbacks pass a `PeakRef`, never a bare id. Test fixtures with both sources should come from `peakList` (`test/fixtures/peaks.ts`), which numbers each source from 1 like the DB.

## Query layer gotchas

- **`queryKeys` namespace isolation — a single-entity key must NOT share a collection prefix.** Collection keys live under `["exposure", id, …]` (`"peaks"`, `"indices"`, `"assignment"`, `"trace"`), and an `invalidateQueries(["exposure", id])` prefix-matches all of them. There are no single-entity keys today; if one comes back, give it its own root (e.g. `["exposure-entity", id]`) so a collection invalidation cannot blast it, and add a reader in the same PR.


## Multi-layer contract testing

Every reconciliation contract has six layers (route emit → SSE payload → `applyRemoteToCache` merge → cache row → `onMutate` → `onSuccess`). When fixing a bug at one layer, add a regression row at every other layer where the same class can manifest. See `docs/contract-testing.md` for canonical paired test files (`cache-shape.test.ts`, `sseEventPayload.contract.test.ts`, `rollbackSymmetry.test.ts`, `authHeaders.test.ts`, `test_route_response_shapes.jl`, `test_idempotency_replay_invariant.jl`).

**Layout measurements go in primitive state.** A hook that measures the DOM in a layout effect and sets state (scroll fades, fit counts) keeps each value in its own `useState<boolean | number>` and sets it directly. Don't merge them into one object with an updater that returns `s` when nothing changed. When a scroll or resize update is still waiting at a lower priority, React replays queued updaters from an older base state on every urgent render, so the updater hands back an older object, React sees a change, commits, the layout effect measures again, and the loop ends in "Maximum update depth exceeded". `Overlay`'s scroll fade did this when Home was pressed during End's smooth scroll (`test/workspace/Overlay.test.tsx` reproduces it with `startTransition` plus `flushSync`).

## Anti-patterns

- Mint `client_op_id` inside `mutationFn`, not at hook creation time.
- Don't assert on Tailwind class strings in tests — use `data-testid` / `data-*` attributes.
- Don't inline appearance utilities (`text-[…]`, `rounded-[…]`, raw colours, side-stripes) in a consumer — `lint:design` fails the build. Put appearance in a `ui/` primitive; `className` is placement-only (see "Design system").

## Dependency pins

`package.json` `overrides` pins `@joshwooding/vite-plugin-react-docgen-typescript` to `^0.9.1` (no vulnerable glob). Storybook 10.6 still asks for `^0.7`, and the plugin only loads with `reactDocgen: "react-docgen-typescript"`, which `.storybook/main.ts` does not set. Drop the override once `@storybook/react-vite` depends on `^0.9` or later.
