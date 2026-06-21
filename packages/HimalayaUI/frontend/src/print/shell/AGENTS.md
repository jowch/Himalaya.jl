# frontend/src/print/shell — app shell, routing, modals

The live application shell: the single layout shell, the topbar chrome, the route table, and the global onboarding/nav/error overlays. Everything else in the frontend is composed under `src/print/` (`ui/` primitives, `components/` composites, `pages/`, the `plot`/`detector`/`comb`/`waterfall` render layers, `export/` figure utilities) — the app shell itself is `src/print/App.tsx` (`PrintApp`) — plus the non-presentational layers at `src/` root (`state.ts`, `queries.ts`, `api.ts`, `lib/`, `hooks/`).

## Where things live

- **Shell + chrome**: `CorpusShell` (the single layout shell) wraps `CorpusTopbar` (wordmark + stage-tabs + Beamtime chip — composed from `../ui` primitives: `TopBar`, `Wordmark`). No appearance utilities live here — chrome appearance is owned by the `ui/` primitives (closed-look / open-placement contract).
- **Routing**: `AppRoutes` is the route table (samples / loupe / focus / series folio·builder·scoping pages, all imported from `../pages`). `IndexSlugRedirect` and `StaleUrlPage` resolve/repair permalink slugs; `ResolvingFallback` is the in-flight placeholder.
- **Onboarding / modals / errors**: `OnboardingFlow`, `NavModal`, `InfrastructureBanner`. The modals compose `ModalShell` (+ `Button`/`Kicker`/`IconButton`) from `../ui`.

## Import boundary

These files live under `src/print/`, so the `lint:design` import-boundary guard (`scripts/check-design.mjs`, `scanLegacyImports`) forbids them from importing a top-level `src/components/**` or `src/pages/**` — both of which are now retired. Reach for state/data via `../../state`, `../../queries`, `../../api`, `../../hooks`, `../../lib`; reach for presentation via the sibling `../ui` primitives and `../pages` route components.

## Focus trapping in modals

`../../hooks/useFocusTrap.ts` exports `useFocusTrap(containerRef, active)` (APG modal-dialog pattern). On activation it moves focus INTO the container (first focusable, or the container via tabIndex=-1) unless the consumer already placed it inside — `autoFocus` runs before the trap effect and is respected; a parent's own focus effect runs after it (child-before-parent order) and wins, which is how `NavModal`'s input takes focus. It intercepts Tab/Shift+Tab at the **document** level (background controls stay unreachable even if focus escapes), skips Tabs a consumer already `preventDefault`ed (NavModal's Tab-commit), pulls escaped focus back via a `focusin` guard, and restores the previously-focused element on cleanup. `ModalShell` wires it for every modal; `NavModal` and `OnboardingFlow` get it through that.

## Keyboard shortcut registry

`shortcuts.ts` is the single source of truth for all gesture vocabulary. `useShortcuts(bindings)`, `<KbdLegend group=…>`, and the `aria-keyshortcuts` on matching buttons all derive from the registry, so the live handler, the on-screen legend, and the a11y annotation can't drift apart. Load-bearing conventions:

- **DECLINE convention.** A `useShortcuts` handler that returns `false` *declines* the key: the event is left un-`preventDefault`ed so the next listener can act (e.g. `TracePlate`'s own Escape-to-disarm). Returning `undefined`/`void` *claims* the key and triggers `preventDefault`. (`useShortcuts.ts`: `if (handler(e) !== false) e.preventDefault()`.)
- **Two-axis model (rev 2, app-wide — supersedes the 2026-06-13 `[`/`]` lock).** `↑`/`↓` = sample (`prevSample`/`nextSample`); `←`/`→` = detail (`prevDetail`/`nextDetail`) — candidate *preview* on Focus, frame on Loupe. On Focus the detail axis drives page-level `previewIndexId`/`previewedCandidate` state, **not** DOM focus, so the keyboard mirrors mouse-hover on the comb/rings. Exposure stepping is **filmstrip-only** — click `ThumbnailGallery`'s `onSelect`, there is no exposure key (`prevExposure`/`nextExposure` were removed from the registry).
- **`useReorderShortcuts`** resolves the focused row via the `data-reorder-index` attribute on the closest `rowSelector` ancestor, then calls the surface's existing move callback; when focus is not inside a reorderable row it returns `false` (declines), preserving Alt+Arrow's native meaning.
- **Builder undo/redo is snapshot-based** (structural edits — reorder/add/remove; the title field keeps native browser undo). Its `useShortcuts` call must live in the top hooks block, above every early return, or React throws "rendered fewer hooks."
- **Scoping undo is a raw keydown effect binding `⌘Z` only — `⌘⇧Z` redo is NOT wired in Scoping** (Builder has both). This cross-surface inconsistency is known and intentional (out of scope for the shortcut-library pass).

## The floatingDock lane protocol

The `InfrastructureBanner` is fixed bottom-centre by default (most prominent for a "Saving…" / "Couldn't save" status). A page that mounts an opaque, higher-z fixed bottom-centre action bar — the contact sheet's `CullBar` / `ComposeBar` (`bg-ink z-50`) — paints straight over the banner (`z-40`) and occludes it entirely (the banner still renders, it's just hidden underneath). LA-COLLIDE.

Any page mounting such a bar MUST claim the lane: call `setCenterLaneOccupied(true)` on mount and release it (`false`) on unmount. The banner reads the flag and steps aside to the bottom-right corner (free — toasts own top-right) until the lane clears, so both stay visible (`floatingDock.ts:1-30`; publisher e.g. `../pages/SamplesPage.tsx:232-242`; reader `./InfrastructureBanner.tsx:18-23`). This is **pure convention** — a small dedicated Zustand store (`useFloatingDock`), a plain boolean (not a pixel offset, so no brittle coupling to bar heights), with **no type-system enforcement**. Forget the call and the bar silently occludes the banner.

## Key routes on `:id` when they hold ephemeral per-id local state

React reuses the same element instance across param changes on a single route, so page-level LOCAL state (undo/redo stacks, view offset/scale, Confirm-chain refs) bleeds from one id into the next unless the route element is keyed on the param. `SeriesBuilderRoute` wraps `SeriesBuilderPage` with `key={id}` so each series mounts a fresh instance (BU-MODESWITCH-LEAK; `AppRoutes.tsx:58-73`). Any new `/foo/:id` route with ephemeral per-visit local state needs the same wrapper. (The global draft slot is intentionally recoverable on return; the key resets only the per-visit state that should never cross an id boundary.)

## OnboardingFlow defers `setUser` until `closeTutorial`

For a brand-new user who hasn't seen the tutorial, `onSubmitName` stashes the created user in local `pending` state and switches to the tutorial phase **instead of** calling `setUser()` immediately (`OnboardingFlow.tsx:118-124`). The reason: `setUser()` makes `username !== undefined`, which trips the `if (username !== undefined && phase !== "tutorial") return null` unmount guard (`:96`) and skips the tutorial. The pending user is committed only in `closeTutorial()` (`:130-135`, behind `setTutorialSeen(true)`). Calling `setUser()` at creation time silently breaks the tutorial for new users.

## Loupe stepper's `location.state.sampleOrder` contract

The loupe prev/next stepper (`CorpusTopbar`, the `loupeStepper`) derives its traversal order from `location.state.sampleOrder` — the serialized contact-sheet order, passed as React Router history state at navigate time (`CorpusTopbar.tsx:237-245`). When absent (a direct URL visit, or an entry-point that forgets to pass it), `resolveSampleOrder` falls back to a fresh sort, and the stepper walks a *different* order than the sheet the user came from. Any code navigating to the loupe must pass `navigate(url, { state: { sampleOrder: ordered } })` (the stepper's own `goto` already does, threading the order forward step-to-step).

## Anti-patterns

- Don't introspect tooltip text or Tailwind classes in tests — use `data-*` attributes. Tests stay behavioral.
- Don't add inline appearance utilities (`text-[…]`, `rounded-[…]`, raw colors, side-stripes) here — they fail the `lint:design` guard. Appearance lives in `../ui` primitives; a consumer's `className` is placement-only. Full contract: `../../AGENTS.md` §"Design system"; visual reference `docs/design-system.html`; system `DESIGN.md` (root).
