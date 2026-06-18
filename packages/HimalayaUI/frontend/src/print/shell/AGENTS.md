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
- **Focus two-axis model.** On Focus: `[`/`]` = sample, `←`/`→` = exposure (over the `acceptableExposures()` pool only, **not** the full list), `↑`/`↓` = candidate *preview* — these drive page-level `previewIndexId`/`previewedCandidate` state, **not** DOM focus, so the keyboard mirrors mouse-hover on the comb/rings.
- **`useReorderShortcuts`** resolves the focused row via the `data-reorder-index` attribute on the closest `rowSelector` ancestor, then calls the surface's existing move callback; when focus is not inside a reorderable row it returns `false` (declines), preserving Alt+Arrow's native meaning.
- **Builder undo/redo is snapshot-based** (structural edits — reorder/add/remove; the title field keeps native browser undo). Its `useShortcuts` call must live in the top hooks block, above every early return, or React throws "rendered fewer hooks."
- **Scoping undo is a raw keydown effect binding `⌘Z` only — `⌘⇧Z` redo is NOT wired in Scoping** (Builder has both). This cross-surface inconsistency is known and intentional (out of scope for the shortcut-library pass).

## Anti-patterns

- Don't introspect tooltip text or Tailwind classes in tests — use `data-*` attributes. Tests stay behavioral.
- Don't add inline appearance utilities (`text-[…]`, `rounded-[…]`, raw colors, side-stripes) here — they fail the `lint:design` guard. Appearance lives in `../ui` primitives; a consumer's `className` is placement-only. Full contract: `../../AGENTS.md` §"Design system"; visual reference `docs/design-system.html`; system `DESIGN.md` (root).
