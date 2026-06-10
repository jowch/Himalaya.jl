# frontend/src/print/shell — app shell, routing, modals

The live application shell: the single layout shell, the topbar chrome, the route table, and the global onboarding/nav/error overlays. Everything else in the frontend is composed under `src/print/` (`ui/` primitives, `components/` composites, `pages/`, the `plot`/`detector`/`comb`/`waterfall` render layers, `export/` figure utilities) plus the non-presentational layers at `src/` root (`App.tsx`, `state.ts`, `queries.ts`, `api.ts`, `lib/`, `hooks/`).

## Where things live

- **Shell + chrome**: `CorpusShell` (the single layout shell) wraps `CorpusTopbar` (wordmark + stage-tabs + Beamtime chip — composed from `../ui` primitives: `TopBar`, `Wordmark`, `Kicker`, `IconButton`). No appearance utilities live here — chrome appearance is owned by the `ui/` primitives (closed-look / open-placement contract).
- **Routing**: `AppRoutes` is the route table (samples / loupe / focus / series folio·builder·scoping pages, all imported from `../pages`). `IndexSlugRedirect` and `StaleUrlPage` resolve/repair permalink slugs; `ResolvingFallback` is the in-flight placeholder.
- **Onboarding / modals / errors**: `OnboardingFlow`, `NavModal`, `InfrastructureBanner`. The modals compose `ModalShell` (+ `Button`/`Kicker`/`IconButton`) from `../ui`.

## Import boundary

These files live under `src/print/`, so the `lint:design` import-boundary guard (`scripts/check-design.mjs`, `scanLegacyImports`) forbids them from importing a top-level `src/components/**` or `src/pages/**` — both of which are now retired. Reach for state/data via `../../state`, `../../queries`, `../../api`, `../../hooks`, `../../lib`; reach for presentation via the sibling `../ui` primitives and `../pages` route components.

## Focus trapping in modals

`../../hooks/useFocusTrap.ts` exports `useFocusTrap(containerRef, active)` (APG modal-dialog pattern). On activation it moves focus INTO the container (first focusable, or the container via tabIndex=-1) unless the consumer already placed it inside — `autoFocus` runs before the trap effect and is respected; a parent's own focus effect runs after it (child-before-parent order) and wins, which is how `NavModal`'s input takes focus. It intercepts Tab/Shift+Tab at the **document** level (background controls stay unreachable even if focus escapes), skips Tabs a consumer already `preventDefault`ed (NavModal's Tab-commit), pulls escaped focus back via a `focusin` guard, and restores the previously-focused element on cleanup. `ModalShell` wires it for every modal; `NavModal` and `OnboardingFlow` get it through that.

## Anti-patterns

- Don't introspect tooltip text or Tailwind classes in tests — use `data-*` attributes. Tests stay behavioral.
- Don't add inline appearance utilities (`text-[…]`, `rounded-[…]`, raw colors, side-stripes) here — they fail the `lint:design` guard. Appearance lives in `../ui` primitives; a consumer's `className` is placement-only. Full contract: `../../AGENTS.md` §"Design system"; visual reference `docs/design-system.html`; system `DESIGN.md` (root).
