# frontend/src/print/shell — app shell, routing, modals

The live application shell: the single layout shell, the topbar chrome, the route table, and the global onboarding/nav/error overlays. Everything else in the frontend is composed under `src/print/` (`ui/` primitives, `components/` composites, `pages/`, the `plot`/`detector` render layers), the triage workspace under `src/workspace/` — the app shell itself is `src/print/App.tsx` (`PrintApp`) — plus the non-presentational layers at `src/` root (`state.ts`, `queries.ts`, `api.ts`, `lib/`, `hooks/`).

## Where things live

- **Shell + chrome**: `AppShell` (the single layout shell) wraps `TopNav` (wordmark + stage-tabs + Experiment chip — composed from `../ui` primitives: `TopBar`, `Wordmark`). No appearance utilities live here — chrome appearance is owned by the `ui/` primitives (closed-look / open-placement contract).
- **Routing**: `AppRoutes` is the route table (experiments home and new experiment, imported from `../pages`). `/experiments/:id` renders `ExperimentEntry`, which forwards to the workspace once the experiment loads, whatever its ingest status; the workspace shows scan progress and failures itself. The triage workspace (`/experiments/:id/workspace/:sampleId?`, `src/workspace/`) mounts **outside** `AppShell` with its own header and its own keydown listener; `?exposure=<id>` opens it on that frame. Build links with `workspaceHref` from `lib/url/workspaceHref.ts`. Every other URL renders `LinkGone`, which toasts `PAGE_GONE` and redirects to `/experiments`; `ExperimentEntry` does the same for a deleted experiment.
- **Onboarding / modals / errors**: `OnboardingFlow`, `ShellPalette`, `InfrastructureBanner`. The modals compose `ModalShell` (+ `Button`/`Kicker`/`IconButton`) from `../ui`.

## Import boundary

These files live under `src/print/`, so the `lint:design` import-boundary guard (`scripts/check-design.mjs`, `scanLegacyImports`) forbids them from importing a top-level `src/components/**` or `src/pages/**` — both of which are now retired. Reach for state/data via `../../state`, `../../queries`, `../../api`, `../../hooks`, `../../lib`; reach for presentation via the sibling `../ui` primitives and `../pages` route components.

## Focus trapping in modals

`../../hooks/useFocusTrap.ts` exports `useFocusTrap(containerRef, active)` (APG modal-dialog pattern). On activation it moves focus INTO the container (first focusable, or the container via tabIndex=-1) unless the consumer already placed it inside — `autoFocus` runs before the trap effect and is respected; a parent's own focus effect runs after it (child-before-parent order) and wins. It intercepts Tab/Shift+Tab at the **document** level (background controls stay unreachable even if focus escapes), skips Tabs a consumer already `preventDefault`ed, pulls escaped focus back via a `focusin` guard, and restores the previously-focused element on cleanup. `ModalShell` wires it for every modal; `OnboardingFlow` gets it through that.

## Keyboard

Pages inside `AppShell` get three shell keys from `useShellKeys` (one window listener, mounted in `AppShell`): `⌘K`/`Ctrl+K` and `/` open the palette (`ShellPalette`: the workspace's `Palette` engine with `shellPalette(...)` from `shellPalette.ts`, open state in the store as `paletteOpen` so `StaleUrlPage` can open it), `?` opens the shortcuts list. Its App group (Theme, Switch user, Keyboard shortcuts) is `appItems`/`themePage` from `src/workspace/palette.ts`, shared with the workspace palette; add app-wide commands there so both palettes get them. The list reuses the workspace's `HelpOverlay`, with sections chosen per route by `shellHelp(pathname)`; add a section there when a page gains keys. Every shell key goes through `suppressGlobalKeys` (`lib/keys.ts`), so nothing fires while typing, inside a dialog, or after a closer handler `preventDefault`ed. The one exception is `?` in an empty text field (`helpKeyFromEmptyField`), so help is reachable from the autofocused funnel field. The experiments home moves focus between cards with the arrow keys (`gridStep` in `lib/gridFocus.ts`: ←/→ reading order, ↑/↓ the nearest card in the next row). The workspace has its own listener and keymap (see `src/workspace/`) but the same guards: `suppressGlobalKeys` for page keys, `helpKeyFromEmptyField` for `?`, and Esc in a field leaves the field.

The app-wide rule for Tab (approved 2026-11-21): in a list of choices Tab/Shift+Tab move the highlight with wrap while focus stays put, Enter commits, Esc backs out; forms (the new-experiment funnel, confirm dialogs, onboarding) keep native Tab. The funnel's folder field keeps Tab-to-complete.

## Key routes on `:id` when they hold ephemeral per-id local state

React reuses the same element instance across param changes on a single route, so page-level LOCAL state (undo/redo stacks, view offset/scale, Confirm-chain refs) bleeds from one id into the next unless the route element is keyed on the param. Any new `/foo/:id` route with ephemeral per-visit local state needs the same wrapper. (The global draft slot is intentionally recoverable on return; the key resets only the per-visit state that should never cross an id boundary.)

## OnboardingFlow defers `setUser` until `closeTutorial`

For a brand-new user who hasn't seen the tutorial, `onSubmitName` stashes the created user in local `pending` state and switches to the tutorial phase **instead of** calling `setUser()` immediately (`OnboardingFlow.tsx:118-124`). The reason: `setUser()` makes `username !== undefined`, which trips the `if (username !== undefined && phase !== "tutorial") return null` unmount guard (`:96`) and skips the tutorial. The pending user is committed only in `closeTutorial()` (`:130-135`, behind `setTutorialSeen(true)`). Calling `setUser()` at creation time silently breaks the tutorial for new users.

## Anti-patterns

- Don't introspect tooltip text or Tailwind classes in tests — use `data-*` attributes. Tests stay behavioral.
- Don't add inline appearance utilities (`text-[…]`, `rounded-[…]`, raw colors, side-stripes) here — they fail the `lint:design` guard. Appearance lives in `../ui` primitives; a consumer's `className` is placement-only. Full contract: `../../AGENTS.md` §"Design system"; visual reference Storybook (`npm run storybook`); system `DESIGN.md` (root).
