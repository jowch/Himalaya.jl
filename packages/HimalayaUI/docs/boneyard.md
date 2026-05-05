# Skeleton loading via boneyard-js

HimalayaUI uses [`boneyard-js`](https://www.npmjs.com/package/boneyard-js)
to render layout-aware skeletons during cold data fetches. Each load-gated
card or list wraps its content in `<Skeleton>` from `boneyard-js/react`.

This doc covers the non-obvious rules. All of them are load-bearing — each
exists to fix a specific failure mode that a "naive" use of the library
hits in this codebase.

---

## Where it's used

| Component | Skeleton name |
|---|---|
| `PlotCard` | trace plot |
| `PhasePanel` | indices list |
| `ChatCard` | message history |
| `DetectorImageCard` | detector image |
| `SampleMetadataCard` | metadata fields |
| `NavModal` | nav modal contents |

`<Skeleton>` is the only consumer-facing primitive. Each call site supplies
a `name` (matches a captured `*.bones.json`), a `loading` boolean, a
`fixture` (JSX rendered during headless capture), a `fallback` (rendered
when no bones exist for `name`), and a layout `className`.

---

## Rule 1: gate on `query.isLoading`, not `query.isPending`

```tsx
<Skeleton name="trace" loading={traceQuery.isLoading} fixture={…} fallback={…}>
  <TraceViewer trace={traceQuery.data} … />
</Skeleton>
```

TanStack v5 distinguishes `isLoading` and `isPending`:

```
isLoading  = isPending && isFetching   // true cold fetch
isPending  = data === undefined        // also true for disabled queries
```

A query becomes "disabled" when its selection key is `undefined` (e.g.
`useExposure(undefined)` while no exposure is selected). Disabled queries
have `data === undefined` and `isPending === true` forever, but they
shouldn't render a skeleton — there's no fetch happening. Background
refetches over already-cached data also have `isPending === false` (data
is defined) but `isFetching === true`. We don't want skeleton flicker on
every refetch either.

`isLoading` is the conjunction that captures *only* true cold fetches —
the only state where a skeleton is the right UX.

Wrong gating ⇒ skeleton flicker on every background refetch, or persistent
skeleton on disabled queries.

---

## Rule 2: `fixture` is `ReactNode`, not raw data

```tsx
<Skeleton
  name="trace"
  fixture={<TraceViewer trace={MOCK_TRACE} peaks={MOCK_PEAKS} onSelect={() => {}} />}
  …
/>
```

The `fixture` prop is JSX rendered by boneyard's headless capture CLI to
*measure layout*. The CLI walks the rendered DOM, captures bounding-box
geometry per element, and emits `*.bones.json`. The skeleton drawn at
runtime mirrors that geometry.

Pass real components with mock props and no-op handlers. Don't pass raw
data shapes (mock objects, arrays) — there's no rendering pipeline at
capture time that would turn those into a DOM tree.

---

## Rule 3: always set `fallback`

```tsx
<Skeleton name="trace" fallback={<HintText italic>Loading trace…</HintText>} …>
```

Without bones for `name` (e.g. on first dev-server boot, or in prod before
the first capture has been committed), the runtime renders `fallback`
during loading. With no `fallback`, the area is blank.

Mirror the original component's loading-state hint text. Italic
`HintText` is the project convention.

---

## Rule 4: `className` on `<Skeleton>` is load-bearing

```tsx
<Skeleton name="chat" className="flex-1 min-h-0 flex flex-col" …>
  <MessageList … />
</Skeleton>
```

Boneyard wraps children in two extra `<div>`s. Without an explicit
`className`, the outer wrapper has no flex semantics — and any parent that
relied on the child being a flex item silently breaks. ChatCard's message
list collapsing to ~60px in the original PR was this exact failure.

Pass the layout role the original child *would have had*:
`flex-1 min-h-0 flex flex-col` for vertically-stretching content,
`h-full w-full` for free-fill.

Companion CSS in `styles.css` makes the *inner* wrapper transparent to
layout:

```css
[data-boneyard-content] { display: contents; }
```

`display: contents` removes the element from the box tree without removing
its children, so the children's flex behaviour passes through to the
outer wrapper unchanged.

---

## Rule 5: `configureBoneyard()` lives in `main.tsx`

Not in `bones/registry.ts`. The Vite HMR plugin regenerates `registry.ts`
on every capture and would wipe any config call placed there.

The values must mirror `boneyard.config.json` (which the capture CLI
reads but the runtime does not). When the card background colour or
animation parameters change, **update both files together**:

- `boneyard.config.json` — read by the capture CLI to colour the bones
- `main.tsx::configureBoneyard()` — read by the runtime to animate them

If they drift, captured bones won't match the runtime palette and the
skeleton looks wrong against the live theme.

---

## Rule 6: bones are committed

`src/bones/*.bones.json` and the auto-generated `src/bones/registry.ts`
are checked into git. They're required for prod builds to render
skeletons; without them, `<Skeleton>` falls through to `fallback` and the
user sees the italic hint text instead of a skeleton.

Capture *organically* during dev — the Vite HMR plugin re-captures on
every update. Commit deliberately to widen prod coverage. The bones do
not change with code unless layout or fixture content changes; review the
diff before committing to make sure no spurious geometry shifts crept in.

---

## Rule 7: JSDOM lacks `window.matchMedia`

Boneyard's dark-mode detection calls `window.matchMedia` on mount. JSDOM
doesn't ship it. The stub in `frontend/test/setup.ts` keeps unit tests
honest:

```ts
// test/setup.ts
Object.defineProperty(window, "matchMedia", {
  value: () => ({ matches: false, addEventListener: () => {}, … }),
});
```

Same pattern as the `ResizeObserver` stub used by `DetectorImage`. If you
add a new test environment (Vitest workspace, Storybook), copy the stubs
or component-tests using `<Skeleton>` will throw on import.

---

## Adding a new skeleton

1. Wrap the load-gated content in `<Skeleton>` with a unique `name`,
   `loading={query.isLoading}`, a real-component `fixture`, an italic
   `fallback`, and the right `className`.
2. Run `npm run dev` and trigger the loading state — the Vite plugin
   captures bones automatically and writes `src/bones/<name>.bones.json`.
3. Verify the captured geometry matches the original layout.
4. Commit the new `*.bones.json` together with the component change.

If the skeleton is for a card that already has a fixed-width layout (e.g.
detector image's 1:1 aspect-ratio frame), still set `className` so the
outer wrapper inherits the framing.

---

## Files

- `frontend/src/main.tsx` — `configureBoneyard()` call
- `frontend/src/bones/` — committed `*.bones.json` + auto-generated `registry.ts`
- `frontend/boneyard.config.json` — capture-CLI config (mirror of `configureBoneyard()`)
- `frontend/src/styles.css` — `[data-boneyard-content] { display: contents }`
- `frontend/test/setup.ts` — `matchMedia` and `ResizeObserver` stubs
