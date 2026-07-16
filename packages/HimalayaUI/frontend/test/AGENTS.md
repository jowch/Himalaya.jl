# frontend/test — Vitest + React Testing Library

Unit tests for `frontend/src/`. JSDOM environment.

## Capture once when slicing a slow run

```bash
npm test > /tmp/vitest.out 2>&1
grep -E "FAIL|passed|failed" /tmp/vitest.out
```

Single-file run:

```bash
node_modules/.bin/vitest run test/print-detector/DetectorImage.test.tsx
```

## RTL over `document.querySelector`

For DOM assertions use `screen.getByText("X").closest("li")` + `toHaveAttribute` rather than `document.querySelector` — the latter bypasses RTL's async-aware retry logic and produces flaky tests under React 18 concurrent rendering.

## `data-testid` / `data-*` over class strings

Never assert on Tailwind class strings — they change when styling evolves. Use stable `data-*` attributes that the component intentionally exposes (`data-active`, `data-orient`, etc.).

## JSDOM `fetch` interceptor pattern

`new Request(input, init)` **throws under JSDOM** (no constructor in the test environment). When spying on global `fetch` to assert headers/URLs:

```ts
const url = typeof input === "string" ? input : String(input);
const headers = init?.headers;
```

Don't construct a `Request` to introspect the call. See `test/queue/authHeaders.test.ts` for the canonical pattern.

## ResizeObserver stub

JSDOM doesn't ship `ResizeObserver`. `test/setup.ts` provides a no-op stub so components like `DetectorImage` (whose auto-rotate logic depends on it) can render. Don't strip the stub even if a passing test seems not to need it — the production component code path still references the global.

## JSDOM ignores `inert`

The hidden floating bars (CullBar, ComposeBar) hide via the `inert` attribute. JSDOM applies the attribute but implements none of its semantics: RTL role queries still match inert content, and `fireEvent` happily clicks inert buttons. Scope queries with `within(...)` when both bars are mounted, assert the attribute (`toHaveAttribute("inert")`) for the hide contract, and leave focus/AT semantics to the e2e real-browser pin (`e2e/corpus-culling.spec.ts`).

## `DetectorImage` object-URL stub

`DetectorImage` recolors via a GPU SVG filter on a plain `<img>`: it fetches the image blob and shows it through `URL.createObjectURL`. JSDOM implements neither object-URL method, so `test/setup.ts` stubs both globally (any component rendering a detector image would otherwise throw). Tests asserting the value re-stub `URL.createObjectURL` locally (`() => "blob:mock-url"`) and assert `revokeObjectURL` for the unmount leak-freedom contract.

## Anti-patterns

- Don't import production singletons (`useAppState`, query client) without resetting them — leakage between tests causes order-dependent failures.
- Don't synthesize SSE events by hand-rolling JSON — use the `applyRemoteToCache` test helpers in `src/lib/queue/testHelpers.ts`.
