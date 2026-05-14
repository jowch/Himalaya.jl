# frontend/test — Vitest + React Testing Library

Unit tests for `frontend/src/`. JSDOM environment.

## Capture once when slicing a slow run

```bash
npm test > /tmp/vitest.out 2>&1
grep -E "FAIL|passed|failed" /tmp/vitest.out
```

Single-file run:

```bash
node_modules/.bin/vitest run test/DetectorImage.test.tsx
```

## RTL over `document.querySelector`

For DOM assertions use `screen.getByText("X").closest("li")` + `toHaveAttribute` rather than `document.querySelector` — the latter bypasses RTL's async-aware retry logic and produces flaky tests under React 18 concurrent rendering.

## `data-testid` / `data-*` over class strings

Never assert on Tailwind class strings — they change when styling evolves. Use stable `data-*` attributes that the component intentionally exposes (`data-index-id`, `data-alternative-id`, `data-active`, `data-low-r2`, `data-speculative`, etc.).

## JSDOM `fetch` interceptor pattern

`new Request(input, init)` **throws under JSDOM** (no constructor in the test environment). When spying on global `fetch` to assert headers/URLs:

```ts
const url = typeof input === "string" ? input : String(input);
const headers = init?.headers;
```

Don't construct a `Request` to introspect the call. See `test/queue/authHeaders.test.ts` for the canonical pattern.

## ResizeObserver stub

JSDOM doesn't ship `ResizeObserver`. `test/setup.ts` provides a no-op stub so components like `DetectorImage` (whose auto-rotate logic depends on it) can render. Don't strip the stub even if a passing test seems not to need it — the production component code path still references the global.

## `ImageBitmap.close()` regression

A getter-based mock in `test/DetectorImage.test.tsx` simulates the spec-mandated width/height-zeroing on `close()`. Keep it green — production code captures `{ width, height }` before closing.

## Anti-patterns

- Don't import production singletons (`useAppState`, query client) without resetting them — leakage between tests causes order-dependent failures.
- Don't synthesize SSE events by hand-rolling JSON — use the `applyRemoteToCache` test helpers in `src/lib/queue/testHelpers.ts`.
