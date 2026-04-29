import "@testing-library/jest-dom/vitest";

// JSDOM does not implement ResizeObserver. Provide a no-op stub so components
// that observe layout changes (e.g. DetectorImage auto-rotate) don't throw.
class ResizeObserverStub {
  observe(): void {}
  unobserve(): void {}
  disconnect(): void {}
}
if (typeof globalThis.ResizeObserver === "undefined") {
  // @ts-expect-error JSDOM lacks ResizeObserver
  globalThis.ResizeObserver = ResizeObserverStub;
}
