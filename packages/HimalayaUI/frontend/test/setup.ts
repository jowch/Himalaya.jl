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

// JSDOM does not implement window.matchMedia. Stub it so boneyard-js dark-mode
// detection doesn't throw during unit tests.
if (typeof window.matchMedia === "undefined") {
  Object.defineProperty(window, "matchMedia", {
    writable: true,
    value: (query: string): MediaQueryList => ({
      matches: false,
      media: query,
      onchange: null,
      addListener: () => {},
      removeListener: () => {},
      addEventListener: () => {},
      removeEventListener: () => {},
      dispatchEvent: () => false,
    }),
  });
}
