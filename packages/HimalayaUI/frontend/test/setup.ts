import "@testing-library/jest-dom/vitest";

// JSDOM does not implement ResizeObserver. Provide a no-op stub so components
// that observe layout changes (e.g. DetectorImage auto-rotate) don't throw.
class ResizeObserverStub {
  observe(): void {}
  unobserve(): void {}
  disconnect(): void {}
}
if (typeof globalThis.ResizeObserver === "undefined") {
  globalThis.ResizeObserver = ResizeObserverStub;
}

// JSDOM does not implement EventSource. Provide a no-op stub so App.tsx's SSE
// subscriber effect doesn't throw when components that render App are tested.
class EventSourceStub {
  static CONNECTING = 0;
  static OPEN = 1;
  static CLOSED = 2;
  CONNECTING: 0 = 0;
  OPEN: 1 = 1;
  CLOSED: 2 = 2;
  readyState: number = 1;
  url: string;
  withCredentials: boolean = false;
  onopen: ((this: EventSource, ev: Event) => unknown) | null = null;
  onmessage: ((this: EventSource, ev: MessageEvent) => unknown) | null = null;
  onerror: ((this: EventSource, ev: Event) => unknown) | null = null;
  constructor(url: string) { this.url = url; }
  addEventListener(): void {}
  removeEventListener(): void {}
  dispatchEvent(): boolean { return true; }
  close(): void { this.readyState = 2; }
}
if (typeof globalThis.EventSource === "undefined") {
  // @ts-expect-error JSDOM lacks EventSource
  globalThis.EventSource = EventSourceStub;
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
