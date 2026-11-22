const CLOSED = 2;

export class FakeEventSource {
  static instances: FakeEventSource[] = [];

  static install(): void {
    FakeEventSource.instances = [];
    (globalThis as unknown as { EventSource: unknown }).EventSource = FakeEventSource;
  }

  static current(): FakeEventSource {
    const live = FakeEventSource.instances.filter((s) => s.readyState !== CLOSED);
    if (live.length !== 1) {
      throw new Error(`expected one open EventSource, found ${live.length}`);
    }
    return live[0]!;
  }

  readyState = 0;
  private listeners: Record<string, ((e: Event) => void)[]> = {};

  constructor(public url: string) {
    FakeEventSource.instances.push(this);
  }

  addEventListener(type: string, fn: (e: Event) => void): void {
    (this.listeners[type] ??= []).push(fn);
  }

  removeEventListener(): void {}

  close(): void {
    this.readyState = CLOSED;
    this.listeners = {};
  }

  fire(type: string): void {
    for (const fn of this.listeners[type] ?? []) fn(new Event(type));
  }

  emit(type: string, data: unknown): void {
    const event = new MessageEvent(type, { data: JSON.stringify(data) });
    for (const fn of this.listeners[type] ?? []) fn(event);
  }
}
