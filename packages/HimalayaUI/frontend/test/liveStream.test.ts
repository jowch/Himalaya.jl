import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { openLiveStream } from "../src/lib/liveStream";

class FakeES {
  static all: FakeES[] = [];
  readyState = 0;
  closed = false;
  listeners: Record<string, ((e: MessageEvent) => void)[]> = {};
  constructor(public url: string) { FakeES.all.push(this); }
  addEventListener(t: string, fn: (e: MessageEvent) => void) { (this.listeners[t] ??= []).push(fn); }
  close() { this.closed = true; this.readyState = 2; }
  fire(type: string, data = "{}") {
    if (type === "open") this.readyState = 1;
    for (const fn of this.listeners[type] ?? []) fn({ data } as MessageEvent);
  }
  failForGood() { this.readyState = 2; this.fire("error"); }
}

const last = () => FakeES.all[FakeES.all.length - 1]!;

describe("openLiveStream", () => {
  let reconnects: number;
  let frames: string[];
  let stop: () => void;

  beforeEach(() => {
    vi.useFakeTimers();
    FakeES.all = [];
    reconnects = 0;
    frames = [];
    stop = openLiveStream("/api/events", {
      listeners: { curation: (e) => frames.push(e.data as string) },
      onReconnect: () => { reconnects += 1; },
      createSource: (url) => new FakeES(url) as unknown as EventSource,
    });
  });
  afterEach(() => { stop(); vi.useRealTimers(); });

  it("passes frames through and does not count the first open as a reconnect", () => {
    last().fire("open");
    last().fire("curation", "x");
    expect(frames).toEqual(["x"]);
    expect(reconnects).toBe(0);
  });

  it("counts the browser's own retry as a reconnect", () => {
    last().fire("open");
    last().fire("error");
    last().fire("open");
    expect(FakeES.all).toHaveLength(1);
    expect(reconnects).toBe(1);
  });

  it("opens a new stream with backoff when the browser gives up after a non-200 answer", () => {
    last().fire("open");
    last().failForGood();
    expect(FakeES.all).toHaveLength(1);
    vi.advanceTimersByTime(1000);
    expect(FakeES.all).toHaveLength(2);
    last().failForGood();
    vi.advanceTimersByTime(1000);
    expect(FakeES.all).toHaveLength(2);
    vi.advanceTimersByTime(1000);
    expect(FakeES.all).toHaveLength(3);
    last().fire("open");
    expect(reconnects).toBe(1);
  });

  it("caps the wait between attempts at 5 s", () => {
    last().fire("open");
    for (const wait of [1000, 2000, 4000, 5000, 5000]) {
      const before = FakeES.all.length;
      last().failForGood();
      vi.advanceTimersByTime(wait - 1);
      expect(FakeES.all).toHaveLength(before);
      vi.advanceTimersByTime(1);
      expect(FakeES.all).toHaveLength(before + 1);
    }
  });

  it("starts the wait from 1 s again after a stream opens", () => {
    last().fire("open");
    last().failForGood();
    vi.advanceTimersByTime(1000);
    last().failForGood();
    vi.advanceTimersByTime(2000);
    expect(FakeES.all).toHaveLength(3);
    last().fire("open");
    last().failForGood();
    vi.advanceTimersByTime(999);
    expect(FakeES.all).toHaveLength(3);
    vi.advanceTimersByTime(1);
    expect(FakeES.all).toHaveLength(4);
  });

  it("counts live updates as signs of life, so a busy stream is never replaced", () => {
    last().fire("open");
    for (let i = 0; i < 6; i++) {
      vi.advanceTimersByTime(30_000);
      last().fire("curation", `f${i}`);
    }
    expect(FakeES.all).toHaveLength(1);
    expect(frames).toHaveLength(6);
  });

  it("replaces a stream that has gone silent, since a proxy can hold it open after the backend dies", () => {
    last().fire("open");
    vi.advanceTimersByTime(20_000);
    last().fire("heartbeat");
    vi.advanceTimersByTime(30_000);
    expect(FakeES.all).toHaveLength(1);
    vi.advanceTimersByTime(15_000);
    expect(FakeES.all[0]!.closed).toBe(true);
    expect(FakeES.all).toHaveLength(2);
    last().fire("open");
    expect(reconnects).toBe(1);
  });

  it("stops everything when disposed", () => {
    last().fire("open");
    stop();
    expect(last().closed).toBe(true);
    vi.advanceTimersByTime(120_000);
    expect(FakeES.all).toHaveLength(1);
  });
});
