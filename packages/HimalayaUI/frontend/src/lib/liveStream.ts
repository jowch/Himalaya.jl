const CLOSED = 2;

export interface LiveStreamOptions {
  listeners: Record<string, (e: MessageEvent) => void>;
  onReconnect: () => void;
  createSource?: (url: string) => EventSource;
  silentAfterMs?: number;
  checkEveryMs?: number;
  maxBackoffMs?: number;
}

export function openLiveStream(url: string, opts: LiveStreamOptions): () => void {
  const create = opts.createSource ?? ((u: string) => new EventSource(u));
  const silentAfterMs = opts.silentAfterMs ?? 40_000;
  const maxBackoffMs = opts.maxBackoffMs ?? 5_000;
  let source: EventSource | null = null;
  let openedBefore = false;
  let attempt = 0;
  let lastHeard = Date.now();
  let retryTimer: ReturnType<typeof setTimeout> | undefined;
  let disposed = false;

  const heard = () => { lastHeard = Date.now(); };

  const replace = (delayMs: number) => {
    source?.close();
    source = null;
    clearTimeout(retryTimer);
    if (delayMs === 0) connect();
    else retryTimer = setTimeout(connect, delayMs);
  };

  function connect(): void {
    if (disposed) return;
    const es = create(url);
    source = es;
    heard();
    es.addEventListener("open", () => {
      heard();
      attempt = 0;
      if (openedBefore) opts.onReconnect();
      openedBefore = true;
    });
    es.addEventListener("heartbeat", heard);
    for (const [type, fn] of Object.entries(opts.listeners)) {
      es.addEventListener(type, (e) => { heard(); fn(e as MessageEvent); });
    }
    es.addEventListener("error", () => {
      if (source !== es || es.readyState !== CLOSED) return;
      replace(Math.min(1000 * 2 ** attempt, maxBackoffMs));
      attempt += 1;
    });
  }

  const watchdog = setInterval(() => {
    if (source && source.readyState !== CLOSED && Date.now() - lastHeard > silentAfterMs) replace(0);
  }, opts.checkEveryMs ?? 5_000);

  connect();
  return () => {
    disposed = true;
    clearInterval(watchdog);
    clearTimeout(retryTimer);
    source?.close();
  };
}
