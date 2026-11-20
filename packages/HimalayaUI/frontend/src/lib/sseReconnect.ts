export function onReconnect(source: Pick<EventSource, "addEventListener">, run: () => void): void {
  let opened = false;
  source.addEventListener("open", () => {
    if (opened) run();
    opened = true;
  });
}
