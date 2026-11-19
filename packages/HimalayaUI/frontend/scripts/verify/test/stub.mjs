import http from "node:http";

export async function listen(handler) {
  const server = http.createServer(handler);
  await new Promise((r) => server.listen(0, "127.0.0.1", r));
  return { server, port: server.address().port, close: () => new Promise((r) => { server.closeAllConnections(); server.close(r); }) };
}

export async function stubBackend(routes) {
  const seen = [];
  const api = await listen((req, res) => {
    seen.push(req.url);
    const body = routes[req.url];
    if (body === undefined) {
      res.statusCode = 404;
      res.end(JSON.stringify({ error: `no stub for ${req.url}` }));
      return;
    }
    res.setHeader("Content-Type", "application/json");
    res.end(JSON.stringify(typeof body === "function" ? body() : body));
  });
  const ui = await listen((req, res) => {
    if (req.url === "/api/health") {
      res.setHeader("Content-Type", "application/json");
      res.end('{"status":"ok"}');
      return;
    }
    res.setHeader("Content-Type", "text/html");
    res.end('<!doctype html><div id="app"></div>');
  });
  return {
    seen,
    env: { DEMO_API_PORT: String(api.port), DEMO_UI_PORT: String(ui.port) },
    close: async () => {
      await api.close();
      await ui.close();
    },
  };
}
