import http from "node:http";

export async function listen(handler) {
  const server = http.createServer(handler);
  await new Promise((r) => server.listen(0, "127.0.0.1", r));
  return { server, port: server.address().port, close: () => new Promise((r) => { server.closeAllConnections(); server.close(r); }) };
}

export async function stubBackend(routes) {
  const seen = [];
  const sent = [];
  const api = await listen(async (req, res) => {
    const key = req.method === "GET" ? req.url : `${req.method} ${req.url}`;
    seen.push(key);
    if (req.method !== "GET") {
      let raw = "";
      for await (const chunk of req) raw += chunk;
      sent.push({ key, body: raw ? JSON.parse(raw) : null });
    }
    const body = routes[key];
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
    sent,
    env: { DEMO_API_PORT: String(api.port), DEMO_UI_PORT: String(ui.port) },
    close: async () => {
      await api.close();
      await ui.close();
    },
  };
}
