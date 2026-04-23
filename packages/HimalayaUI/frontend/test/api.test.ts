import { describe, it, expect, vi, beforeEach } from "vitest";
import * as api from "../src/api";
import { ApiError } from "../src/api";

describe("api", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("sends X-Username header on mutating calls when opts.username provided", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({ id: 1, username: "alice" }), { status: 200 }),
    );
    await api.createUser("alice", { username: "alice" });
    const init = fetchSpy.mock.calls[0]![1] as RequestInit;
    expect((init.headers as Record<string, string>)["X-Username"]).toBe("alice");
  });

  it("omits X-Username on GET requests", async () => {
    const fetchSpy = vi.spyOn(global, "fetch").mockResolvedValue(
      new Response("[]", { status: 200 }),
    );
    await api.listUsers();
    const init = fetchSpy.mock.calls[0]![1] as RequestInit;
    const headers = (init.headers as Record<string, string>) ?? {};
    expect(headers["X-Username"]).toBeUndefined();
  });

  it("does not export setUsername/getUsername", () => {
    expect((api as unknown as Record<string, unknown>).setUsername).toBeUndefined();
    expect((api as unknown as Record<string, unknown>).getUsername).toBeUndefined();
  });

  it("throws ApiError on non-2xx responses with server error message", async () => {
    vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({ error: "nope" }), { status: 400 }),
    );
    await expect(api.createUser("bob")).rejects.toMatchObject({
      name: "ApiError",
      status: 400,
      message: "nope",
    });
  });
});
