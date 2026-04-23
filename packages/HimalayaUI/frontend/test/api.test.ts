import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import * as api from "../src/api";

function mockFetchJson(body: unknown, status = 200): void {
  global.fetch = vi.fn(async () => new Response(JSON.stringify(body), {
    status, headers: { "Content-Type": "application/json" },
  })) as unknown as typeof fetch;
}

describe("api client", () => {
  beforeEach(() => { api.setUsername(undefined); });
  afterEach(() => { vi.restoreAllMocks(); });

  it("listUsers returns [] for empty response", async () => {
    mockFetchJson([]);
    expect(await api.listUsers()).toEqual([]);
  });

  it("listUsers parses user rows", async () => {
    mockFetchJson([{ id: 1, username: "alice" }]);
    const users = await api.listUsers();
    expect(users.length).toBe(1);
    expect(users[0]!.username).toBe("alice");
  });

  it("createUser POSTs JSON body with Content-Type", async () => {
    const spy = vi.fn(async () => new Response(JSON.stringify({ id: 1, username: "bob" }), {
      status: 201, headers: { "Content-Type": "application/json" },
    }));
    global.fetch = spy as unknown as typeof fetch;
    const user = await api.createUser("bob");
    expect(user.username).toBe("bob");
    const init = spy.mock.calls[0]![1];
    expect(init?.method).toBe("POST");
    expect((init?.headers as Record<string, string>)["Content-Type"]).toBe("application/json");
    expect(init?.body).toBe(JSON.stringify({ username: "bob" }));
  });

  it("mutating calls send X-Username when set", async () => {
    const spy = vi.fn(async () => new Response("null", {
      status: 200, headers: { "Content-Type": "application/json" },
    }));
    global.fetch = spy as unknown as typeof fetch;
    api.setUsername("alice");
    await api.updateSample(1, { notes: "hi" });
    const init = spy.mock.calls[0]![1];
    expect((init?.headers as Record<string, string>)["X-Username"]).toBe("alice");
  });

  it("non-2xx throws ApiError with status and message", async () => {
    mockFetchJson({ error: "not found" }, 404);
    await expect(api.getExperiment(99)).rejects.toMatchObject({
      status: 404,
      message: expect.stringContaining("not found"),
    });
  });
});
