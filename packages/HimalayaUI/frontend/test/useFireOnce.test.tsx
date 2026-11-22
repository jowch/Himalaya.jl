import { describe, it, expect, vi } from "vitest";
import { renderHook } from "@testing-library/react";
import { useFireOnce } from "../src/hooks/useFireOnce";

type Props = { token: { path: string } | undefined; ready: boolean };

describe("useFireOnce", () => {
  it("fires a token once, even when a render that still holds it flips ready again", () => {
    const fire = vi.fn();
    const token = { path: "/b/day1/setup1/" };
    const { rerender } = renderHook(({ token, ready }: Props) => useFireOnce(token, ready, fire), {
      initialProps: { token, ready: false },
    });
    expect(fire).not.toHaveBeenCalled();
    rerender({ token, ready: true });
    rerender({ token, ready: false });
    rerender({ token, ready: true });
    rerender({ token, ready: true });
    expect(fire).toHaveBeenCalledTimes(1);
    expect(fire).toHaveBeenCalledWith(token);
  });

  it("fires again for a new token and never for none", () => {
    const fire = vi.fn();
    const { rerender } = renderHook(({ token, ready }: Props) => useFireOnce(token, ready, fire), {
      initialProps: { token: undefined, ready: true } as Props,
    });
    expect(fire).not.toHaveBeenCalled();
    rerender({ token: { path: "/a" }, ready: true });
    rerender({ token: { path: "/a" }, ready: true });
    expect(fire).toHaveBeenCalledTimes(2);
  });
});
