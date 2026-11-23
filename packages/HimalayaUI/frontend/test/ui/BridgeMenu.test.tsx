import { describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen, within } from "@testing-library/react";
import { BridgeMenu } from "../../src/ui/BridgeMenu";

const mates = [{ id: 12, title: "A3 in Long pipe" }, { id: 13, title: "A4 in Long pipe" }];

describe("BridgeMenu", () => {
  it("shows a bridge icon with the count, named after every same sample", () => {
    render(<BridgeMenu mates={mates} editable onGo={vi.fn()} onRemove={vi.fn()} />);
    const button = screen.getByRole("button", { name: "Same as A3 in Long pipe, A4 in Long pipe" });
    expect(button).toHaveTextContent("⇄2");
    expect(button).toHaveAttribute("aria-expanded", "false");
  });

  it("a click lists the same samples; picking one goes there and closes the list", () => {
    const onGo = vi.fn();
    render(<BridgeMenu mates={mates} editable onGo={onGo} onRemove={vi.fn()} />);
    fireEvent.click(screen.getByRole("button", { name: /^Same as / }));
    const list = screen.getByRole("dialog", { name: "Same sample as" });
    expect(within(list).getByRole("button", { name: "Go to A3 in Long pipe" })).toHaveFocus();
    fireEvent.click(within(list).getByRole("button", { name: "Go to A4 in Long pipe" }));
    expect(onGo).toHaveBeenCalledWith(13);
    expect(screen.queryByRole("dialog")).toBeNull();
  });

  it("takes one sample out of the bridge, and offers that only when editable", () => {
    const onRemove = vi.fn();
    const { rerender } = render(<BridgeMenu mates={mates} editable onGo={vi.fn()} onRemove={onRemove} />);
    fireEvent.click(screen.getByRole("button", { name: /^Same as / }));
    fireEvent.click(screen.getByRole("button", { name: "Not the same sample as A3 in Long pipe" }));
    expect(onRemove).toHaveBeenCalledWith(12);
    rerender(<BridgeMenu mates={mates} editable={false} onGo={vi.fn()} onRemove={onRemove} />);
    expect(screen.queryByRole("button", { name: /^Not the same sample/ })).toBeNull();
  });

  it("Esc closes the list and returns focus to the icon", () => {
    render(<BridgeMenu mates={mates} editable onGo={vi.fn()} onRemove={vi.fn()} />);
    const button = screen.getByRole("button", { name: /^Same as / });
    fireEvent.click(button);
    fireEvent.keyDown(screen.getByRole("dialog", { name: "Same sample as" }), { key: "Escape" });
    expect(screen.queryByRole("dialog")).toBeNull();
    expect(button).toHaveFocus();
  });

  it("after × keeps focus in the list on the next same sample, so Esc still closes it", () => {
    const { rerender } = render(<BridgeMenu mates={mates} editable onGo={vi.fn()} onRemove={vi.fn()} />);
    fireEvent.click(screen.getByRole("button", { name: /^Same as / }));
    const remove = screen.getByRole("button", { name: "Not the same sample as A3 in Long pipe" });
    remove.focus();
    fireEvent.click(remove);
    rerender(<BridgeMenu mates={[mates[1]!]} editable onGo={vi.fn()} onRemove={vi.fn()} />);
    expect(screen.getByRole("button", { name: "Go to A4 in Long pipe" })).toHaveFocus();
    fireEvent.keyDown(document.activeElement!, { key: "Escape" });
    expect(screen.queryByRole("dialog")).toBeNull();
  });

  it("closes when focus leaves it, as when Tab moves on", () => {
    render(
      <>
        <BridgeMenu mates={mates} editable onGo={vi.fn()} onRemove={vi.fn()} />
        <button type="button">elsewhere</button>
      </>,
    );
    fireEvent.click(screen.getByRole("button", { name: /^Same as / }));
    const inside = screen.getByRole("button", { name: "Go to A3 in Long pipe" });
    fireEvent.blur(inside, { relatedTarget: screen.getByRole("button", { name: "Not the same sample as A3 in Long pipe" }) });
    expect(screen.getByRole("dialog", { name: "Same sample as" })).toBeInTheDocument();
    fireEvent.blur(inside, { relatedTarget: screen.getByRole("button", { name: "elsewhere" }) });
    expect(screen.queryByRole("dialog")).toBeNull();
  });
});
