import { render, screen, fireEvent } from "@testing-library/react";
import { KeptCell } from "../src/print/components/KeptCell";

test("KeptCell shows a Restore button only when something is dropped", () => {
  const onRestore = vi.fn();
  const { rerender } = render(
    <KeptCell kept={3} dropped={0} total={3} onRestore={onRestore} />,
  );
  expect(screen.queryByRole("button", { name: /restore/i })).toBeNull();
  rerender(<KeptCell kept={2} dropped={1} total={3} onRestore={onRestore} />);
  fireEvent.click(screen.getByRole("button", { name: /restore/i }));
  expect(onRestore).toHaveBeenCalledOnce();
});
