import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { MemberList, type MemberDatum } from "../../src/print/components/MemberList";

const members: MemberDatum[] = [
  { key: "a", phases: ["Pn3m"], variableValue: "1:0", dataLine: "a = 205 Å · q₁ 0.061 Å⁻¹" },
  { key: "b", phases: ["Pn3m", "Lamellar"], variableValue: "1:0.5", dataLine: "a 195 · d 60 Å" },
  { key: "c", phases: [], variableValue: "1:1.5", dataLine: "no Bragg peaks · q₁ —" },
];

describe("MemberList", () => {
  it("renders one row per member in order", () => {
    const { getAllByTestId } = render(<MemberList members={members} />);
    expect(getAllByTestId("series-member-row")).toHaveLength(3);
  });
  it("marks the hovered key hot and reports hover/leave", () => {
    const onHoverMember = vi.fn();
    const { getAllByTestId } = render(
      <MemberList members={members} hoveredKey="b" onHoverMember={onHoverMember} />,
    );
    const rows = getAllByTestId("series-member-row");
    expect(rows[1]!.getAttribute("data-hot")).toBe("true");
    expect(rows[0]!.getAttribute("data-hot")).toBeNull();
    fireEvent.mouseEnter(rows[0]!);
    expect(onHoverMember).toHaveBeenCalledWith("a");
    fireEvent.mouseLeave(rows[0]!);
    expect(onHoverMember).toHaveBeenCalledWith(undefined);
  });
});
