import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { makeAxis } from "../../src/print/plot/projection";
import { Axis } from "../../src/print/plot/Axis";

describe("Axis", () => {
  it("renders one tick label per tick value (bottom)", () => {
    const a = makeAxis([0.01, 1], [0, 300], "log");
    const ticks = a.ticks(6);
    const { container } = render(
      <svg>
        <Axis
          axis={a}
          orientation="bottom"
          plotWidth={300}
          plotHeight={150}
        />
      </svg>,
    );
    const texts = container.querySelectorAll('[data-role="axis-bottom"] text');
    expect(texts.length).toBe(ticks.length);
  });

  it("renders a label when provided (left)", () => {
    const a = makeAxis([1, 1000], [150, 0], "log");
    const { container } = render(
      <svg>
        <Axis
          axis={a}
          orientation="left"
          plotWidth={300}
          plotHeight={150}
          label="I"
        />
      </svg>,
    );
    const labels = [...container.querySelectorAll("text")].filter(
      (t) => t.textContent === "I",
    );
    expect(labels.length).toBe(1);
  });
});
