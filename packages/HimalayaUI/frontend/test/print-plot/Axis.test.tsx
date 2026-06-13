import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { makeAxis, axisTicks } from "../../src/print/plot/projection";
import { Axis } from "../../src/print/plot/Axis";

describe("Axis", () => {
  it("renders one tick label per labelled tick value (bottom)", () => {
    const a = makeAxis([0.01, 1], [0, 300], "log");
    const labelled = axisTicks(a).filter((t) => t.kind !== "minor");
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
    expect(texts.length).toBe(labelled.length);
  });

  it("paints tick labels in ink-soft, the AA-normal data tone (CC-AXIS-TICK)", () => {
    // The q-values ARE the data; at 10.5px they are small text and must clear
    // AA-normal (4.5:1). ink-faint (3.16:1) is AA-large/non-text only, so tick
    // numbers must use ink-soft.
    const a = makeAxis([0.01, 1], [0, 300], "log");
    const { container } = render(
      <svg>
        <Axis axis={a} orientation="bottom" plotWidth={300} plotHeight={150} />
      </svg>,
    );
    const tick = container.querySelector('[data-role="axis-bottom"] text') as SVGTextElement;
    expect(tick).toBeTruthy();
    expect(tick.style.fill).toBe("var(--color-ink-soft)");
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

  it("renders a spine and decade gridlines on the bottom axis", () => {
    const ax = makeAxis([0.02, 0.23], [0, 400], "log");
    const { container } = render(
      <svg><Axis axis={ax} orientation="bottom" plotWidth={400} plotHeight={200} label="q (Å⁻¹)" /></svg>,
    );
    expect(container.querySelector('[data-role="axis-spine"]')).toBeTruthy();
    expect(container.querySelectorAll('[data-role="gridline"]').length).toBeGreaterThan(0);
    const texts = container.querySelectorAll('[data-role="axis-bottom"] text');
    expect(texts.length).toBeGreaterThan(0);
  });
});
