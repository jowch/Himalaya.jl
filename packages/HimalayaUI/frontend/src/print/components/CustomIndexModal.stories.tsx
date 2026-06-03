import { useState } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { CustomIndexModal } from "./CustomIndexModal";
import { IM3M } from "../comb/comb.fixtures";

const SYMS: Record<string, { param: string; def: number; min: number; max: number }> = {
  Pn3m:       { param: "a", def: 197, min: 120, max: 320 },
  Im3m:       { param: "a", def: 252, min: 120, max: 360 },
  Ia3d:       { param: "a", def: 218, min: 120, max: 360 },
  Lamellar:   { param: "d", def: 60,  min: 30,  max: 130 },
  Hexagonal:  { param: "a", def: 70,  min: 40,  max: 160 },
};

const SYMMETRIES = ["Pn3m", "Im3m", "Ia3d", "Lamellar", "Hexagonal"] as const;

const noop = () => {};

const meta = {
  title: "components/CustomIndexModal",
  component: CustomIndexModal,
  // Required props — overridden by the interactive render below.
  args: {
    open: true,
    onClose: noop,
    onCancel: noop,
    onAdd: noop,
    symmetries: SYMMETRIES,
    symmetry: "Im3m",
    paramName: "a",
    paramValue: "252",
    paramMin: 120,
    paramMax: 360,
    paramStep: 1,
    onSymmetryChange: noop,
    onParamChange: noop,
    previewSeries: IM3M,
    observed: [0.028, 0.04, 0.049, 0.057],
    fit: { landed: 3, total: 4, snapped: false },
  },
} satisfies Meta<typeof CustomIndexModal>;

export default meta;
type Story = StoryObj<typeof meta>;

function OpenDemo() {
  const [open, setOpen] = useState(true);
  const [symmetry, setSymmetry] = useState<string>("Im3m");
  const [paramValue, setParamValue] = useState<string>(String(SYMS["Im3m"]!.def));

  const cfg = SYMS[symmetry]!;

  const handleSymmetryChange = (s: string) => {
    setSymmetry(s);
    setParamValue(String(SYMS[s]!.def));
  };

  return (
    <div>
      <button onClick={() => setOpen(true)}>Open custom index</button>
      <CustomIndexModal
        open={open}
        onClose={() => setOpen(false)}
        onCancel={() => setOpen(false)}
        onAdd={() => { setOpen(false); }}
        symmetries={SYMMETRIES}
        symmetry={symmetry}
        onSymmetryChange={handleSymmetryChange}
        paramName={cfg.param}
        paramValue={paramValue}
        paramMin={cfg.min}
        paramMax={cfg.max}
        paramStep={1}
        onParamChange={setParamValue}
        unit="Å"
        previewSeries={IM3M}
        observed={[0.028, 0.04, 0.049, 0.057]}
        fit={{ landed: 3, total: 4, snapped: false }}
      />
    </div>
  );
}

export const Open: Story = {
  render: () => <OpenDemo />,
};
