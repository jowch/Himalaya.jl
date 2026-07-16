import { phaseColor } from "../../phases";
import type { CombSeries } from "./combModel";

export const PN3M: CombSeries = {
  phase: "Pn3m", color: phaseColor("Pn3m"), latticeLabel: "a = 197 Å", rSquared: 0.998,
  teeth: [
    { q: 0.0712, label: "√2", observed: true, residual: 0.004 },
    { q: 0.0872, label: "√3", observed: true, residual: -0.006 },
    { q: 0.1007, label: "√4", observed: true, residual: 0.001 },
    { q: 0.1126, label: "√5", observed: false },
    { q: 0.1234, label: "√6", observed: true, residual: 0.012 },
    { q: 0.1333, label: "√7", observed: false },
  ],
};
export const IM3M: CombSeries = {
  phase: "Im3m", color: phaseColor("Im3m"), latticeLabel: "a = 142 Å", rSquared: 0.991,
  teeth: [
    { q: 0.0626, label: "√2", observed: true, residual: -0.003 },
    { q: 0.0885, label: "√4", observed: true, residual: 0.007 },
    { q: 0.1084, label: "√6", observed: false },
    { q: 0.1252, label: "√8", observed: true, residual: 0.041 },
  ],
};
export const LEFTOVER: number[] = [0.156, 0.205];
