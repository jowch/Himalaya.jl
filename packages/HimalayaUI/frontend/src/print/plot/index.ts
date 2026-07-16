export { TracePlot } from "./TracePlot";
export type {
  TracePlotProps,
  TraceModel,
  PlotContext,
  TracePlotInteraction,
} from "./TracePlot";
export type { PlotPeak } from "./marks/PlotPeaks";
export { TraceLine } from "./marks/TraceLine";
export { PlotPeaks } from "./marks/PlotPeaks";
export { Axis } from "./Axis";
export { PlotFrame } from "./PlotFrame";
export type { Margins, PlotDims } from "./PlotFrame";
export {
  makeProjection,
  makeAxis,
  positiveExtent,
  type Projection,
  type Axis1D,
  type ScaleType,
} from "./projection";
export {
  hitTestPeaks,
  zoomXDomain,
  PEAK_HIT_PX,
  BRUSH_DEADZONE_PX,
  MIN_SPAN_FRAC,
} from "./interaction";
