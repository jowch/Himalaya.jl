// types.ts — export spec passed from adapters to the renderer.
import * as Plot from "@observablehq/plot";

export interface LegendRow {
  /** CSS color used for the swatch / line stroke. */
  color: string;
  /** Mark style — small filled square for trace key, line for predicted-q phases. */
  symbol: "swatch" | "line";
  label: string;
}

export interface LegendSpec {
  rows: LegendRow[];
}

export interface ExportSpec {
  title: { primary: string; secondary?: string };
  width: number;   // logical (CSS) pixels
  height: number;
  /** Marks + x/y configs handed to Plot.plot(). Adapters MUST NOT set
   *  `title`, `caption`, or `figure: true` — the renderer wraps the SVG
   *  itself, and those fields make Plot return HTMLElement instead of
   *  SVGSVGElement (breaks the rest of the pipeline). */
  plot: Plot.PlotOptions;
  legend?: LegendSpec;
  /** Override the figure font family (Plan E E-8 clean preset → Arial). When
   *  omitted the renderer uses its default Print sans-serif. */
  fontFamily?: string;
  /** A centered footnote line under the plot (Plan E E-8 clean preset). */
  footnote?: string;
}
