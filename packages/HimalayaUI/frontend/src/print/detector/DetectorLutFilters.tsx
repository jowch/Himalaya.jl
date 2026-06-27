import { detectorLutTableValues } from "./detectorLut";

// Built once at module load from the shared LUT source (buildPrintDetectorLut).
const NEUTRAL = detectorLutTableValues("neutral");
const WARM = detectorLutTableValues("warm");

/**
 * Shared SVG `<filter>` defs for the detector colormap, referenced by
 * `DetectorImage` via `filter: url(#detector-lut-<variant>)`. The browser
 * recolors the grayscale detector PNG on the GPU — no per-pixel JS loop, no
 * canvas readback, and switching `lutVariant` is a pure style swap (no refetch).
 *
 * MOUNT ONCE, app-wide (App.tsx) — and in the Storybook preview, so detector
 * stories resolve the same `#id`. A CSS `filter` that references a MISSING `#id`
 * hides the element in some browsers, so an unmounted defs set isn't just "no
 * colour" — it can blank the image. Same-document `#id` refs are universally
 * supported (data-URI filter refs are NOT, in Chrome — hence a real mount).
 *
 * `color-interpolation-filters="sRGB"` is load-bearing: the default `linearRGB`
 * would gamma-shift every value, since the LUT was authored in sRGB 8-bit space.
 */
export function DetectorLutFilters(): JSX.Element {
  return (
    <svg
      aria-hidden="true"
      focusable="false"
      style={{ position: "absolute", width: 0, height: 0 }}
    >
      <defs>
        <filter id="detector-lut-neutral" colorInterpolationFilters="sRGB">
          <feComponentTransfer>
            <feFuncR type="table" tableValues={NEUTRAL.r} />
            <feFuncG type="table" tableValues={NEUTRAL.g} />
            <feFuncB type="table" tableValues={NEUTRAL.b} />
          </feComponentTransfer>
        </filter>
        <filter id="detector-lut-warm" colorInterpolationFilters="sRGB">
          <feComponentTransfer>
            <feFuncR type="table" tableValues={WARM.r} />
            <feFuncG type="table" tableValues={WARM.g} />
            <feFuncB type="table" tableValues={WARM.b} />
          </feComponentTransfer>
        </filter>
      </defs>
    </svg>
  );
}
