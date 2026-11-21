import { DETECTOR_LUT_VARIANTS, detectorLutFilterId, detectorLutTableValues } from "./detectorLut";

const TABLES = DETECTOR_LUT_VARIANTS.map((variant) => ({ id: detectorLutFilterId(variant), ...detectorLutTableValues(variant) }));

export function DetectorLutFilters(): JSX.Element {
  return (
    <svg
      aria-hidden="true"
      focusable="false"
      style={{ position: "absolute", width: 0, height: 0 }}
    >
      <defs>
        {TABLES.map((t) => (
          <filter key={t.id} id={t.id} colorInterpolationFilters="sRGB">
            <feComponentTransfer>
              <feFuncR type="table" tableValues={t.r} />
              <feFuncG type="table" tableValues={t.g} />
              <feFuncB type="table" tableValues={t.b} />
            </feComponentTransfer>
          </filter>
        ))}
      </defs>
    </svg>
  );
}
