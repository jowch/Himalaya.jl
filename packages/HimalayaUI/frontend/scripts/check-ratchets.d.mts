export type Counts = Record<string, number>;

export interface Measure {
  appearance?: Counts;
  comments?: Counts;
}

export interface Change {
  kind: "appearance" | "comments";
  file: string;
  was: number;
  now: number;
}

export function classNameLiterals(content: string): string[];
export function countAppearance(relPath: string, content: string): number;
export function countTsComments(content: string): number;
export function countJuliaComments(content: string): number;
export function compareToBaseline(baseline: Measure, current: Measure): { increased: Change[]; decreased: Change[] };
export function canUpdate(baseline: Measure, current: Measure, exists?: (file: string) => boolean): boolean;
export function countCssComments(content: string): number;
export function trackedSources(): string[];
export const TS_EXTENSIONS: readonly string[];
export const JULIA_EXTENSIONS: readonly string[];
export const CSS_EXTENSIONS: readonly string[];
export const EXCLUDED: readonly string[];
export const DIRECTIVE_MAX: number;
export function measure(): Measure;
