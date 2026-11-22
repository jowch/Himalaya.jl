export interface Violation {
  rule: string;
  file: string;
  line: number;
  text: string;
}

export function scanContent(relPath: string, content: string): Violation[];
export function scanLegacyImports(relPath: string, content: string): Violation[];
export function scanUnpairedFontSize(relPath: string, content: string): Violation[];
