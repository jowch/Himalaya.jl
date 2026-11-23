export interface ContractRow {
  id: string;
  behaviour: string;
  status: string;
  evidence: string;
  decided?: string;
  line: number;
}

export interface Contract {
  version: string | undefined;
  sha256: string | undefined;
  rows: ContractRow[];
}

export interface ContractDoc {
  name: string;
  pinned: boolean;
  dated: boolean;
  root: string;
}

export const MOCKUP: ContractDoc;
export const DECISIONS: ContractDoc;
export function parseContract(text: string): Contract;
export function checkContract(contract: Contract, read: (path: string) => string | undefined, doc?: ContractDoc): string[];
