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
  floor: number;
  supersededById?: boolean;
}

export const MOCKUP: ContractDoc;
export const DECISIONS: ContractDoc;
export function parseContract(text: string): Contract;
export function decisionRefs(evidence: string): string[];
export function deletedRows(base: Contract, contract: Contract): string[];
export function checkRowCount(contract: Contract, doc: ContractDoc): string[];
export function checkContract(contract: Contract, read: (path: string) => string | undefined, doc?: ContractDoc, decisionIds?: ReadonlySet<string>): string[];
