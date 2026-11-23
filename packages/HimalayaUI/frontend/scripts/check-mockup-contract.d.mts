export interface ContractRow {
  id: string;
  behaviour: string;
  status: string;
  evidence: string;
  line: number;
}

export interface Contract {
  version: string | undefined;
  sha256: string | undefined;
  rows: ContractRow[];
}

export function parseContract(text: string): Contract;
export function checkContract(contract: Contract, read: (path: string) => string | undefined): string[];
