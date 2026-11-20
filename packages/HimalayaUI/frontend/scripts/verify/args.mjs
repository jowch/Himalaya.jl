export class UsageError extends Error {
  constructor(message, hint) {
    super(message);
    this.name = "UsageError";
    this.hint = hint;
  }
}

export function parseArgs(argv, spec = {}) {
  const flags = {};
  const positional = [];
  const booleans = new Set(spec.booleans ?? []);
  const strings = new Set(spec.strings ?? []);
  for (let i = 0; i < argv.length; i++) {
    const arg = argv[i];
    if (arg === "--") {
      positional.push(...argv.slice(i + 1));
      break;
    }
    if (!arg.startsWith("--")) {
      positional.push(arg);
      continue;
    }
    const eq = arg.indexOf("=");
    const name = eq === -1 ? arg.slice(2) : arg.slice(2, eq);
    if (booleans.has(name)) {
      if (eq !== -1) throw new UsageError(`--${name} takes no value`);
      flags[name] = true;
    } else if (strings.has(name)) {
      let value;
      if (eq !== -1) value = arg.slice(eq + 1);
      else if (i + 1 < argv.length) value = argv[++i];
      else throw new UsageError(`--${name} needs a value`);
      flags[name] = value;
    } else {
      const known = [...booleans, ...strings].map((f) => `--${f}`).join(", ") || "none";
      throw new UsageError(`unknown option --${name}`, `options for this command: ${known}`);
    }
  }
  return { flags, positional };
}

export function intFlag(flags, name, fallback) {
  if (flags[name] === undefined) return fallback;
  const n = Number(flags[name]);
  if (!Number.isInteger(n) || n < 0) throw new UsageError(`--${name} must be a non-negative integer, got ${flags[name]}`);
  return n;
}
