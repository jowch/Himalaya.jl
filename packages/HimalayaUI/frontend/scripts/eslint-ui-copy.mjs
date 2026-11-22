const COPY_ATTRS = new Set([
  "title", "aria-label", "aria-description", "aria-roledescription", "placeholder", "alt",
  "label", "content", "hint", "description", "message", "what", "gesture", "heading", "note", "tail",
]);
const COPY_KEYS = new Set(["what", "gesture", "hint", "label", "title", "description", "message", "placeholder", "content", "note"]);

const VOCABULARY = [
  { re: /(?<!\btotal\s)\bexposures?\b(?!\s+times?\b|\s*\(s\))/i, word: "exposure", use: "image", why: "The UI calls one detector image of a sample an \"image\". \"Exposure\" stays only for time: \"exposure time\", \"Exposure (s)\", \"total exposure\"." },
  { re: /\bframes?\b/i, word: "frame", use: "image", why: "The UI calls one detector image of a sample an \"image\"." },
  { re: /\btag(s|ged|ging)?\b/i, word: "tag", use: "label", why: "Tags were renamed labels." },
  { re: /\bmut(e|es|ed|ing)\b/i, word: "mute", use: "Drop", why: "Use the plain literal verb: an image is dropped or restored." },
];

const HINT_VERBS = [
  "add", "remove", "hide", "show", "open", "close", "cancel", "confirm", "subtract", "toggle", "drop", "restore",
  "select", "move", "zoom", "reset", "commit", "clear", "apply", "create", "delete", "undo", "redo", "enlarge",
  "nudge", "compare", "go", "jump", "expand", "collapse", "rename", "accept", "dismiss",
];
const VERB_FORMS = HINT_VERBS.map((v) => (/(sh|ch|x|ss|z|o)$/.test(v) ? `${v}es` : /[^aeiou]y$/.test(v) ? `${v.slice(0, -1)}ies` : `${v}s`));
const GESTURE = new RegExp(`(^|[\\s(,;:])(\\S*[Cc]lick|Enter|Esc|Escape|Return|Space|Tab|Drag|Scroll|Backspace)\\s+(${VERB_FORMS.join("|")})\\b`);
const ABBREVIATIONS = /\b(e\.g|i\.e|vs|etc|cf|approx|no)\.(?=\s)/gi;
const SENTENCE_BREAK = /[.!?]\s+\S/;
const INSTRUCTION = /^(\S*[Cc]lick|Press|Drag|Scroll|Hold)\b.*\bto\b/;
const PROSE = /^[A-Z⌘][a-z]*\.?\s+[a-z]/;
const NOT_COPY_ATTRS = /^(className|class|style|data-|id$|key$|to$|href$|src$|role$|type$|name$|d$|viewBox$|fill$|stroke)/;
const TITLE_CASE = /^[A-Z][a-z]+(\s+[A-Z][a-z]+){2,}$/;

function stem(verb) {
  if (/ies$/.test(verb)) return `${verb.slice(0, -3)}y`;
  if (/(sh|ch|x|ss|z|o)es$/.test(verb)) return verb.slice(0, -2);
  return verb.slice(0, -1);
}

function fixFor(gesture, verb) {
  const v = stem(verb);
  return /click/i.test(gesture) ? `${gesture} to ${v}` : /^(Drag|Scroll)$/.test(gesture) ? `${gesture} to ${v}` : `Press ${gesture} to ${v}`;
}

function textOf(node) {
  if (node.type === "Literal" && typeof node.value === "string") return { text: node.value, template: false };
  if (node.type === "TemplateLiteral") return { text: node.quasis.map((q) => q.value.cooked ?? "").join("${}"), template: node.expressions.length > 0 };
  return null;
}

function collect(expr, out) {
  if (!expr) return out;
  if (expr.type === "JSXExpressionContainer") return collect(expr.expression, out);
  if (expr.type === "ConditionalExpression") { collect(expr.consequent, out); return collect(expr.alternate, out); }
  if (expr.type === "LogicalExpression") { collect(expr.left, out); return collect(expr.right, out); }
  if (expr.type === "TSAsExpression" || expr.type === "TSSatisfiesExpression") return collect(expr.expression, out);
  const t = textOf(expr);
  if (t) out.push({ node: expr, ...t });
  if (expr.type === "TemplateLiteral") for (const e of expr.expressions) collect(e, out);
  return out;
}

function attrName(attr) {
  return attr.name.type === "JSXNamespacedName" ? `${attr.name.namespace.name}:${attr.name.name.name}` : attr.name.name;
}

function keyName(prop) {
  if (prop.computed) return null;
  if (prop.key.type === "Identifier") return prop.key.name;
  if (prop.key.type === "Literal") return String(prop.key.value);
  return null;
}

function isSoleText(node) {
  const kids = node.parent?.children ?? [];
  return kids.filter((k) => !(k.type === "JSXText" && k.value.trim() === "")).length === 1;
}

function owner(node) {
  let n = node.parent;
  while (n && ["ConditionalExpression", "LogicalExpression", "JSXExpressionContainer", "TSAsExpression", "TSSatisfiesExpression", "TemplateLiteral"].includes(n.type)) n = n.parent;
  return n;
}

function handled(node) {
  const o = owner(node);
  if (!o) return true;
  if (o.type === "JSXAttribute") return true;
  if (o.type === "Property" && o.value !== node && o.key === node) return true;
  if (o.type === "Property" && COPY_KEYS.has(keyName(o) ?? "")) return true;
  if (["ImportDeclaration", "ExportNamedDeclaration", "ExportAllDeclaration", "TSLiteralType", "TSExternalModuleReference"].includes(o.type)) return true;
  return false;
}

function flatText(node) {
  if (node.type === "JSXText") return node.value;
  if (node.type === "JSXExpressionContainer") {
    const t = node.expression && textOf(node.expression);
    return t ? t.text : "${}";
  }
  if (node.type === "JSXElement" || node.type === "JSXFragment") return node.children.map(flatText).join("");
  return "";
}

function hasElementChild(node) {
  return node.children.some((c) => c.type === "JSXElement");
}

function check(context, node, raw, { template, whole }) {
  const text = raw.replace(/\s+/g, " ").trim();
  if (!/[A-Za-z]{2}/.test(text)) return;

  const g = GESTURE.exec(text);
  if (g) {
    const gesture = g[2].replace("${}", "⌘");
    context.report({ node, messageId: "telegraphic", data: { found: `${g[2]} ${g[3]}`, fix: fixFor(gesture, g[3]) } });
  }

  if (whole) {
    const plain = text.replace(ABBREVIATIONS, "$1");
    const sentence = SENTENCE_BREAK.test(plain) || INSTRUCTION.test(plain);
    if (sentence && !/[.!?…:)]$/.test(plain)) context.report({ node, messageId: "period", data: { text } });
  }

  for (const v of VOCABULARY) {
    if (v.re.test(text)) context.report({ node, messageId: "vocabulary", data: { word: v.word, use: v.use, why: v.why } });
  }

  if (!template && TITLE_CASE.test(text)) context.report({ node, messageId: "titleCase", data: { text } });
}

const rule = {
  meta: {
    type: "suggestion",
    docs: { description: "UI copy style: full imperative hints, sentence case, closing periods and house vocabulary. See src/AGENTS.md, UI copy." },
    schema: [],
    messages: {
      telegraphic: "\"{{found}}\" is telegraphic. Write a full imperative: \"{{fix}}\". A hint says what to do, not what the control does (\"Click to hide.\", \"Press Enter to add it.\").",
      period: "\"{{text}}\" is a sentence, so it ends with a period. Only name-like labels and tooltips (\"Keyboard shortcuts\") go without one.",
      vocabulary: "User-facing copy says \"{{use}}\", not \"{{word}}\". {{why}}",
      titleCase: "\"{{text}}\" is in Title Case. UI copy uses sentence case: capitalise only the first word and proper names.",
    },
  },
  create(context) {
    return {
      JSXElement(node) {
        if (!hasElementChild(node)) return;
        const text = flatText(node).replace(/\s+/g, " ").trim();
        const g = GESTURE.exec(text);
        if (g && !node.children.some((c) => c.type === "JSXText" && GESTURE.test(c.value.replace(/\s+/g, " ").trim()))) {
          const inner = node.children.filter((c) => c.type === "JSXElement").some((c) => GESTURE.test(flatText(c).replace(/\s+/g, " ").trim()));
          if (!inner) context.report({ node, messageId: "telegraphic", data: { found: `${g[2]} ${g[3]}`, fix: fixFor(g[2].replace("${}", "⌘"), g[3]) } });
        }
      },
      JSXText(node) {
        check(context, node, node.value, { template: false, whole: isSoleText(node) });
      },
      JSXAttribute(node) {
        const name = attrName(node);
        if (NOT_COPY_ATTRS.test(name)) return;
        const strict = COPY_ATTRS.has(name);
        for (const s of collect(node.value, [])) {
          if (strict || PROSE.test(s.text)) check(context, s.node, s.text, { template: s.template, whole: name !== "placeholder" });
        }
      },
      Literal(node) {
        if (typeof node.value !== "string" || !PROSE.test(node.value) || handled(node)) return;
        check(context, node, node.value, { template: false, whole: true });
      },
      TemplateLiteral(node) {
        const t = textOf(node);
        if (!t || !PROSE.test(t.text) || handled(node)) return;
        check(context, node, t.text, { template: t.template, whole: true });
      },
      Property(node) {
        const k = keyName(node);
        if (!k || !COPY_KEYS.has(k) || node.parent?.type !== "ObjectExpression") return;
        for (const s of collect(node.value, [])) check(context, s.node, s.text, { template: s.template, whole: true });
      },
    };
  },
};

export default { meta: { name: "himalaya-ui-copy" }, rules: { "ui-copy": rule } };
