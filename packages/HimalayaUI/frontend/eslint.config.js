import tseslint from "typescript-eslint";
import reactHooks from "eslint-plugin-react-hooks";

export default tseslint.config(
  { ignores: ["dist/**", "dist-e2e/**", "storybook-static/**", "node_modules/**"] },
  {
    files: ["src/**/*.{ts,tsx}"],
    languageOptions: { parser: tseslint.parser },
    linterOptions: { reportUnusedDisableDirectives: "error" },
    plugins: { "react-hooks": reactHooks },
    rules: {
      "react-hooks/rules-of-hooks": "error",
      "react-hooks/exhaustive-deps": "error",
      "no-console": ["error", { allow: ["warn", "error"] }],
      "no-restricted-syntax": [
        "error",
        {
          selector: "JSXOpeningElement[name.name=/^([a-z].*|Link|NavLink)$/] > JSXAttribute[name.name='title']",
          message: "Native title tooltips are slow, unstyled and never show on keyboard focus. Wrap the element in <Tooltip content=…> from src/ui, or <Tooltip overflowOnly> when the text only repeats a truncated label.",
        },
      ],
    },
  },
  {
    files: ["src/ui/CompareBar.tsx"],
    rules: { "no-restricted-syntax": "off" },
  },
  {
    files: ["src/**/*.stories.tsx"],
    rules: { "no-console": "off" },
  },
);
