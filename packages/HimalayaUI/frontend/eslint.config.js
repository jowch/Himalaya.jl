import tseslint from "typescript-eslint";
import reactHooks from "eslint-plugin-react-hooks";
import uiCopy from "./scripts/eslint-ui-copy.mjs";

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
    },
  },
  {
    files: ["src/**/*.{ts,tsx}"],
    ignores: ["src/**/*.stories.tsx"],
    plugins: { copy: uiCopy },
    rules: { "copy/ui-copy": "error" },
  },
  {
    files: ["src/**/*.stories.tsx"],
    rules: { "no-console": "off" },
  },
);
