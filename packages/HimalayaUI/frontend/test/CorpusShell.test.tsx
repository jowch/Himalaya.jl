import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { CorpusShell } from "../src/components/CorpusShell";
import { SamplesPage } from "../src/pages/SamplesPage";

describe("CorpusShell", () => {
  it("renders the topbar and the matched child route in its Outlet", () => {
    render(
      <MemoryRouter initialEntries={["/samples"]}>
        <Routes>
          <Route element={<CorpusShell />}>
            <Route path="/samples" element={<SamplesPage />} />
          </Route>
        </Routes>
      </MemoryRouter>,
    );
    expect(screen.getByTestId("corpus-shell")).toBeInTheDocument();
    expect(screen.getByTestId("corpus-topbar")).toBeInTheDocument();
    expect(screen.getByTestId("samples-page")).toBeInTheDocument();
  });
});

describe("SamplesPage", () => {
  it("renders the placeholder body", () => {
    render(<SamplesPage />);
    expect(screen.getByTestId("samples-page")).toBeInTheDocument();
  });
});
