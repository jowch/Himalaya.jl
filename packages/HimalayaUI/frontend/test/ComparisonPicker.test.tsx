/**
 * ComparisonPicker (Plan §Phase 5, Task 5.2 frontend half).
 *
 * Multi-select exposure picker mounted at body root. Surfaces a search box,
 * experiment + tag filters, "Recently used" section, and a multi-select
 * `<ul role="listbox">` of `<ExposureListRow>`. Already-added rows are
 * locked so the user can't double-add.
 *
 * Tests cover:
 *   - filters (search, experiment, tag)
 *   - default sort (alphabetical by exposure name)
 *   - multi-select toggling
 *   - already-added locks
 *   - empty state
 *   - focus trap (Tab cycles), Esc closes
 *   - restored focus on close
 */
import { describe, it, expect, beforeEach, vi } from "vitest";
import { screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { ComparisonPicker } from "../src/components/ComparisonPicker";
import { useAppState } from "../src/state";
import * as api from "../src/api";

const EXPERIMENTS: api.Experiment[] = [
  {
    id: 1,
    name: "Exp Alpha",
    path: "/data/alpha",
    data_dir: "/data/alpha/data",
    analysis_dir: "/data/alpha/analysis",
    manifest_path: null,
    created_at: "2026-04-01",
    q_units: null,
  },
  {
    id: 2,
    name: "Exp Beta",
    path: "/data/beta",
    data_dir: "/data/beta/data",
    analysis_dir: "/data/beta/analysis",
    manifest_path: null,
    created_at: "2026-04-15",
    q_units: null,
  },
];

const SAMPLES_E1: api.Sample[] = [
  {
    id: 10,
    experiment_id: 1,
    label: "JC001",
    name: "Sample A1",
    notes: "DOPC + chol",
    tags: [{ id: 1, key: "lipid", value: "DOPC", source: "manifest" }],
  },
  {
    id: 11,
    experiment_id: 1,
    label: "JC002",
    name: "Sample A2",
    notes: "POPC control",
    tags: [{ id: 2, key: "lipid", value: "POPC", source: "manifest" }],
  },
];

const SAMPLES_E2: api.Sample[] = [
  {
    id: 20,
    experiment_id: 2,
    label: "JC100",
    name: "Sample B1",
    notes: "Buffer",
    tags: [{ id: 3, key: "control", value: "buffer", source: "manifest" }],
  },
];

const EXPOSURES_E1_S10: api.Exposure[] = [
  exp(101, 10, "JC001-101"),
  exp(102, 10, "JC001-102"),
];
const EXPOSURES_E1_S11: api.Exposure[] = [exp(110, 11, "JC002-110")];
const EXPOSURES_E2_S20: api.Exposure[] = [exp(201, 20, "JC100-201")];

function exp(id: number, sample_id: number, name: string): api.Exposure {
  return {
    id,
    sample_id,
    filename: `${name}.dat`,
    kind: "file",
    selected: false,
    status: "accepted",
    image_path: null,
    image_version: "",
    tags: [],
    sources: [],
    trace_hash: null,
    analysis_inputs_hash: null,
  };
}

function resetStore(): void {
  localStorage.clear();
  sessionStorage.clear();
  useAppState.setState({
    username: "alice",
    activeExperimentId: 1,
    activeDraft: {
      id: undefined,
      baseHash: undefined,
      title: "",
      description: "",
      members: [],
    },
  });
}

beforeEach(() => {
  resetStore();
  vi.spyOn(api, "listExperiments").mockResolvedValue(EXPERIMENTS);
  vi.spyOn(api, "listSamples").mockImplementation((expId: number) =>
    Promise.resolve(expId === 1 ? SAMPLES_E1 : SAMPLES_E2),
  );
  vi.spyOn(api, "listExposures").mockImplementation((sampleId: number) => {
    if (sampleId === 10) return Promise.resolve(EXPOSURES_E1_S10);
    if (sampleId === 11) return Promise.resolve(EXPOSURES_E1_S11);
    if (sampleId === 20) return Promise.resolve(EXPOSURES_E2_S20);
    return Promise.resolve([]);
  });
  vi.spyOn(api, "getRecentlyPickedExposures").mockResolvedValue([]);
  vi.spyOn(api, "getSampleTags").mockResolvedValue([
    { key: "lipid", value: "DOPC" },
    { key: "lipid", value: "POPC" },
  ]);
});

describe("<ComparisonPicker>", () => {
  it("renders nothing when isOpen is false", () => {
    const { container } = renderWithProviders(
      <ComparisonPicker
        isOpen={false}
        onClose={() => {}}
        experimentId={1}
      />,
    );
    expect(container.querySelector('[data-testid="comparison-picker"]')).toBeNull();
  });

  it("renders dialog with role + aria-label when isOpen", async () => {
    renderWithProviders(
      <ComparisonPicker isOpen onClose={() => {}} experimentId={1} />,
    );
    const dialog = await screen.findByTestId("comparison-picker");
    expect(dialog.getAttribute("role")).toBe("dialog");
    expect(dialog).toHaveAttribute("aria-modal", "true");
  });

  it("filters by search across exposure name + sample name + notes", async () => {
    const user = userEvent.setup();
    renderWithProviders(
      <ComparisonPicker isOpen onClose={() => {}} experimentId={1} />,
    );
    // Both initial rows present
    await screen.findByText("JC001-101");
    await screen.findByText("JC001-102");
    await screen.findByText("JC002-110");

    await user.type(screen.getByTestId("comparison-picker-search"), "POPC");
    // "POPC" is in JC002-110's sample notes — should remain visible.
    await waitFor(() => {
      expect(screen.queryByText("JC001-101")).toBeNull();
      expect(screen.queryByText("JC001-102")).toBeNull();
    });
    expect(screen.getByText("JC002-110")).toBeInTheDocument();
  });

  it("rows are sorted alphabetically by exposure name (default)", async () => {
    renderWithProviders(
      <ComparisonPicker isOpen onClose={() => {}} experimentId={1} />,
    );
    const rows = await screen.findAllByTestId("picker-row");
    const names = rows.map((r) => r.querySelector("span")?.textContent ?? "");
    // Filtered to exposure-name first column; alphabetical:
    // JC001-101, JC001-102, JC002-110
    expect(names[0]).toBe("JC001-101");
    expect(names[1]).toBe("JC001-102");
    expect(names[2]).toBe("JC002-110");
  });

  it("multi-select: clicking checkboxes toggles selection state", async () => {
    const user = userEvent.setup();
    renderWithProviders(
      <ComparisonPicker isOpen onClose={() => {}} experimentId={1} />,
    );
    const rows = await screen.findAllByTestId("picker-row");
    expect(rows.length).toBeGreaterThanOrEqual(2);
    const cb1 = rows[0]!.querySelector(
      '[data-testid="exposure-list-row-checkbox"]',
    ) as HTMLInputElement;
    const cb2 = rows[1]!.querySelector(
      '[data-testid="exposure-list-row-checkbox"]',
    ) as HTMLInputElement;
    expect(cb1.checked).toBe(false);
    expect(cb2.checked).toBe(false);
    await user.click(cb1);
    await user.click(cb2);
    expect(cb1.checked).toBe(true);
    expect(cb2.checked).toBe(true);

    // Add count surfaces in the action button label.
    const addBtn = screen.getByTestId("comparison-picker-add");
    expect(addBtn).toHaveTextContent(/Add 2 selected/i);
  });

  it("already-added rows render locked and cannot be toggled", async () => {
    const user = userEvent.setup();
    // Pre-populate the active draft with exposure 101 — picker should lock that row.
    useAppState.setState({
      activeDraft: {
        id: undefined,
        baseHash: undefined,
        title: "",
        description: "",
        members: [
          {
            id: undefined,
            exposure_id: 101,
            display_order: 0,
            band_height: 1,
            y_offset: 0,
            normalization: "none",
            color_override: undefined,
            label_override: undefined,
            q_window_min: undefined,
            q_window_max: undefined,
            peak_display: undefined,
            snapshot: undefined,
          },
        ],
      },
    });
    renderWithProviders(
      <ComparisonPicker isOpen onClose={() => {}} experimentId={1} />,
    );
    const lockedRow = (await screen.findAllByTestId("picker-row")).find(
      (r) => r.getAttribute("data-exposure-id") === "101",
    );
    expect(lockedRow).toBeTruthy();
    expect(lockedRow!.getAttribute("data-locked")).toBe("true");
    const cb = lockedRow!.querySelector(
      '[data-testid="exposure-list-row-checkbox"]',
    ) as HTMLInputElement;
    expect(cb.disabled).toBe(true);

    await user.click(cb);
    // Disabled checkbox doesn't toggle in jsdom; assert state.
    expect(cb.checked).toBe(true); // already-added shows as checked
  });

  it("clicking 'Add selected' invokes addMember for each picked exposure and closes", async () => {
    const user = userEvent.setup();
    const onClose = vi.fn();
    renderWithProviders(
      <ComparisonPicker isOpen onClose={onClose} experimentId={1} />,
    );
    const rows = await screen.findAllByTestId("picker-row");
    const cb1 = rows[0]!.querySelector(
      '[data-testid="exposure-list-row-checkbox"]',
    ) as HTMLInputElement;
    const cb2 = rows[1]!.querySelector(
      '[data-testid="exposure-list-row-checkbox"]',
    ) as HTMLInputElement;
    await user.click(cb1);
    await user.click(cb2);

    await user.click(screen.getByTestId("comparison-picker-add"));

    await waitFor(() => {
      expect(onClose).toHaveBeenCalledTimes(1);
    });
    const draft = useAppState.getState().activeDraft;
    expect(draft).not.toBeNull();
    expect(draft!.members).toHaveLength(2);
  });

  it("'Add selected' is disabled when nothing is selected", async () => {
    renderWithProviders(
      <ComparisonPicker isOpen onClose={() => {}} experimentId={1} />,
    );
    const addBtn = await screen.findByTestId("comparison-picker-add");
    expect((addBtn as HTMLButtonElement).disabled).toBe(true);
  });

  it("empty state surfaces a hint when filters yield zero rows", async () => {
    const user = userEvent.setup();
    renderWithProviders(
      <ComparisonPicker isOpen onClose={() => {}} experimentId={1} />,
    );
    await screen.findByText("JC001-101");
    await user.type(screen.getByTestId("comparison-picker-search"), "ZZZ-no-match");
    expect(await screen.findByTestId("comparison-picker-empty")).toBeInTheDocument();
  });

  it("Esc key closes the modal", async () => {
    const user = userEvent.setup();
    const onClose = vi.fn();
    renderWithProviders(
      <ComparisonPicker isOpen onClose={onClose} experimentId={1} />,
    );
    await screen.findByTestId("comparison-picker");
    await user.keyboard("{Escape}");
    await waitFor(() => expect(onClose).toHaveBeenCalled());
  });

  it("focus trap: Tab cycles within modal (last → first)", async () => {
    const user = userEvent.setup();
    renderWithProviders(
      <ComparisonPicker isOpen onClose={() => {}} experimentId={1} />,
    );
    const dialog = await screen.findByTestId("comparison-picker");
    const focusable = dialog.querySelectorAll(
      'button:not([disabled]),input:not([disabled]),[tabindex]:not([tabindex="-1"])',
    );
    expect(focusable.length).toBeGreaterThan(1);
    const first = focusable[0] as HTMLElement;
    const last = focusable[focusable.length - 1] as HTMLElement;
    last.focus();
    expect(document.activeElement).toBe(last);
    await user.keyboard("{Tab}");
    expect(document.activeElement).toBe(first);
  });

  it("focus trap: Shift+Tab from first cycles to last", async () => {
    const user = userEvent.setup();
    renderWithProviders(
      <ComparisonPicker isOpen onClose={() => {}} experimentId={1} />,
    );
    const dialog = await screen.findByTestId("comparison-picker");
    const focusable = dialog.querySelectorAll(
      'button:not([disabled]),input:not([disabled]),[tabindex]:not([tabindex="-1"])',
    );
    const first = focusable[0] as HTMLElement;
    const last = focusable[focusable.length - 1] as HTMLElement;
    first.focus();
    await user.keyboard("{Shift>}{Tab}{/Shift}");
    expect(document.activeElement).toBe(last);
  });

  it("restores focus to the previously focused element on close", async () => {
    const user = userEvent.setup();
    function Harness(): JSX.Element {
      const [open, setOpen] = (require as any)("react").useState(false);
      return (
        <div>
          <button data-testid="trigger" onClick={() => setOpen(true)}>
            open
          </button>
          <ComparisonPicker
            isOpen={open}
            onClose={() => setOpen(false)}
            experimentId={1}
          />
        </div>
      );
    }
    renderWithProviders(<Harness />);
    const trigger = screen.getByTestId("trigger");
    trigger.focus();
    expect(document.activeElement).toBe(trigger);
    await user.click(trigger);
    await screen.findByTestId("comparison-picker");
    await user.keyboard("{Escape}");
    await waitFor(() => {
      expect(document.activeElement).toBe(trigger);
    });
  });

  it("renders the experiment filter chip seeded to the current experiment", async () => {
    renderWithProviders(
      <ComparisonPicker isOpen onClose={() => {}} experimentId={1} />,
    );
    const chip = await screen.findByTestId("picker-filter-experiment");
    expect(chip).toBeInTheDocument();
  });

  it("tag filter: selecting a tag pair filters the list to matching samples", async () => {
    const user = userEvent.setup();
    renderWithProviders(
      <ComparisonPicker isOpen onClose={() => {}} experimentId={1} />,
    );
    // Wait for tags to load.
    await screen.findByTestId("picker-filter-tag");

    // Select the lipid:DOPC tag chip — selects only Sample A1 (id 10),
    // hiding A2's exposure (110) which carries lipid:POPC.
    const dopcChip = await within(
      screen.getByTestId("picker-filter-tag"),
    ).findByTestId("picker-tag-option-lipid:DOPC");
    await user.click(dopcChip);

    await waitFor(() => {
      expect(screen.queryByText("JC002-110")).toBeNull();
    });
    expect(screen.getByText("JC001-101")).toBeInTheDocument();
    expect(screen.getByText("JC001-102")).toBeInTheDocument();
  });

  it("Recently used section shows when the user has prior picks", async () => {
    vi.spyOn(api, "getRecentlyPickedExposures").mockResolvedValue([102, 101]);
    // Make sure listUsers can resolve "alice" → an id.
    vi.spyOn(api, "listUsers").mockResolvedValue([
      { id: 7, username: "alice", first_name: null, last_name: null },
    ]);
    renderWithProviders(
      <ComparisonPicker isOpen onClose={() => {}} experimentId={1} />,
    );
    const recents = await screen.findByTestId("comparison-picker-recents");
    expect(within(recents).getByText("JC001-102")).toBeInTheDocument();
    expect(within(recents).getByText("JC001-101")).toBeInTheDocument();
  });
});
