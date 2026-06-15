// builderAdapters.ts — pure view-model adapters for the Series Builder page.
// No React, no side effects. Maps API shapes → composite prop shapes.

import type { SeriesMember, SeriesSample, CorpusSample } from "../../api";
import type { MemberDatum } from "../components/MemberList";
import { dominantPhase } from "../waterfall/waterfallModel";

// ── membersToMemberData ───────────────────────────────────────────────────────

/**
 * Map API `SeriesMember[]` → `MemberDatum[]` for `MemberList`.
 *
 * - `key`           = stable member id string
 * - `phases`        = all phases for the swatch + label row:
 *                     form_factor / null state → `[]`;
 *                     indexed with confirmed_phases → that list's phase names;
 *                     indexed with only confirmed_index → `[dominant]`.
 * - `variableValue` = label_override ?? exposure-id fallback label
 * - `dataLine`      = "a = <lattice_d> Å · q₁ <q₁> Å⁻¹" when available, else "".
 *                     Matches the tone shown in SeriesMemberRow's JSDoc.
 */
export function membersToMemberData(members: SeriesMember[]): MemberDatum[] {
  return members.map((m): MemberDatum => {
    const snap = m.snapshot;
    const dominant = dominantPhase(m);

    // phases: coexistence-aware list
    let phases: string[];
    if (dominant === null) {
      phases = [];
    } else {
      const cp = snap?.confirmed_phases;
      if (cp && cp.length > 0) {
        phases = cp.map((p) => p.phase);
      } else {
        phases = [dominant];
      }
    }

    // variableValue: label_override wins; fall back to the same label
    // toWaterfallRows uses so the two panels always agree.
    const variableValue: string =
      m.label_override !== null
        ? m.label_override
        : m.exposure_id != null
          ? `exp ${m.exposure_id}`
          : "";

    // dataLine: lattice + first indexed peak q₁, formatted mono.
    const confirmedIndex = snap?.confirmed_index ?? null;
    let dataLine: string;
    if (confirmedIndex !== null && confirmedIndex.lattice_d != null) {
      // Find q₁ — the smallest q among the confirmed peak ids' effective peaks.
      const peakById = new Map(
        (snap?.effective_peaks ?? []).map((p) => [p.id, p]),
      );
      const indexedQs = confirmedIndex.peak_ids
        .map((id) => peakById.get(id)?.q)
        .filter((q): q is number => q !== undefined);
      const q1 = indexedQs.length > 0 ? Math.min(...indexedQs) : null;
      const latticeStr = confirmedIndex.lattice_d.toFixed(0);
      dataLine =
        q1 !== null
          ? `a = ${latticeStr} Å · q₁ ${q1.toFixed(3)} Å⁻¹`
          : `a = ${latticeStr} Å`;
    } else {
      dataLine = "";
    }

    return { key: String(m.id), phases, variableValue, dataLine };
  });
}

// ── recipeRowView ─────────────────────────────────────────────────────────────

/** View-model for one editable recipe row (a `SeriesSample` in the builder rail). */
export interface RecipeRowView {
  name: string;
  position: number;
}

/**
 * Derive the display name + position for one recipe row.
 * `sampleNameById` is the caller's pre-built lookup (id → display_name ?? name ?? "Sample <id>").
 * Falls back to `"Sample <sample_id>"` when the id is absent from the map.
 */
export function recipeRowView(
  row: SeriesSample,
  sampleNameById: Record<number, string>,
): RecipeRowView {
  return {
    name: sampleNameById[row.sample_id] ?? `Sample ${row.sample_id}`,
    position: row.position,
  };
}

// ── addableSamples ────────────────────────────────────────────────────────────

/**
 * Corpus samples whose id is NOT already in `draftRecipe` — the add-sample
 * dropdown options. Preserves corpus order.
 */
export function addableSamples(
  corpusSamples: CorpusSample[],
  draftRecipe: SeriesSample[],
): CorpusSample[] {
  const alreadyIn = new Set(draftRecipe.map((r) => r.sample_id));
  return corpusSamples.filter((s) => !alreadyIn.has(s.id));
}

// ── legendPhasesOf ────────────────────────────────────────────────────────────

/**
 * Collect distinct phase names across all members for the plate legend.
 * Includes both the dominant phase and any coexistence phases from
 * `confirmed_phases`. Null dominants (form_factor / unindexed) are excluded.
 * Order is first-appearance order across the member list.
 */
export function legendPhasesOf(members: SeriesMember[]): string[] {
  const seen = new Set<string>();
  const out: string[] = [];
  for (const m of members) {
    const snap = m.snapshot;
    const dominant = dominantPhase(m);
    if (dominant === null) continue;
    // Coexistence-aware: collect all phases from confirmed_phases when present.
    const cp = snap?.confirmed_phases;
    const phases: string[] =
      cp && cp.length > 0 ? cp.map((p) => p.phase) : [dominant];
    for (const phase of phases) {
      if (!seen.has(phase)) {
        seen.add(phase);
        out.push(phase);
      }
    }
  }
  return out;
}
