/**
 * AUTO-GENERATED real fixtures — DO NOT hand-edit.
 *
 * Captured from the himalaya-devdata dev DB (exposures 65/66/67/93/37) and the
 * SSRL `tot_files`. Peak q-values, phases, and lattices are REAL; the two
 * empty-state members are SYNTHETIC (the pre-Plan-A dev DB predates the
 * form_factor / null assignment states). Measured traces live alongside in
 * `traces/<exposureId>.json`, thumbnails in `thumbs/<exposureId>.png`.
 *
 * Regenerate via tmp/generate.py against a mounted /Volumes/data SSRL share.
 */
import type { SeriesMember } from "../../api";

/** The real Sample-9 phase transition + two stress members, in display order:
 *  Ia3d -> (Im3m + Ia3d coexistence) -> Im3m, then a dense cubic + sparse lamellar. */
export const realMembers: SeriesMember[] = [
  {
    "id": 9000,
    "series_id": 9001,
    "exposure_id": 65,
    "display_order": 0,
    "band_height": 1.0,
    "y_offset": 0,
    "normalization": "none",
    "color_override": null,
    "label_override": "C06-5",
    "q_window_min": null,
    "q_window_max": null,
    "peak_display": null,
    "snapshot": {
      "effective_peaks": [
        {
          "id": 2658,
          "q": 0.01102,
          "intensity": 404.8658,
          "sharpness": 33.8794,
          "source": "auto"
        },
        {
          "id": 2659,
          "q": 0.01318,
          "intensity": 206.8743,
          "sharpness": 3.8251,
          "source": "auto"
        },
        {
          "id": 2660,
          "q": 0.01621,
          "intensity": 86.01688,
          "sharpness": 3.6161,
          "source": "auto"
        },
        {
          "id": 2661,
          "q": 0.04559,
          "intensity": 138.1744,
          "sharpness": 3.4818,
          "source": "auto"
        },
        {
          "id": 2662,
          "q": 0.05078,
          "intensity": 151.2155,
          "sharpness": 11.2778,
          "source": "auto"
        },
        {
          "id": 2663,
          "q": 0.07195,
          "intensity": 118.0144,
          "sharpness": 5.1362,
          "source": "auto"
        },
        {
          "id": 2664,
          "q": 0.08794,
          "intensity": 138.0274,
          "sharpness": 7.4163,
          "source": "auto"
        }
      ],
      "confirmed_index": {
        "id": 761,
        "phase": "Ia3d",
        "lattice_d": 335.284,
        "r_squared": 0.99268,
        "ngc": 4,
        "peak_ids": [
          2661,
          2662,
          2663,
          2664
        ]
      },
      "assignment_state": "indexed",
      "confirmed_phases": [
        {
          "phase": "Ia3d",
          "lattice_d": 335.284
        }
      ],
      "analysis_inputs_hash": "3db0041dba40c65d5c237d49ba0138d6cc5709c7e14fb0d14e2c185c13a0f182"
    },
    "is_stale": false,
    "created_by": null,
    "created_at": "2026-05-31T00:00:00Z"
  },
  {
    "id": 9001,
    "series_id": 9001,
    "exposure_id": 66,
    "display_order": 1,
    "band_height": 1.0,
    "y_offset": 0,
    "normalization": "none",
    "color_override": null,
    "label_override": "C06-6",
    "q_window_min": null,
    "q_window_max": null,
    "peak_display": null,
    "snapshot": {
      "effective_peaks": [
        {
          "id": 3233,
          "q": 0.04559,
          "intensity": 63.86096,
          "sharpness": 1.6306,
          "source": "auto"
        },
        {
          "id": 90047,
          "q": 0.04943,
          "intensity": null,
          "sharpness": 0.0,
          "source": "manual"
        },
        {
          "id": 3234,
          "q": 0.06417,
          "intensity": 44.72523,
          "sharpness": 0.5681,
          "source": "auto"
        },
        {
          "id": 90051,
          "q": 0.06768,
          "intensity": null,
          "sharpness": 0.0,
          "source": "manual"
        },
        {
          "id": 3235,
          "q": 0.07022,
          "intensity": 45.76375,
          "sharpness": 0.7366,
          "source": "auto"
        },
        {
          "id": 3236,
          "q": 0.07843,
          "intensity": 57.73309,
          "sharpness": 0.679,
          "source": "auto"
        },
        {
          "id": 3237,
          "q": 0.08535,
          "intensity": 51.75329,
          "sharpness": 0.6316,
          "source": "auto"
        },
        {
          "id": 90048,
          "q": 0.11134,
          "intensity": null,
          "sharpness": 0.0,
          "source": "manual"
        },
        {
          "id": 3238,
          "q": 0.12121,
          "intensity": 32.00312,
          "sharpness": 0.0847,
          "source": "auto"
        },
        {
          "id": 90049,
          "q": 0.13088,
          "intensity": null,
          "sharpness": 0.0,
          "source": "manual"
        },
        {
          "id": 90050,
          "q": 0.14936,
          "intensity": null,
          "sharpness": 0.0,
          "source": "manual"
        }
      ],
      "confirmed_index": {
        "id": 766,
        "phase": "Im3m",
        "lattice_d": 179.148,
        "r_squared": 0.99986,
        "ngc": 7,
        "peak_ids": [
          47,
          48,
          49,
          50,
          3235,
          3237,
          3238
        ]
      },
      "assignment_state": "indexed",
      "confirmed_phases": [
        {
          "phase": "Im3m",
          "lattice_d": 179.148
        },
        {
          "phase": "Ia3d",
          "lattice_d": 227.448
        },
        {
          "phase": "Im3m",
          "lattice_d": 194.135
        }
      ],
      "analysis_inputs_hash": "853d6fbb0904e8b97a53d3ed3f24757049af518d58c4bd9b0ca0ae0c52e22367"
    },
    "is_stale": false,
    "created_by": null,
    "created_at": "2026-05-31T00:00:00Z"
  },
  {
    "id": 9002,
    "series_id": 9001,
    "exposure_id": 67,
    "display_order": 2,
    "band_height": 1.0,
    "y_offset": 0,
    "normalization": "none",
    "color_override": null,
    "label_override": "C06-7",
    "q_window_min": null,
    "q_window_max": null,
    "peak_display": null,
    "snapshot": {
      "effective_peaks": [
        {
          "id": 2671,
          "q": 0.01102,
          "intensity": 63.63788,
          "sharpness": 0.8391,
          "source": "auto"
        },
        {
          "id": 2672,
          "q": 0.04516,
          "intensity": 61.05825,
          "sharpness": 2.0892,
          "source": "auto"
        },
        {
          "id": 2673,
          "q": 0.04948,
          "intensity": 61.56813,
          "sharpness": 3.0013,
          "source": "auto"
        },
        {
          "id": 2674,
          "q": 0.06417,
          "intensity": 43.46361,
          "sharpness": 0.463,
          "source": "auto"
        },
        {
          "id": 2675,
          "q": 0.06979,
          "intensity": 47.89203,
          "sharpness": 1.1687,
          "source": "auto"
        },
        {
          "id": 2676,
          "q": 0.07843,
          "intensity": 56.25599,
          "sharpness": 0.5875,
          "source": "auto"
        },
        {
          "id": 2677,
          "q": 0.08578,
          "intensity": 56.25473,
          "sharpness": 1.8563,
          "source": "auto"
        },
        {
          "id": 2678,
          "q": 0.12121,
          "intensity": 32.12883,
          "sharpness": 0.1365,
          "source": "auto"
        },
        {
          "id": 2679,
          "q": 0.13115,
          "intensity": 31.5565,
          "sharpness": 0.1209,
          "source": "auto"
        }
      ],
      "confirmed_index": {
        "id": 772,
        "phase": "Im3m",
        "lattice_d": 179.478,
        "r_squared": 0.99998,
        "ngc": 5,
        "peak_ids": [
          2673,
          2675,
          2677,
          2678,
          2679
        ]
      },
      "assignment_state": "indexed",
      "confirmed_phases": [
        {
          "phase": "Im3m",
          "lattice_d": 179.478
        }
      ],
      "analysis_inputs_hash": "27b8426320b879d1d2ab5b214a27974cab113d81c3178938af67a875c0d4c7e3"
    },
    "is_stale": false,
    "created_by": null,
    "created_at": "2026-05-31T00:00:00Z"
  },
  {
    "id": 9003,
    "series_id": 9001,
    "exposure_id": 93,
    "display_order": 3,
    "band_height": 1.0,
    "y_offset": 0,
    "normalization": "none",
    "color_override": null,
    "label_override": "C10-1",
    "q_window_min": null,
    "q_window_max": null,
    "peak_display": null,
    "snapshot": {
      "effective_peaks": [
        {
          "id": 2705,
          "q": 0.00584,
          "intensity": 1117.069,
          "sharpness": 43.079,
          "source": "auto"
        },
        {
          "id": 2706,
          "q": 0.01016,
          "intensity": 505.7177,
          "sharpness": 19.5963,
          "source": "auto"
        },
        {
          "id": 2707,
          "q": 0.01102,
          "intensity": 464.7073,
          "sharpness": 20.0226,
          "source": "auto"
        },
        {
          "id": 2708,
          "q": 0.01318,
          "intensity": 325.4491,
          "sharpness": 6.5233,
          "source": "auto"
        },
        {
          "id": 2709,
          "q": 0.01405,
          "intensity": 293.2132,
          "sharpness": 2.8122,
          "source": "auto"
        },
        {
          "id": 2710,
          "q": 0.01491,
          "intensity": 287.5872,
          "sharpness": 7.9864,
          "source": "auto"
        },
        {
          "id": 2711,
          "q": 0.01837,
          "intensity": 153.0102,
          "sharpness": 2.3586,
          "source": "auto"
        },
        {
          "id": 2712,
          "q": 0.01966,
          "intensity": 122.0701,
          "sharpness": 1.8966,
          "source": "auto"
        },
        {
          "id": 2713,
          "q": 0.02053,
          "intensity": 128.2819,
          "sharpness": 3.5831,
          "source": "auto"
        }
      ],
      "confirmed_index": {
        "id": 382,
        "phase": "Ia3d",
        "lattice_d": 1202.345,
        "r_squared": 0.93113,
        "ngc": 4,
        "peak_ids": [
          2708,
          2710,
          2712,
          2713
        ]
      },
      "assignment_state": "indexed",
      "confirmed_phases": [
        {
          "phase": "Ia3d",
          "lattice_d": 1202.345
        }
      ],
      "analysis_inputs_hash": "devdb-fixture"
    },
    "is_stale": false,
    "created_by": null,
    "created_at": "2026-05-31T00:00:00Z"
  },
  {
    "id": 9004,
    "series_id": 9001,
    "exposure_id": 37,
    "display_order": 4,
    "band_height": 1.0,
    "y_offset": 0,
    "normalization": "none",
    "color_override": null,
    "label_override": "C05-1",
    "q_window_min": null,
    "q_window_max": null,
    "peak_display": null,
    "snapshot": {
      "effective_peaks": [
        {
          "id": 2628,
          "q": 0.11733,
          "intensity": 934.1186,
          "sharpness": 168.3602,
          "source": "auto"
        },
        {
          "id": 2629,
          "q": 0.23486,
          "intensity": 31.54261,
          "sharpness": 2.1336,
          "source": "auto"
        }
      ],
      "confirmed_index": {
        "id": 747,
        "phase": "Lamellar",
        "lattice_d": 53.515,
        "r_squared": 1.0,
        "ngc": 2,
        "peak_ids": [
          2628,
          2629
        ]
      },
      "assignment_state": "indexed",
      "confirmed_phases": [
        {
          "phase": "Lamellar",
          "lattice_d": 53.515
        }
      ],
      "analysis_inputs_hash": "944953027ef9e00cca55cc532036f426bb9b161370ceb4eba4a2c2de0bedcafa"
    },
    "is_stale": false,
    "created_by": null,
    "created_at": "2026-05-31T00:00:00Z"
  }
];

/** Ia3d -> Im3m+Ia3d -> Im3m. The hero multi-trace series (exps 65/66/67). */
export const transitionSeries: SeriesMember[] = realMembers.slice(0, 3);

/** exp 66 — real repeated-phase 3-index coexistence group [Im3m, Ia3d, Im3m]. */
export const coexistenceMember: SeriesMember = realMembers[1]!;
/** exp 93 — dense Ia3d (9 peaks). */
export const denseCubicMember: SeriesMember = realMembers[3]!;
/** exp 37 — sparse Lamellar (2 peaks). */
export const sparseLamellarMember: SeriesMember = realMembers[4]!;

/** SYNTHETIC: a measured trace with auto peaks but no confirmed index. */
export const unindexedMember: SeriesMember = {
  "id": 9100,
  "series_id": 9001,
  "exposure_id": null,
  "display_order": 100,
  "band_height": 1.0,
  "y_offset": 0,
  "normalization": "none",
  "color_override": null,
  "label_override": "no-index",
  "q_window_min": null,
  "q_window_max": null,
  "peak_display": null,
  "snapshot": {
    "effective_peaks": [
      {
        "id": 99001,
        "q": 0.018,
        "intensity": null,
        "sharpness": 0,
        "source": "auto"
      },
      {
        "id": 99002,
        "q": 0.026,
        "intensity": null,
        "sharpness": 0,
        "source": "auto"
      },
      {
        "id": 99003,
        "q": 0.041,
        "intensity": null,
        "sharpness": 0,
        "source": "auto"
      }
    ],
    "confirmed_index": null,
    "confirmed_phases": [],
    "analysis_inputs_hash": "synthetic-unindexed"
  },
  "is_stale": false,
  "created_by": null,
  "created_at": "2026-05-31T00:00:00Z"
};
/** SYNTHETIC: a real trace with no Bragg peaks (form-factor only). */
export const formFactorMember: SeriesMember = {
  "id": 9101,
  "series_id": 9001,
  "exposure_id": null,
  "display_order": 101,
  "band_height": 1.0,
  "y_offset": 0,
  "normalization": "none",
  "color_override": null,
  "label_override": "form-factor",
  "q_window_min": null,
  "q_window_max": null,
  "peak_display": null,
  "snapshot": {
    "effective_peaks": [],
    "confirmed_index": null,
    "assignment_state": "form_factor",
    "confirmed_phases": [],
    "analysis_inputs_hash": "synthetic-form-factor"
  },
  "is_stale": false,
  "created_by": null,
  "created_at": "2026-05-31T00:00:00Z"
};

/** Every PhaseStrip cell state in one array (real indexed + synthetic empties). */
export const allStripStates: SeriesMember[] = [
  ...transitionSeries,
  denseCubicMember,
  sparseLamellarMember,
  formFactorMember,
  unindexedMember,
];
