/**
 * comparePath helper — deep-link URL builder for the Compare page.
 *
 * Plan §Phase 4 follow-up: the global `/compare/all/:id` deep-link routes
 * were missing in the original PR. This helper centralizes URL shape so
 * the six navigation sites (sidebar, edit/fork buttons, conflict modal,
 * lineage badge, forks popover, warm-add menu, needs-review badge) stay
 * in sync.
 *
 * Compare UX Phase B-2: the `/edit` URL segment is gone — `comparePath`
 * no longer accepts an `edit` option. These cases lock the bare-URL shape.
 */
import { describe, it, expect } from "vitest";
import { comparePath } from "../src/lib/comparison/routes";

describe("comparePath", () => {
  describe("scope=experiment", () => {
    it("returns the list path for scope+eid only", () => {
      expect(comparePath({ scope: "experiment", eid: 7 })).toBe(
        "/experiments/7/compare",
      );
    });

    it("returns the review path for scope+eid+id", () => {
      expect(comparePath({ scope: "experiment", eid: 7, id: 42 })).toBe(
        "/experiments/7/compare/42",
      );
    });

    it("returns the create path for scope+eid+isNew", () => {
      expect(comparePath({ scope: "experiment", eid: 7, isNew: true })).toBe(
        "/experiments/7/compare/new",
      );
    });
  });

  describe("scope=all", () => {
    it("returns /compare/all for scope only", () => {
      expect(comparePath({ scope: "all" })).toBe("/compare/all");
    });

    it("returns /compare/all/:id for scope+id (the deep-link gap this fixes)", () => {
      expect(comparePath({ scope: "all", id: 99 })).toBe("/compare/all/99");
    });

    it("returns /compare/all/new for scope+isNew", () => {
      expect(comparePath({ scope: "all", isNew: true })).toBe(
        "/compare/all/new",
      );
    });

    it("ignores eid when scope=all (global is global)", () => {
      expect(comparePath({ scope: "all", eid: 7, id: 99 })).toBe(
        "/compare/all/99",
      );
    });
  });

  describe("defensive fallbacks", () => {
    it("falls through to /compare/all when scope=experiment but eid is undefined", () => {
      // TypeScript can't always prove eid is set; the runtime fallback
      // keeps the user on a coherent path rather than `/experiments//compare`.
      expect(comparePath({ scope: "experiment", eid: undefined, id: 5 })).toBe(
        "/compare/all/5",
      );
    });
  });
});
