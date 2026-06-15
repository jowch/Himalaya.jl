# Product

## Register

product

## Users

Structural scientists working with small-angle X-ray scattering (SAXS) data: soft-matter chemists, structural biologists, and beamline researchers analyzing liquid-crystalline and lipid mesophases. They arrive with a beamtime's worth of 1D integration traces and a question, not a dataset to admire: *what phase is this, and can I trust the assignment?*

Their context is focused and iterative. They are not browsing; they are forming a hypothesis about one diffraction pattern, testing it against known phase-ratio series (Pn3m, Im3m, Ia3d, Fm3m, Fd3m, Hexagonal, Lamellar, Square), and often comparing several patterns side by side. They are domain experts who can overrule the software, and will, so the software has to show its work.

**Primary persona: Dr. Maya Okonkwo, structural scientist.** A soft-matter chemist back from a synchrotron beamtime with a few dozen 1D integration traces of a lipid mesophase series. She knows the phase-ratio series cold and trusts her own eye on a clean pattern, so she uses HimalayaUI to do the arithmetic she would otherwise do by hand: propose a phase, score the fit, and let her curate peaks and overrule the call in seconds. She wants the assignment, the reasoning behind it, and a publication-ready figure, with the least chrome between her and the data.

## Product Purpose

HimalayaUI is a thinking tool for forming and testing phase hypotheses about SAXS diffraction patterns, both deep within a single pattern and across several. It wraps the core Himalaya peak-finding and indexing engine in an interactive workspace: find Bragg peaks in a trace, fit their q-values to candidate phases, curate the result, and compare patterns.

Success is a scientist trusting an assignment they did not have to compute by hand, and being able to overrule it in seconds when their judgment disagrees. The product earns trust by being transparent (every score, ratio, and peak is inspectable) and fast (curation is a live edit, not a batch job).

## Brand Personality

Assured, precise, unshowy. The interface is a **confident expert**: it makes the call, proposes the best phase and the best index, and shows the reasoning behind it rather than dumping raw options and asking the scientist to sort them out. Confidence lives in the *judgments*, not the chrome: the UI itself is quiet, dense, and fast, so the diffraction data stays the subject. It speaks like a knowledgeable colleague, not a wizard and not a lab manual: direct, specific, never tentative, never cute.

## Anti-references

This product should not resemble any of its four common neighbors:

- **Legacy scientific software** (ImageJ, Origin, dated beamline tools): cluttered toolbars, grey system chrome, widget soup. HimalayaUI is modern and composed, not a pile of controls.
- **Generic SaaS dashboard** (the Stripe/Vercel template look): hero-metric tiles, identical card grids, decorative gradient accents. HimalayaUI's cards are few, deliberate, and never repeated as a grid.
- **Consumer / playful app**: rounded-everything, mascots, bright primary colors, gamified flourishes. The tone is professional, not friendly-for-its-own-sake.
- **Bare academic utility** (a Jupyter notebook shipped as a product): unstyled, function-only, no visual care. Craft is expected; it signals the analysis can be trusted.

The target sits in the gap between all four: visually composed and modern, but serious, dense, and instrument-like.

## Design Principles

1. **The tool makes the call.** Propose the best phase and the best index, ranked, with a clear primary choice. Never present raw possibilities and make the scientist do the sorting. A confident default that can be overruled beats a neutral menu.
2. **Show the reasoning.** This is a thinking tool, not an oracle. Every assertion is backed by inspectable evidence: the score, the q-ratios, the peak persistence. Trust comes from transparency, not from authority.
3. **Color is semantic, never decorative.** Hue carries meaning: peak provenance (auto vs. manual), phase identity, status. A color is a label. Decorative color is a bug, and every encoding must survive color blindness.
4. **Recede so the data leads.** The trace, the detector image, and the Miller plot are the subject. Chrome is quiet, dense, and quick. Animation is polish measured in milliseconds, never spectacle.
5. **Earned craft over template reflex.** Resist the nearest cliché (SaaS card grid, beamline toolbar). Extend the existing system (the type scale, the token set) rather than scattering one-offs. Consistency is how a tool reads as trustworthy.

## Accessibility & Inclusion

- **Color-blind safety is a hard requirement.** Phase and peak colors carry meaning, so every color encoding must be paired with a second channel: shape, label, position, or pattern. An index identified only by hue is broken. Verify encodings against deuteranopia and protanopia.
- **WCAG AA contrast** on all text and interactive elements. The product ships a single light "Print" identity (warm paper and ink); every text and interaction colour, including the eight phase hues, is tuned to clear AA on the `plate` surface. There is no dark theme to maintain a second contrast budget for.
- Honor `prefers-reduced-motion`. The system's motion is deliberately minimal: 120ms ease-out transitions on colour, background, border, and opacity, plus a few short functional animations (overlay fade-ins, the skeleton-loading pulse). Under reduced-motion these near-zero out so the interface settles with no perceptible spatial movement, while functional feedback that conveys state (spinners, progress fills, focus rings) is preserved.
