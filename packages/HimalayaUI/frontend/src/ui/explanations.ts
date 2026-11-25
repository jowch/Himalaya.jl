export const IMAGE_INFO_TERMS: Readonly<Record<string, string>> = {
  Kind: "One image as measured, or the mean of the selected images, each usually already divided by its own I₂. Averaging N images cuts the noise by about √N.",
  Exposure: "How long the detector counted. Scattered counts grow in proportion, so counting longer reveals weaker peaks, at the cost of more radiation damage to the sample.",
  "Flight tube": "The evacuated tube between sample and detector. A longer tube puts the detector farther away, so it reaches lower q at finer resolution but covers a narrower q range.",
  Energy: "The X-ray photon energy. It fixes the wavelength, λ (Å) = 12.398 / E (keV), which turns each scattering angle 2θ into q = 4π sin θ / λ.",
  "I₀": "The ion chamber upstream of the sample, reading the incident beam. I₂ / I₀ is proportional to the sample's transmission.",
  "I₂": "The photodiode in the beamstop, reading the beam that got through the sample: incident flux times transmission. Traces are usually divided by it, which corrects for beam intensity and absorption but not for how much material the beam crossed.",
};

export const SCORE =
  "How well this lattice explains the peaks, from 0 to 1. Low-order reflections are the strongest and appear first, so each found peak counts by 1 / its order. Peaks from one lattice are similarly sharp, so mixed sharpness lowers the score.";

export const CURVATURE =
  "The average Gaussian curvature of the minimal surface the bilayer midplane follows, from Gauss–Bonnet: 2πχ / (A₀a²), with χ the Euler characteristic and A₀a² the surface area per unit cell of lattice parameter a. Cubic phases related by the Bonnet transformation share it, so Pn3m, Im3m and Ia3d that coexist in equilibrium should show the same κ.";

export const R_SQUARED =
  "How closely the peaks follow this lattice. Each peak's q against its reflection's ratio should fall on a straight line through the origin, whose slope sets the lattice parameter. 1 is a perfect fit.";

export const NOT_FOUND =
  "No lattice of this phase fits enough of the peaks in the selected images to be offered. The call belongs to the sample, so it stays, with lines drawn from the lattice it was called on.";

export const RAW = "The measured trace, with its background's scattering still in it.";

export const minusExplanation = (background: string): string =>
  `The current trace minus ${background}, scaled as its chip shows. What remains is the sample's excess scattering over its background.`;
