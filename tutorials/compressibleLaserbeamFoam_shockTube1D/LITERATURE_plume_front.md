# Plume and shock fronts in ns-laser ablation at 1 atm: literature and blast-wave estimates for our case

Our case: 10 mJ, 1064 nm, 5 ns, laser radius 50 um (spot 100 um; mean fluence 127 J/cm2, mean irradiance
2.5e10 W/cm2), 316L, ambient "air" phase with argon properties at 1 bar, 300 K (rho0 = 1.60 kg/m3, gamma = 5/3).

## Blast-wave scaling used in the literature
Sedov-Taylor point explosion, for a geometry exponent alpha = 3 (sphere), 2 (cylinder) or 1 (plane):

    R(t) = xi(gamma) (E/rho0)^(1/(2+alpha)) t^(2/(2+alpha)),     v = dR/dt = (2/(2+alpha)) R/t

Notation varies: some papers write (E/(K rho0))^(1/5), with K = 0.6-0.8 for laser shocks in gases (Gatti 1988).
- **Sphere:** R ~ (E/rho0)^(1/5) t^(2/5). xi = 1.033 for gamma 1.4 (air) and 1.152 for gamma 5/3 (argon).
  - A hemispherical blast on a surface carries E in half the volume: use 2E in the sphere formula.
- **Cylinder:** R ~ (E_l/rho0)^(1/4) t^(1/2), with E_l per unit length. Use 2 E_l for a half-cylinder on a
  surface (our planar 2D case).
- **Fitted energies:** E fitted to measured trajectories is usually 30-70% of the pulse energy (the rest is
  reflected, conducted, spent on ablation and radiated).
- **Valid window.** It runs from when the swept-up gas mass exceeds the ablated mass,
  R >> (3m/(2 pi rho0))^(1/3), to when the shock weakens (Mach ~2, R << (E/p0)^(1/3)). Before that window the
  front moves faster than t^(2/5): measured exponents are 0.55-0.65 in the first ~100 ns (Palanco 2014).
- **During the pulse,** above ~1e9 W/cm2, a laser-supported detonation (LSD) wave can run up the beam. Its
  speed is D ~ [2(gamma^2-1) I/rho0]^(1/3).

## Published front positions and speeds (ns pulses, 1 atm)
| # | Reference | Target | lambda, tau | Energy, spot | Ambient | Measured | Front position / speed, scaling |
|---|---|---|---|---|---|---|---|
| 1 | Kwapis, Hewitt, Hartig, Opt. Express 31, 10694 (2023), [doi:10.1364/OE.483055](https://doi.org/10.1364/OE.483055) (open access) | Ce metal | 1064 nm, 7 ns | 20.8 / 60 / 100 mJ, ~250 um spot (40-200 J/cm2) | air and **argon**, 760 and 100 Torr | focused shadowgraphy, 30 ns - 10 us | Hemispherical shock. At 1 atm no linear free-expansion phase is seen after 30-50 ns. Sedov-Taylor fits 60-250 ns within <2%; Mach > 10 early. Air fronts are farther and faster than argon (lower density). Power law with exponent > 2/5 later; ionisation front behind the shock 0.4-0.9 us. (Numbers in Figs. 4-5.) |
| 2 | Palanco et al., Opt. Express 22, 16552 (2014), [doi:10.1364/OE.22.016552](https://doi.org/10.1364/OE.22.016552) (open access) | Cu | 1064 nm, 6 ns | (a) 500 mJ / 0.8 mm (100 J/cm2), (b) 50 mJ / 0.8 mm (10 J/cm2), (c) 500 mJ / 2.5 mm | air, 1 atm | Schlieren, 0-20 us; emission | Initial shock speed 5.0e4 m/s (a) and 2.9e4 m/s (b). Early power-law exponents 0.62 / 0.65 / 0.55: faster than t^(2/5) at first, then decelerates. Sedov-Taylor energies at late times 331 / 33 / 318 mJ (64-67% of the pulse). Plume (emission) stopping distance ~5 mm (a), ~4 mm (b). |
| 3 | Nagel et al., Opt. Express 29, 33481 (2021), [doi:10.1364/OE.440212](https://doi.org/10.1364/OE.440212) (open access) | Al | 1064 nm, 10 ns | 40 mJ, 1.5 mm spot | air, 1 atm | shadowgraphy, 0.1-4 us | Sedov fit with free dimensionality: alpha = 2.75 +- 0.11 (nearly hemispherical); E_d = 28.5 +- 18.2 mJ (~70% of 40 mJ). Assuming alpha = 3 underestimates E_d. No data before 100 ns. |
| 4 | Harilal, Miloshevsky, Diwakar, LaHaye, Hassanein, Phys. Plasmas 19, 083504 (2012), [OSTI record](https://www.osti.gov/biblio/22086083); [author PDF](https://engineering.purdue.edu/CMUXE/Publications/Harilal/2012POP-Hari.pdf) | Al | 1064 nm, 6 ns | (not checked; the same group's fs comparison quotes ~100 mJ for its ns pulse) | **argon, 1 atm** | focused shadowgraphy + 2 ns gated imaging, with CFD | Plume front from fast photography agrees with the shock front from shadowgraphy up to <= 400 ns, then the plume is confined while the shock continues. Secondary shock behind the primary at 100-500 ns. CFD matches the early expansion. (Positions in figures, not checked here: the PDF could not be fetched from this machine.) |
| 5 | Hussein, Diwakar, Harilal, Hassanein, J. Appl. Phys. 113, 143305 (2013), [author PDF](https://engineering.purdue.edu/CMUXE/Publications/Harilal/2013-JAP-Amina.pdf) | Al | 1064 / 532 / 266 nm, ns | (not checked) | air, 1 atm | shadowgraphy + fast photography, 0-500 ns | Wavelength dependence of plume generation, expansion and confinement in the first 500 ns at 1 atm (abstract only; values and conclusions not checked). |
| 6 | Corsi et al., Appl. Spectrosc. 57, 715 (2003), [doi:10.1366/000370203322005436](https://doi.org/10.1366/000370203322005436) (paywalled) | **steel** | Nd:YAG | two energy regimes | air, 1 atm | shadowgraphy + optical time-of-flight | The only steel shadowgraphy study found. The plume has an air/continuum region driven by shock ionisation ahead of the steel-vapour region. Front speeds not verified (abstract only). |
| 7 | Hermann, Axente, Craciun, Taleb, Pelascini, Spectrochim. Acta B 143, 63 (2018) | **steel** | UV, ns | - | **argon**, ~1 atm | spectroscopy, 200 ns - 12 us | No front positions. The plume pressure drops below the ambient within microseconds; vapour-gas interdiffusion is negligible on that time scale. This supports a sharp vapour/argon contact, as in our VoF model. |
| - | Campanella et al., Atoms 7, 57 (2019), [doi:10.3390/atoms7020057](https://doi.org/10.3390/atoms7020057) (review) | - | - | - | - | - | Sedov-Taylor for laser shocks. K(gamma) = 0.6-0.8 fits Nd:YAG shocks in air/N2/He (Gatti et al., Opt. Commun. 69, 141 (1988)). Holds down to Mach ~2. |

I found no open-access paper giving numerical shock or plume positions for steel or iron in 1 atm argon at
10-100 ns. The closest in energy, spot size and ambient gas is #1 (Ce, 20.8 mJ, ~250 um spot, 1 atm argon,
1064 nm, 7 ns). Next is #2(b) for the early speed at a similar irradiance (Cu, 50 mJ, 2.9e4 m/s initial shock
speed), then #3 for the energy fraction and dimensionality. For steel specifically, #6 has the shadowgrams but
is paywalled. Iron is heavier than Al and Cu, so the vapour front lags the shock more and decelerates
earlier. The shock itself depends only on E_d and rho0.

## What this means for our case (10 mJ, 316L, argon 1 bar)
Hemispherical Sedov-Taylor, R = 1.152 (2 f E/rho0)^(1/5) t^(2/5), E = 10 mJ, rho0 = 1.60 kg/m3, f = fraction of
the pulse that ends up in the blast. Times are from the energy release (about the pulse peak, ~2.5 ns after
the start of a 5 ns pulse).

| t | f = 0.3 | f = 0.5 | f = 0.7 | f = 1.0 |
|---|---|---|---|---|
| 5 ns | 180 um, 14 km/s | 200 um, 16 km/s | 213 um, 17 km/s | 229 um, 18 km/s |
| 10 ns | 238 um, 9.5 km/s | 263 um, 10.5 km/s | 282 um, 11 km/s | 303 um, 12 km/s |
| 20 ns | 314 um, 6.3 km/s | 348 um, 7.0 km/s | 372 um, 7.4 km/s | 399 um, 8.0 km/s |
| 50 ns | 453 um, 3.6 km/s | 501 um, 4.0 km/s | 536 um, 4.3 km/s | 576 um, 4.6 km/s |
| 100 ns | 597 um, 2.4 km/s | 662 um, 2.6 km/s | 708 um, 2.8 km/s | 760 um, 3.0 km/s |

The literature fractions are 0.33-0.7 (#2, #3), so the expected **shock front is about 250-280 um at 10 ns,
330-370 um at 20 ns and 450-540 um at 50 ns, moving at 10, 7 and 4 km/s**. This assumes the Sedov-Taylor
regime already holds, which is optimistic at these times:
- **Ablated mass:** about 0.1-0.3 ug for a 100 um spot at ~130 J/cm2 (a few um deep).
- **Mass condition:** that mass needs R >> 310-450 um before the swept-up argon dominates. So at 10-50 ns
  the flow is still in the piston (vapour-driven) phase.
- **Shape of the early trajectory:** measured fronts in this phase run faster than t^(2/5) (exponents
  0.55-0.65, #2). At 1 atm they lose the linear free-expansion stage within 30-50 ns (#1).
- **During the pulse:** at 2.5e10 W/cm2 an LSD wave is possible, D ~ 80 km/s for full absorption. Measured
  early speeds at comparable irradiance are 3-5e4 m/s (#2).
- **Plume (vapour/contact) front:** it moves with the shock until detachment, which happens at < 50 ns at
  1 atm (#1, #4). After that it lags behind at roughly 0.5-0.8 R_shock and stalls within ~1 us. The
  stopping distance is mm-scale for 50 mJ (#2), so less for 10 mJ.

**So for 10-50 ns the simulation should show a hemispherical shock at roughly 250-550 um and the vapour front a
little behind it (about 150-400 um), both decelerating from >10 km/s to ~4 km/s.** The pressure-equilibrium
radius (E/p0)^(1/3) is about 3.7 mm (f = 0.5), so the shock stays strong (Mach > 5 in argon, c0 = 322 m/s)
well past 100 ns. The 3D domain (1.5 mm in each direction) should be reached by the shock at about 0.8 us
(f = 0.5).

Use `tutorials/tools/plumeFront.py <case> --E 0.005 --rho0 1.60` (E = f x 10 mJ) to draw this line on the
simulated plume and shock heights.

For the planar 2D plume-edge test case (tutorials/compressibleLaserbeamFoam_plumeEdge2D), the initial vapour
half-disc holds E_l = p V/(gamma-1) ~ 0.21 J/m. The matching half-cylinder blast is
R ~ (2 E_l/rho0)^(1/4) t^(1/2): about 51 um at 5 ns, 72 um at 10 ns and 101 um at 20 ns (xi ~ 1, an
estimate only: the initial disc of radius 30 um is not a point source, so this holds only once R >> 30 um).
Use `--El 0.21`.
