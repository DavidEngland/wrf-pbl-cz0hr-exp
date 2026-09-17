# WRF C-Z0HR PBL Regularization — Experimental Release (`wrf-pbl-cz0hr-exp`)

**Status:** Experimental research code. Not for operational use. Baseline behavior is preserved exactly when the feature is switched off.

This repository packages **Curvature-Aware Zero-Offset Hyperbolic Regularization (C-Z0HR)** for the **MYNN-EDMF** planetary boundary layer (PBL) scheme in the Weather Research and Forecasting (WRF) model. C-Z0HR is an *upstream* diagnostic regularization of the gradient Richardson number profile, designed to mitigate spurious quenching, numerical stiffness, and runaway nocturnal surface cold-pool biases in stable boundary layer (SBL) simulations.

---

## What this is (and is not)

- **Upstream scalar regularization.** C-Z0HR maps the raw gradient Richardson number profile $Ri_g \rightarrow Ri_g^{\mathrm{reg}}$ using grid cell thickness and profile curvature, *before* the MYNN stability closures ($S_m$, $S_h$) consume it. The Mellor–Yamada closure algebra itself is **not modified**.
- **A namelist toggle, not a fork.** A single switch, `bl_mynn_cz0hr` (`0` = off / baseline, `1` = on), enables A/B evaluation in one build. Setting it to `0` is intended to reproduce baseline results bit-for-bit.
- **Targets WRF ≥ v4.6.** The MYNN scheme lives in the `phys/MYNN-EDMF` submodule in modern WRF. These patches were generated against MYNN-EDMF commit `90f36c25259ec1960b24325f5b29ac7c5adeac73` (the commit pinned by WRF v4.7.1).

## Repository layout

```
patches/
  mynn_edmf_cz0hr.patch   # MYNN-EDMF scheme changes (apply inside phys/MYNN-EDMF)
  wrf_core_cz0hr.patch    # WRF plumbing: Registry + namelist + driver chain (apply at WRF root)
docs/
```

## Quick start

```bash
cd $WRF_DIR                                  # your WRF >= v4.6 checkout
git submodule update --init phys/MYNN-EDMF   # populate the submodule

# 1. Scheme changes (inside the submodule)
git -C phys/MYNN-EDMF apply /path/to/wrf-pbl-cz0hr-exp/patches/mynn_edmf_cz0hr.patch

# 2. WRF namelist plumbing (WRF root)
git apply /path/to/wrf-pbl-cz0hr-exp/patches/wrf_core_cz0hr.patch

# 3. Rebuild
./clean -a && ./configure   # then ./compile em_real (or your case)
```

### Running an A/B test

Add to `&physics` in `namelist.input`:

```
bl_pbl_physics = 5,          ! MYNN-EDMF
bl_mynn_cz0hr  = 0,          ! 0 = baseline, 1 = C-Z0HR regularized Ri_g
```

Run once with `bl_mynn_cz0hr = 0` (baseline) and once with `bl_mynn_cz0hr = 1` (regularized) and compare. Because the toggle is read at runtime, **no recompilation is needed between the two runs**, and the `0` case is the control for bit-identical verification.

## What to report back

We are seeking collaborators to run comparative A/B evaluations. The most useful diagnostics:

- Near-surface skin temperature ($T_s$) and 2 m temperature ($T_{2m}$) cooling rates through the nocturnal transition.
- Low-Level Jet nose height ($z_{\mathrm{LLJ}}$) and peak wind speed evolution.
- Vertical sensible heat flux ($H_0$) and momentum stress profiles $\tau(z)$.
- Solver/iteration behavior in stable regimes; any NaNs or runtime failures.
- Confirmation that `bl_mynn_cz0hr = 0` reproduces your baseline.

Suggested test configurations: idealized GABLS3 diurnal cycle (SCM), CASES-99 flat-terrain nocturnal LLJ, and a polar/cryospheric case (e.g., SHEBA).

## Algorithm summary

At each wall level $k$ the raw gradient Richardson number is $Ri_g = N^2/S^2 = -\mathrm{gh}/\max(\mathrm{gm}, \epsilon)$. C-Z0HR regularizes the *profile* (not the per-level value independently):

1. Non-uniform central second derivative $Ri_{zz} = \partial^2 Ri_g / \partial z^2$.
2. Grid-resolved curvature metric $C(z, \Delta z) = \tfrac{1}{2}|Ri_{zz}|(\Delta z)^2 + \epsilon_c Ri_c$.
3. $C^\infty$ hyperbolic distance $\Phi_\epsilon(x) = \sqrt{x^2 + \epsilon^2}$, with $x = Ri_g - Ri_c$.
4. Mapping $Ri_g^{\mathrm{reg}} = Ri_c + x\,\dfrac{\Phi_\epsilon + C}{\Phi_\epsilon + (1+\alpha) C}$.

With $\alpha = 2$ the slope at the critical threshold is attenuated by $1/(1+\alpha) = 1/3$. Hyperparameters (`cz0hr_alpha`, `cz0hr_eps_c`, `cz0hr_eps`) are compile-time parameters in `module_bl_mynnedmf` (kind-phys precision); `Ri_c = 0.20`.

### Where it hooks in

- `mym_level2`: the raw $Ri_g$ profile is built in a first pass, regularized by `cz0hr_regularize_ri` when `bl_mynn_cz0hr = 1`, and the regularized profile `rig` drives the flux-Richardson / $S_m$ / $S_h$ algebra in a second pass.
- `mym_turbulence`: consumes the same `rig` profile for the Canuto–Kitamura `a2fac` and the Kondo Prandtl limiter, keeping both Ri consumers consistent.

### Note on the Canuto–Kitamura interaction

MYNN-EDMF already attenuates closures at high $Ri$ via `CKmod`/`a2fac`. C-Z0HR is **complementary upstream profile shaping** — it reduces threshold sensitivity and curvature-driven noise amplification (the $\mathcal{O}(S^{-6})$ problem near zero-shear LLJ cores) before that attenuation is applied. Testers should evaluate the combination, not assume independence.

## Provenance of the benchmark numbers

The performance figures quoted in the announcement (e.g., GABLS3 2 m cold bias reduced from $\approx -3.2$ K to $\approx -0.1$ K, LLJ height bias from $\approx -45$ m to $\approx +2$ m, false-bifurcation rate from $\approx 38\%$ to $0\%$) were produced by the **SBLToolkit.jl Julia single-column framework**, not by WRF+C-Z0HR. They are the scientific motivation for this WRF port. Obtaining equivalent WRF results is the goal of this collaboration call; please treat the Julia SCM numbers as *expected direction/magnitude*, not as WRF validation.

## License & citation

- The C-Z0HR patches: see [LICENSE](LICENSE).
- MYNN-EDMF is developed by NOAA/GSL, NCAR, and collaborators (see the [MYNN-EDMF repository](https://github.com/NCAR/MYNN-EDMF)); WRF is UCAR/NCAR. Their respective licenses apply to the upstream code these patches modify.
- If you use this in a publication, please cite this repository (see [CITATION.cff](CITATION.cff)) and the MYNN-EDMF reference (Olson et al., NOAA Tech. Memo.).
