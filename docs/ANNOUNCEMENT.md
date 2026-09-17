# 🚀 Experimental Release & Call for Collaboration: WRF C-Z0HR PBL Regularization (`wrf-pbl-cz0hr-exp`)

We are inviting community testing and evaluation of an experimental modification to the **MYNN-EDMF** Planetary Boundary Layer (PBL) scheme in the Weather Research and Forecasting (WRF) model.

This release implements **Curvature-Aware Zero-Offset Hyperbolic Regularization (C-Z0HR)**, an upstream diagnostic scalar mapping of the gradient Richardson number profile, designed to resolve **spurious quenching, numerical stiffness, and runaway nocturnal surface cold-pool biases** in stable boundary layer (SBL) simulations.

> **Status:** Experimental research code. A single namelist switch, `bl_mynn_cz0hr` (`0` = baseline, `1` = regularized), enables A/B testing in a single build with no recompilation between runs.

---

## 1. Scientific & Operational Motivation

Under strongly stable nocturnal conditions or polar surface cooling, unregularized PBL schemes encounter two numerical pathologies:

- **LLJ noise amplification.** Near Low-Level Jet cores where vertical wind shear vanishes ($S \equiv |\partial \mathbf{V}/\partial z| \to 0$), the diagnostic Richardson gain scales as $\mathcal{O}(S^{-2})$; differentiating raw velocity gradients amplifies grid/measurement noise as $\mathcal{O}(S^{-6})$, injecting high-frequency noise into the solver.
- **Runaway cold-pool divergence.** When under-resolved grid curvature pushes the discrete $Ri_g$ across the critical threshold ($Ri_c \approx 0.20\text{–}0.25$), eddy diffusivities collapse ($K_h \to 0$), halting downward sensible heat transport and decoupling the surface energy balance, driving spurious cold pools.

Rather than ad-hoc "long-tail" over-mixing, C-Z0HR shifts regularization **upstream** into the diagnostic scalar input ($Ri_g \to Ri_g^{\mathrm{reg}}$) that feeds the MYNN stability closures.

---

## 2. Algorithmic & Practical Advantages

1. **Zero retuning of baseline physics.** C-Z0HR modifies only the diagnostic argument passed into the stability functions; the Mellor–Yamada closure algebra and empirical constants remain untouched.
2. **Grid-resolved curvature bounding.** Local profile curvature is evaluated across the vertical grid via a non-uniform central second derivative, expanding at LLJ noses to attenuate threshold sensitivity by $1/(1+\alpha)$ (a 66.7% reduction for $\alpha = 2$).
3. **Smooth ($C^\infty$) formulation.** The piecewise conditional is replaced by the hyperbolic norm $\Phi_\epsilon(x) = \sqrt{x^2 + \epsilon^2}$, giving continuous derivatives and avoiding branch divergence in vectorized kernels.
4. **Runtime A/B toggle.** `bl_mynn_cz0hr` in `&physics` switches the feature; `0` is intended to be bit-identical to baseline, enabling clean control runs without a second build.

### Motivating results (provenance: Julia SCM, not WRF)

The following figures come from the **SBLToolkit.jl Julia single-column framework** and are the scientific motivation for this WRF port. They are *not* WRF results; obtaining them in WRF is the goal of this call.

- GABLS3 24 h: 2 m cold bias $\approx -3.2$ K (unregularized MYJ-class baseline) → $\approx -0.1$ K; LLJ height bias $\approx -45$ m → $\approx +2$ m; false-bifurcation trigger rate $\approx 38\%$ → $0\%$.

---

## 3. How to Test

**Target:** WRF ≥ v4.6 with the `phys/MYNN-EDMF` submodule (patches generated against MYNN-EDMF `90f36c25`, the commit pinned by WRF v4.7.1).

```bash
cd $WRF_DIR
git submodule update --init phys/MYNN-EDMF

# Scheme changes (inside the submodule)
git -C phys/MYNN-EDMF apply <repo>/patches/mynn_edmf_cz0hr.patch

# WRF namelist plumbing (WRF root)
git apply <repo>/patches/wrf_core_cz0hr.patch

./clean -a && ./configure && ./compile em_real   # or your case
```

In `namelist.input` (`&physics`):

```
bl_pbl_physics = 5,    ! MYNN-EDMF
bl_mynn_cz0hr  = 0,    ! 0 = baseline (control), 1 = C-Z0HR
```

**Suggested evaluations:** idealized GABLS3 diurnal SCM; CASES-99 nocturnal LLJ; a polar/cryospheric case (e.g., SHEBA). Compare baseline vs. regularized on $T_{2m}$/$T_s$ cooling, $z_{\mathrm{LLJ}}$, $H_0$, $\tau(z)$, and solver behavior.

---

## 4. Call for Collaboration

We are seeking research partners to:

- Run A/B comparisons across MYNN-EDMF configurations and other schemes.
- Test 3D real-case regional forecasts and polar simulations.
- Assess Tangent Linear / Adjoint behavior in 4D-Var workflows.

Please open an issue or discussion in the `wrf-pbl-cz0hr-exp` repository with your configuration and metrics.
