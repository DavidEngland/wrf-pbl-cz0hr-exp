**Algorithm Summary**

At each interior wall level $k$, the raw gradient Richardson number $Ri_g$ is extracted before stability closure evaluations. C-Z0HR regularizes the full profile across vertical columns:

1. **Non-Uniform Grid Hessian:** Evaluates pointwise second derivatives with exact $\mathcal{O}(h^2)$ accuracy on stretched coordinate layers:

$$Ri_{zz} = \frac{2 \left[ \Delta z_{\text{dn}} Ri_{k+1} - (\Delta z_{\text{dn}} + \Delta z_{\text{up}}) Ri_k + \Delta z_{\text{up}} Ri_{k-1} \right]}{\Delta z_{\text{dn}} \Delta z_{\text{up}} (\Delta z_{\text{dn}} + \Delta z_{\text{up}})}$$


2. **Upstream Curvature Metric:** Scales local curvature using the upstream cell thickness $\Delta z_{\text{dn}}$ supplying interface wall $k$, with a background noise floor $C_{\min} = \epsilon_c Ri_c$:

$$C(z, \Delta z) = \frac{1}{2} \vert{}Ri_{zz}\vert{} (\Delta z_{\text{dn}})^2 + \epsilon_c Ri_c$$


3. **$C^\infty$ Hyperbolic Distance:** Measures distance $x = Ri_g - Ri_c$ smoothly relative to critical threshold $Ri_c = 0.20$:

$$\Phi_\epsilon(x) = \sqrt{x^2 + \epsilon^2}$$


4. **Upstream Mapping:** Regularizes the raw input profile $Ri_g \to Ri_g^{\text{reg}}$ in-place:

$$Ri_g^{\text{reg}} = Ri_c + x \cdot \frac{\Phi_\epsilon + C}{\Phi_\epsilon + (1 + \alpha) C}$$



Setting $\alpha = 2.0$ yields asymptotic threshold slope attenuation of $\frac{1}{1 + \alpha} = \frac{1}{3}$ ($66.7\%$ reduction in threshold sensitivity at $Ri_c$). Hyperparameters (`cz0hr_alpha = 2.0`, `cz0hr_eps_c = 0.05`, `cz0hr_eps = 1.0e-4`, `cz0hr_ric = 0.20`) are declared as `kind_phys` module parameters in `module_bl_mynnedmf.F90`.

**Code Hook Locations**

* **`mym_level2` Interception:** Called immediately after `bl_mynn_ess` populates the raw `ri` array when `bl_mynn_cz0hr == 1`. The regularized profile directly feeds downstream flux Richardson number ($Rf$) calculations and stability closures ($S_m, S_h$).
* **`mym_turbulence` Consistency:** Passes the regularized `ri` profile through to Canuto–Kitamura `a2fac` scaling and the Kondo Prandtl number limiter, ensuring identical profile inputs across all diagnostic consumers.
* **Namelist Integration:** Controlled via `bl_mynn_cz0hr` in `namelist.input` (wired into `Registry/Registry.EM_COMMON`), allowing runtime switching without re-compilation and bit-identical baseline verification when set to `0`.

**Canuto–Kitamura (`CKmod`) Interaction**

MYNN-EDMF applies downstream closure damping at high $Ri$ via `CKmod`/`a2fac`. C-Z0HR serves as an **upstream pre-conditioning operator**—it suppresses noise amplification ($\mathcal{O}(S^{-6})$ near zero-shear Low-Level Jet noses) and smoothes threshold transitions before `CKmod` evaluates. The two mechanisms act complementarily; community testers should evaluate their combined impact during nocturnal SBL benchmarks.