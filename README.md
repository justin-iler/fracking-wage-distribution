# Beyond the Mean: Nonparametric and Quantile Extensions of Feyrer, Mansur & Sacerdote (2017)

**David Mears, Justin Iler, and James Canning** | Spring 2026

This repository contains the replication and extension code for our empirical economics project. We extend [Feyrer, Mansur & Sacerdote (2017)](https://doi.org/10.1257/aer.20151326) — hereafter FMS — which estimates the average causal effect of new oil and gas production on local US wages using IV regression. Our extensions ask three questions the original paper cannot address:

1. **Is the wage–fracking relationship linear?** Local polynomial regression reveals a hump-shaped treatment response.
2. **Is the causal effect homogeneous across the wage-growth distribution?** IV quantile regression shows the contemporaneous effect is slightly larger at Q90 (\~$94,000) than Q10 (\~$86,000), but the lagged give-back is near zero at Q10 and \~−$37,000 at Q90 — so the net effect is largest for low-growth counties: a boom–bust asymmetry the mean conceals.
3. **Does fracking widen spatial wage inequality?** Unconditional RIF regression shows the 90th-percentile effect ($31,000 per million dollars of production) is over six times the median effect ($5,000).

---

## Repository Structure

```
.
├── code/
│   ├── quantile_regression_fracking.do          # QR at τ = {0.10, 0.50, 0.90}, with lag — Table 1 / Table A2
│   ├── quantile_regression_lagged_byindustry.do # QR by industry with lag — Table A2
│   ├── nonparametric_fracking.do                # KDE, CDF, local polynomial — Figures 1–3
│   ├── nonparametric_fracking_lagged_byindustry.do  # Nonparametric extension by industry — Figure 3 panels
│   └── figure3_nonparametric.do                 # Bandwidth robustness checks — Figure 3
├── results/
│   ├── quantile_regression_results.csv          # QR coefficient table (no lag)
│   ├── quantile_lagged_results_byindustry.csv   # QR by industry with lag
│   └── quantile_lagged_results_total.csv        # QR total with lag
├── figures/
│   ├── fig1_kdensity_fracking_groups.png
│   ├── fig2_cdf_fracking_groups.png
│   ├── fig3_lpoly_linear.png
│   ├── fig4_lpoly_degree_comparison.png
│   ├── fig_np_1_distributional.png
│   ├── fig_np_2_linearity.png
│   ├── fig_np_3a_lpoly_total.png
│   ├── fig_np_3b_lpoly_mining.png
│   ├── fig_np_3c_lpoly_transport.png
│   ├── fig_np_3d_lpoly_construction.png
│   ├── fig_np_3e_lpoly_manufacturing.png
│   ├── fig_np_3f_lpoly_educ.png
│   ├── fig_np_3g_lpoly_services.png
│   ├── fig_np_3h_lpoly_government.png
│   ├── fig_np_4a_kde_mining.png
│   ├── fig_np_4b_kde_transport.png
│   ├── fig_np_4c_kde_manufacturing.png
│   └── fig_np_4d_kde_government.png
├── paper/
│   └── Mears_Iler_Canning_Fracking.pdf          # Final report
├── slides/
│   └── Mears_Iler_Canning_Slides.pdf            # Presentation slides
└── README.md
```

---

## Data

**We do not redistribute the underlying data.** The working datasets (`BLS_IRS_fossil_working.dta`, `BLS_IRS_fossil_cz_working.dta`, `BLS_IRS_fossil_state_working.dta`, `BLS_Distance_working.dta`) are part of the official FMS replication package, archived by the American Economic Association:

> Feyrer, James, Mansur, Erin T., and Sacerdote, Bruce. *Replication data for: Geographic Dispersion of Economic Shocks: Evidence from the Fracking Revolution.* AEA / ICPSR, 2019. https://doi.org/10.3886/E113098V1

To replicate our results:
1. Download the replication package from the link above.
2. Place the `.dta` files in a local `data/` directory (not tracked by this repo — see `.gitignore`).
3. Update the `cd` / `use` paths at the top of each `.do` file to point to your local `data/` directory.

---

## Replication

All code is written in **Stata**. The do-files are self-contained except for the data dependency above. Recommended run order:

```stata
// 1. Nonparametric analysis (Figures 1–3)
do code/nonparametric_fracking.do

// 2. Bandwidth robustness (Figure 3 panels)
do code/figure3_nonparametric.do

// 3. Nonparametric by industry
do code/nonparametric_fracking_lagged_byindustry.do

// 4. Quantile regression — total industry (Table 1)
do code/quantile_regression_fracking.do

// 5. Quantile regression — by industry with lag (Table A2)
do code/quantile_regression_lagged_byindustry.do
```

Key packages required: `bsqreg`, `ivqreg2` (or `ivqte`), `rifreg`. Install via `ssc install <package>` if not already present.

---

## Methods Summary

All specifications absorb county and year fixed effects via iterative Gauss–Seidel within-group demeaning (Canay, 2011) before estimation.

| Extension | Method | Figures / Tables |
|---|---|---|
| Distributional comparison | Kernel density + KS test | Figures 1–2 |
| Linearity test | Local polynomial regression (Epanechnikov, degree 1 & 2) | Figures 3–4 |
| Conditional quantile effects | Koenker–Bassett QR + Chernozhukov–Hansen IVQR | Table 1, Table A2 |
| Unconditional distributional effects | Firpo–Fortin–Lemieux RIF regression | Table 2 |

---

## Citation

If you use or build on this code, please cite:

> Mears, D., Iler, J., and Canning, J. (2026). "Beyond the Mean: Nonparametric and Quantile Extensions of Feyrer, Mansur & Sacerdote (2017)." Unpublished manuscript, Spring 2026.

And the original paper:

> Feyrer, J., Mansur, E. T., and Sacerdote, B. (2017). "Geographic Dispersion of Economic Shocks: Evidence from the Fracking Revolution." *American Economic Review*, 107(4), 1313–1334. https://doi.org/10.1257/aer.20151326

---

## License

Code in this repository is released under the [MIT License](https://opensource.org/licenses/MIT). The underlying data are subject to the AEA data license; see the [openICPSR deposit](https://doi.org/10.3886/E113098V1) for terms.
