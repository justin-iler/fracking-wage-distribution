/*
Nonparametric Analysis: Geographic Dispersion of Economic Shocks
Feyrer, Mansur & Sacerdote (AER 2015)

Methods:
  1. Kernel density estimation — ΔWages for high vs low fracking counties
  2. Local polynomial regression — ΔWages on NewValue (tests linearity)

Run AFTER quantile_regression_fracking.do, or set working directory first:
  cd "C:\Users\YourName\Desktop\Nonparametric"
*/

clear all
set more off

* ── 1. Load & filter data ────────────────────────────────────────────────────

use "BLS_IRS_fossil_working.dta", clear

keep if industry == "1 Total"
keep if sample   == 1
keep if !missing(d_wages_capita, newvalue_capita)

egen county_id = group(fips)
egen year_id   = group(year)

* ── 2. Two-way demeaning (same as QR script) ─────────────────────────────────

foreach var of varlist d_wages_capita newvalue_capita {
    gen `var'_dm = `var'
}

local tol     = 1e-8
local maxiter = 100

forvalues iter = 1/`maxiter' {
    gen _y_old = d_wages_capita_dm
    gen _x_old = newvalue_capita_dm

    foreach var of varlist d_wages_capita_dm newvalue_capita_dm {
        qui bysort county_id: egen _mean = mean(`var')
        qui replace `var' = `var' - _mean
        drop _mean
    }
    foreach var of varlist d_wages_capita_dm newvalue_capita_dm {
        qui bysort year_id: egen _mean = mean(`var')
        qui replace `var' = `var' - _mean
        drop _mean
    }

    qui gen _diff_y = abs(d_wages_capita_dm - _y_old)
    qui gen _diff_x = abs(newvalue_capita_dm - _x_old)
    qui summarize _diff_y, meanonly
    local dy = r(max)
    qui summarize _diff_x, meanonly
    local dx = r(max)
    drop _y_old _x_old _diff_y _diff_x

    if max(`dy', `dx') < `tol' {
        di "Demeaning converged at iteration `iter'"
        continue, break
    }
}

* ── 3. Define high / low fracking groups ─────────────────────────────────────
* "High fracking" = above 90th percentile of NewValue (among all obs)
* "Low  fracking" = NewValue == 0 (no fossil fuel production)

qui summarize newvalue_capita, detail
local p90 = r(p90)

gen fracking_group = .
replace fracking_group = 0 if newvalue_capita == 0          // no fracking
replace fracking_group = 1 if newvalue_capita >= `p90' & ///
                               !missing(newvalue_capita)    // high fracking

label define fgrp 0 "No fracking (NewValue=0)" 1 "High fracking (top 10%)"
label values fracking_group fgrp

tab fracking_group, missing

* ── 4. Kernel Density: ΔWages by fracking group ──────────────────────────────
* Winsorize at 1st/99th pct of demeaned wages for cleaner plots

qui summarize d_wages_capita_dm, detail
local p1  = r(p1)
local p99 = r(p99)

gen d_wages_plot = d_wages_capita_dm
replace d_wages_plot = `p1'  if d_wages_plot < `p1'
replace d_wages_plot = `p99' if d_wages_plot > `p99'

* --- Figure 1: Kernel density overlay ---
twoway ///
    (kdensity d_wages_plot if fracking_group == 0, ///
        lcolor(navy) lwidth(medthick) lpattern(solid)) ///
    (kdensity d_wages_plot if fracking_group == 1, ///
        lcolor(cranberry) lwidth(medthick) lpattern(dash)), ///
    legend(order(1 "No fracking (NewValue=0)" ///
                 2 "High fracking (top 10%)") ///
           position(1) ring(0)) ///
    xtitle("ΔWages per capita (FE-demeaned, winsorized)", size(small)) ///
    ytitle("Density", size(small)) ///
    title("Distribution of Wage Changes by Fracking Intensity", size(medsmall)) ///
    subtitle("County + year FE removed via within-group demeaning", size(small)) ///
    note("Source: BLS/IRS data, Feyrer et al. (2015). High fracking = NewValue ≥ p90.", ///
         size(vsmall)) ///
    scheme(s2color) graphregion(color(white))

graph export "fig1_kdensity_fracking_groups.png", replace width(1200)
di "Saved: fig1_kdensity_fracking_groups.png"

* --- Figure 2: Cumulative distribution (CDF) overlay ---
* Useful for stochastic dominance — does high fracking shift the whole dist?

cumul d_wages_plot if fracking_group == 0, gen(cdf_no)   equal
cumul d_wages_plot if fracking_group == 1, gen(cdf_high) equal

twoway ///
    (line cdf_no   d_wages_plot if fracking_group == 0, ///
        sort lcolor(navy) lwidth(medthick)) ///
    (line cdf_high d_wages_plot if fracking_group == 1, ///
        sort lcolor(cranberry) lwidth(medthick) lpattern(dash)), ///
    legend(order(1 "No fracking" 2 "High fracking (top 10%)") ///
           position(5) ring(0)) ///
    xtitle("ΔWages per capita (FE-demeaned)", size(small)) ///
    ytitle("Cumulative probability", size(small)) ///
    title("CDF of Wage Changes by Fracking Intensity", size(medsmall)) ///
    yline(0.5, lcolor(gs10) lpattern(dot)) ///
    note("Source: BLS/IRS data, Feyrer et al. (2015).", size(vsmall)) ///
    scheme(s2color) graphregion(color(white))

graph export "fig2_cdf_fracking_groups.png", replace width(1200)
di "Saved: fig2_cdf_fracking_groups.png"

* --- Kolmogorov-Smirnov test: are the two distributions different? ---
di ""
di "════════════════════════════════════════════════════"
di "  Kolmogorov-Smirnov test (no fracking vs high)"
di "════════════════════════════════════════════════════"
ksmirnov d_wages_plot, by(fracking_group)

* ── 5. Local Polynomial Regression ───────────────────────────────────────────
* Strategy: partial out FE from wages (y), then plot residuals against the
* ORIGINAL (non-demeaned) NewValue so the x-axis is interpretable.
* The demeaned y captures within-county/year variation in wages;
* plotting it against raw NewValue shows the dose-response relationship.

* Trim NewValue at 99th pct among fracking counties to avoid sparse-tail noise
qui summarize newvalue_capita if newvalue_capita > 0, detail
local xp99 = r(p99)
local xmed = r(p50)

di "NewValue p50 (fracking): `xmed'"
di "NewValue p99 (fracking): `xp99'"

* --- Figure 3: Local linear regression with CI ---
* y = FE-demeaned wages (within variation), x = original NewValue
* Restricted to fracking counties (NewValue > 0), trimmed at p99

twoway ///
    (lpolyci d_wages_capita_dm newvalue_capita ///
        if newvalue_capita > 0 & newvalue_capita <= `xp99', ///
        degree(1) kernel(epanechnikov) ///
        lcolor(navy) lwidth(medthick) ///
        acolor(navy%20)) ///
    (lfit d_wages_capita_dm newvalue_capita ///
        if newvalue_capita > 0 & newvalue_capita <= `xp99', ///
        lcolor(cranberry) lwidth(medium) lpattern(dash)), ///
    legend(order(2 "Local linear (lpoly)" ///
                 3 "OLS linear fit") ///
           position(11) ring(0)) ///
    xtitle("NewValue per capita (original, $/capita)", size(small)) ///
    ytitle("ΔWages per capita (FE-demeaned residual)", size(small)) ///
    title("Local Polynomial Regression: Wages on New Fossil Fuel Value", ///
          size(medsmall)) ///
    subtitle("Fracking counties only (NewValue > 0), trimmed at p99", size(small)) ///
    note("Epanechnikov kernel, degree 1. Shaded area = 95% CI." ///
         "Dashed line = OLS fit. Source: Feyrer et al. (2015).", size(vsmall)) ///
    scheme(s2color) graphregion(color(white))

graph export "fig3_lpoly_linear.png", replace width(1200)
di "Saved: fig3_lpoly_linear.png"

* --- Figure 4: Compare degree 1 vs degree 2 (test for nonlinearity) ---
twoway ///
    (lpoly d_wages_capita_dm newvalue_capita ///
        if newvalue_capita > 0 & newvalue_capita <= `xp99', ///
        degree(1) kernel(epanechnikov) ///
        lcolor(navy) lwidth(medthick)) ///
    (lpoly d_wages_capita_dm newvalue_capita ///
        if newvalue_capita > 0 & newvalue_capita <= `xp99', ///
        degree(2) kernel(epanechnikov) ///
        lcolor(cranberry) lwidth(medium) lpattern(dash)) ///
    (lfit d_wages_capita_dm newvalue_capita ///
        if newvalue_capita > 0 & newvalue_capita <= `xp99', ///
        lcolor(gs8) lwidth(thin) lpattern(dot)), ///
    legend(order(1 "Local linear (degree 1)" ///
                 2 "Local quadratic (degree 2)" ///
                 3 "OLS") ///
           position(11) ring(0)) ///
    xtitle("NewValue per capita (original, $/capita)", size(small)) ///
    ytitle("ΔWages per capita (FE-demeaned residual)", size(small)) ///
    title("Linearity Check: Local Polynomial Degrees 1 vs 2", size(medsmall)) ///
    note("Epanechnikov kernel. Fracking counties only, trimmed at p99." ///
         "Source: Feyrer et al. (2015).", size(vsmall)) ///
    scheme(s2color) graphregion(color(white))

graph export "fig4_lpoly_degree_comparison.png", replace width(1200)
di "Saved: fig4_lpoly_degree_comparison.png"

* ── 6. Summary statistics by group ──────────────────────────────────────────

di ""
di "════════════════════════════════════════════════════════"
di "  Mean ΔWages by fracking group (demeaned)"
di "════════════════════════════════════════════════════════"
tabstat d_wages_capita_dm, by(fracking_group) ///
    stats(n mean sd p10 p50 p90) format(%10.4f)

di ""
di "All figures saved to working directory."
di "Files: fig1_kdensity_fracking_groups.png"
di "       fig2_cdf_fracking_groups.png"
di "       fig3_lpoly_linear.png"
di "       fig4_lpoly_degree_comparison.png"
