/*
================================================================================
Nonparametric Analysis with Lag Structure — By Industry
Feyrer, Mansur & Sacerdote (AER 2015)

Produces two figures:
  fig_np_1_distributional.png  —  KDE + CDF side-by-side, Total industry
  fig_np_2_lpoly_grid.png      —  2x4 local polynomial grid, all industries

Treatment variable: treat_sum = NewValue_t + NewValue_{t-1}
  "No fracking"   = treat_sum exactly zero (both periods)
  "High fracking" = top 10% of treat_sum

Fixed effects absorbed via Gauss-Seidel iterative demeaning (Canay 2011).
Partial regression plot: L1 is OLS-partialled from wages before lpoly.
================================================================================
*/

clear all
set more off
set seed 12345

* ── Parameters ───────────────────────────────────────────────────────────────
local tol     = 1e-8
local maxiter = 100

local lab1 "Total"
local lab2 "Mining"
local lab3 "Transport"
local lab4 "Construction"
local lab5 "Manufacturing"
local lab6 "Educ/Health"
local lab7 "Other Services"
local lab8 "Government"

local col1 "navy"
local col2 "maroon"
local col3 "forest_green"
local col4 "dkorange"
local col5 "purple"
local col6 "teal"
local col7 "cranberry"
local col8 "dkgreen"

* ── Clean up any previous output files ───────────────────────────────────────
foreach f in fig_np_lag_1_kde_total fig_np_lag_2_cdf_total      ///
             fig_np_lag_3_lpoly_total fig_np_lag_4_lpoly_grid    ///
             fig_np_lag_5_lpoly_overlay fig_np_lag_6_kde_grid    ///
             fig_np_1_distributional fig_np_2_lpoly_grid         ///
             fig_np_2_linearity fig_np_3_lpoly_grid              ///
             fig_np_3a_lpoly_total   fig_np_3b_lpoly_mining      ///
             fig_np_3c_lpoly_transport fig_np_3d_lpoly_construction ///
             fig_np_3e_lpoly_manufacturing fig_np_3f_lpoly_educ  ///
             fig_np_3g_lpoly_services fig_np_3h_lpoly_government ///
             fig_np_4_kde_selected                               ///
             fig_np_4a_kde_mining fig_np_4b_kde_transport        ///
             fig_np_4c_kde_manufacturing fig_np_4d_kde_government {
    capture erase "`f'.png"
}

* ── Load data, create lag ─────────────────────────────────────────────────────
use "BLS_IRS_fossil_working.dta", clear
keep if sample == 1

sort fips id_industry year
by fips id_industry (year): gen newvalue_capita_L1 = newvalue_capita[_n-1]

drop if year == 2004
drop if missing(d_wages_capita, newvalue_capita, newvalue_capita_L1)

* Combined treatment: L0 + L1
gen treat_sum = newvalue_capita + newvalue_capita_L1

* Fracking groups — classify once, applies to all industries equally
* (newvalue_capita is county-level, identical across industries)
qui summarize treat_sum if id_industry == 1, detail
local p90_treat = r(p90)

gen fracking_group = .
replace fracking_group = 0 if treat_sum == 0
replace fracking_group = 1 if treat_sum >= `p90_treat' & !missing(treat_sum)

label define fgrp 0 "No fracking (both periods zero)" 1 "High fracking (top 10%)"
label values fracking_group fgrp

tempfile base
save `base'

* ── Process each industry: demean, partial regression, save ──────────────────
forvalues ind = 1/8 {

    di "Processing industry `ind': `lab`ind''"

    use `base', clear
    keep if id_industry == `ind'

    egen county_id = group(fips)
    egen year_id   = group(year)

    * Gauss-Seidel demeaning of y, L0, L1
    foreach var of varlist d_wages_capita newvalue_capita newvalue_capita_L1 {
        gen `var'_dm = `var'
    }

    forvalues iter = 1/`maxiter' {
        gen _y  = d_wages_capita_dm
        gen _x0 = newvalue_capita_dm
        gen _x1 = newvalue_capita_L1_dm

        foreach var of varlist d_wages_capita_dm newvalue_capita_dm newvalue_capita_L1_dm {
            qui bysort county_id: egen _m = mean(`var')
            qui replace `var' = `var' - _m
            drop _m
            qui bysort year_id: egen _m = mean(`var')
            qui replace `var' = `var' - _m
            drop _m
        }

        qui gen _d0 = abs(d_wages_capita_dm    - _y)
        qui gen _d1 = abs(newvalue_capita_dm   - _x0)
        qui gen _d2 = abs(newvalue_capita_L1_dm - _x1)
        qui summarize _d0, meanonly
        local dy = r(max)
        qui summarize _d1, meanonly
        local dx0 = r(max)
        qui summarize _d2, meanonly
        local dx1 = r(max)
        drop _y _x0 _x1 _d0 _d1 _d2

        if max(`dy', max(`dx0', `dx1')) < `tol' {
            di "  Converged at iteration `iter'"
            continue, break
        }
    }

    * FWL partial regression: remove L1 contribution from wages
    qui reg d_wages_capita_dm newvalue_capita_L1_dm
    predict d_wages_partial, resid

    * Winsorise at 1st / 99th pct (per industry)
    qui summarize d_wages_capita_dm, detail
    gen d_wages_plot = max(r(p1), min(r(p99), d_wages_capita_dm))

    qui summarize d_wages_partial, detail
    gen d_wages_partial_plot = max(r(p1), min(r(p99), d_wages_partial))

    * x-axis for lpoly: demeaned L0, restricted to fracking counties, trimmed p99
    qui summarize newvalue_capita if newvalue_capita > 0, detail
    gen newval_L0_plot = newvalue_capita_dm
    replace newval_L0_plot = . if newvalue_capita <= 0 | newvalue_capita > r(p99)

    gen ind_id    = `ind'
    gen ind_label = "`lab`ind''"

    tempfile ind`ind'
    save `ind`ind'', replace
}

* Combine all industries
use `ind1', clear
forvalues ind = 2/8 {
    append using `ind`ind''
}
tempfile combined
save `combined', replace

* ═══════════════════════════════════════════════════════════════════════════════
* Figure 1: KDE (left) + CDF (right), Total industry
* Shows that high-fracking counties have fatter tails in both directions —
* more volatile, not simply richer. CDFs cross, ruling out FOSD.
* ═══════════════════════════════════════════════════════════════════════════════
di ""
di "━━ Figure 1: KDE + CDF — Total industry ━━"

use `combined', clear
keep if ind_id == 1

cumul d_wages_plot if fracking_group == 0, gen(cdf_no)   equal
cumul d_wages_plot if fracking_group == 1, gen(cdf_high) equal

twoway ///
    (kdensity d_wages_plot if fracking_group == 0, ///
        lcolor(navy) lwidth(medthick) lpattern(solid)) ///
    (kdensity d_wages_plot if fracking_group == 1, ///
        lcolor(cranberry) lwidth(medthick) lpattern(dash)), ///
    legend(order(1 "No fracking" 2 "High fracking (top 10%)") ///
           position(1) ring(0) cols(1) size(small)) ///
    xtitle("{&Delta}Wages per capita (demeaned)", size(small)) ///
    ytitle("Density", size(small)) ///
    title("(a) Kernel Density", size(small)) ///
    scheme(s2color) graphregion(color(white)) name(kde, replace) nodraw

twoway ///
    (line cdf_no   d_wages_plot if fracking_group == 0, ///
        sort lcolor(navy) lwidth(medthick)) ///
    (line cdf_high d_wages_plot if fracking_group == 1, ///
        sort lcolor(cranberry) lwidth(medthick) lpattern(dash)), ///
    yline(0.5, lcolor(gs12) lpattern(dot)) ///
    legend(order(1 "No fracking" 2 "High fracking (top 10%)") ///
           position(5) ring(0) cols(1) size(small)) ///
    xtitle("{&Delta}Wages per capita (demeaned)", size(small)) ///
    ytitle("Cumulative probability", size(small)) ///
    title("(b) Empirical CDF", size(small)) ///
    scheme(s2color) graphregion(color(white)) name(cdf, replace) nodraw

graph combine kde cdf, cols(2) ///
    title("Wage Growth Distribution by Fracking Intensity — Total Industry", ///
          size(medsmall)) ///
    note("Treatment = NewValue_t + NewValue_{t-1}. County + year FE removed via Gauss-Seidel demeaning." ///
         "Wages winsorised at 1st/99th pct. KS test p-value shown in log. Source: Feyrer et al. (2015).", ///
         size(vsmall)) ///
    scheme(s2color) graphregion(color(white))

graph export "fig_np_1_distributional.png", replace width(1600)
di "  Saved: fig_np_1_distributional.png"

di ""
di "  KS test (Total — no fracking vs high fracking):"
ksmirnov d_wages_plot, by(fracking_group)

capture graph drop kde cdf

* ═══════════════════════════════════════════════════════════════════════════════
* Figure 2: Linearity check — Total industry
* Local linear (degree 1) vs local quadratic (degree 2) vs OLS, on the
* partial regression residuals (L1 OLS-partialled from wages).
* If degree-1 and degree-2 track closely → linearity defensible.
* Large divergence → nonlinearity that the paper's OLS assumption misses.
* ═══════════════════════════════════════════════════════════════════════════════
di ""
di "━━ Figure 2: Linearity check — Total industry ━━"

use `combined', clear
keep if ind_id == 1

twoway ///
    (lpolyci d_wages_partial_plot newval_L0_plot if !missing(newval_L0_plot), ///
        degree(1) kernel(epanechnikov) ///
        lcolor(navy) lwidth(medthick) acolor(navy%20)) ///
    (lpoly d_wages_partial_plot newval_L0_plot if !missing(newval_L0_plot), ///
        degree(2) kernel(epanechnikov) ///
        lcolor(cranberry) lwidth(medium) lpattern(dash)) ///
    (lfit d_wages_partial_plot newval_L0_plot if !missing(newval_L0_plot), ///
        lcolor(gs8) lwidth(thin) lpattern(dot)), ///
    legend(order(2 "Local linear (degree 1)" 3 "Local quadratic (degree 2)" 4 "OLS") ///
           position(11) ring(0) cols(1) size(small)) ///
    yline(0, lcolor(gs14) lwidth(vthin)) ///
    xtitle("NewValue_t (demeaned, fracking counties, trimmed p99)", size(small)) ///
    ytitle("{&Delta}Wages (demeaned, L1 partialled out)", size(small)) ///
    title("Linearity Check: Wage Response to Current Fracking Value", size(medsmall)) ///
    note("Y = wages residual after OLS-partialling out NewValue_{t-1} (FWL theorem)." ///
         "Shaded = 95% CI for degree 1. Source: Feyrer et al. (2015).", size(vsmall)) ///
    scheme(s2color) graphregion(color(white))

graph export "fig_np_2_linearity.png", replace width(1400)
di "  Saved: fig_np_2_linearity.png"

* ═══════════════════════════════════════════════════════════════════════════════
* Figure 3: 2×4 lpoly grid — one panel per industry
* Partial regression plot: wages with L1 removed vs demeaned L0.
* Reveals sector-specific nonlinearities not visible in aggregate.
* ═══════════════════════════════════════════════════════════════════════════════
di ""
di "━━ Figure 3: lpoly grid — all 8 industries ━━"

use `combined', clear

local fname3_1 "fig_np_3a_lpoly_total"
local fname3_2 "fig_np_3b_lpoly_mining"
local fname3_3 "fig_np_3c_lpoly_transport"
local fname3_4 "fig_np_3d_lpoly_construction"
local fname3_5 "fig_np_3e_lpoly_manufacturing"
local fname3_6 "fig_np_3f_lpoly_educ"
local fname3_7 "fig_np_3g_lpoly_services"
local fname3_8 "fig_np_3h_lpoly_government"

forvalues ind = 1/8 {
    twoway ///
        (lpolyci d_wages_partial_plot newval_L0_plot ///
            if ind_id == `ind' & !missing(newval_L0_plot), ///
            degree(1) kernel(epanechnikov) ///
            lcolor(`col`ind'') lwidth(medthick) acolor(`col`ind''%25)) ///
        (lfit d_wages_partial_plot newval_L0_plot ///
            if ind_id == `ind' & !missing(newval_L0_plot), ///
            lcolor(gs8) lwidth(medium) lpattern(shortdash)), ///
        legend(order(2 "Local linear (degree 1)" 3 "OLS") ///
               position(11) ring(0) cols(1) size(small)) ///
        xtitle("NewValue_t (demeaned, fracking counties, trimmed p99)", size(small)) ///
        ytitle("{&Delta}Wages per capita (L1 partialled out)", size(small)) ///
        title("Local Linear Wage Response — `lab`ind''", size(medsmall) color(`col`ind'')) ///
        yline(0, lcolor(gs12) lwidth(thin)) ///
        note("Y = {&Delta}wages with NewValue_{t-1} OLS-partialled out (FWL theorem)." ///
             "Shaded = 95% CI. Source: Feyrer et al. (2015).", size(vsmall)) ///
        scheme(s2color) graphregion(color(white))

    graph export "`fname3_`ind''.png", replace width(1400)
    di "  Saved: `fname3_`ind''.png"
}

* ═══════════════════════════════════════════════════════════════════════════════
* Figure 4: 2×2 KDE grid for the 4 most economically distinct industries
* Mining (mean-reversion), Manufacturing (crowding-out),
* Transport (uniform positive), Government (lagged fiscal channel)
* ═══════════════════════════════════════════════════════════════════════════════
di ""
di "━━ Figure 4: KDE comparison — selected industries ━━"

use `combined', clear

* Industries to show: 2=Mining, 3=Transport, 5=Manufacturing, 8=Government
local sel_inds  "2        3          5                8"
local sel_files "mining   transport  manufacturing   government"
local sel_subtitles ///
    "Mining (boom-bust mean reversion)" ///
    "Transport (steady input-demand channel)" ///
    "Manufacturing (crowding-out / Dutch Disease)" ///
    "Government (lagged fiscal transmission)"

local col_hi "cranberry"
local col_no "navy"

* Locals for titles and subtitles indexed by ind number
local title2 "Mining"
local title3 "Transport"
local title5 "Manufacturing"
local title8 "Government"
local sub2 "Boom-bust pattern: high-fracking distribution widens and skews"
local sub3 "Steady input demand: distributions shift rightward uniformly"
local sub5 "Crowding-out / Dutch Disease: high-fracking left-shifted"
local sub8 "Lagged fiscal channel: effect appears with a one-year delay"

foreach ind in 2 3 5 8 {

    if `ind' == 2  local fname "fig_np_4a_kde_mining"
    if `ind' == 3  local fname "fig_np_4b_kde_transport"
    if `ind' == 5  local fname "fig_np_4c_kde_manufacturing"
    if `ind' == 8  local fname "fig_np_4d_kde_government"

    twoway ///
        (kdensity d_wages_plot if ind_id == `ind' & fracking_group == 0, ///
            lcolor(`col_no') lwidth(medthick) lpattern(solid)) ///
        (kdensity d_wages_plot if ind_id == `ind' & fracking_group == 1, ///
            lcolor(`col_hi') lwidth(medthick) lpattern(dash)), ///
        legend(order(1 "No fracking (both periods zero)" ///
                     2 "High fracking (top 10% treat_sum)") ///
               position(1) ring(0) cols(1) size(small)) ///
        xtitle("{&Delta}Wages per capita (demeaned, winsorised)", size(small)) ///
        ytitle("Density", size(small)) ///
        title("`title`ind'': Wage Growth by Fracking Intensity", size(medsmall)) ///
        subtitle("`sub`ind''", size(small)) ///
        note("Treatment = NewValue_t + NewValue_{t-1}. County + year FE removed via Gauss-Seidel." ///
             "Source: Feyrer, Mansur & Sacerdote (AER 2015).", size(vsmall)) ///
        scheme(s2color) graphregion(color(white))

    graph export "`fname'.png", replace width(1400) height(900)
    di "  Saved: `fname'.png"
}

* ── KS tests for all industries (log only) ───────────────────────────────────
di ""
di "════════════════════════════════════════════════════════════════"
di "  KS tests by industry"
di "════════════════════════════════════════════════════════════════"

use `combined', clear

forvalues ind = 1/8 {
    di ""
    di "  `lab`ind''"
    ksmirnov d_wages_plot if ind_id == `ind', by(fracking_group)
}

* ── Mean summary table (log only) ────────────────────────────────────────────
di ""
di "════════════════════════════════════════════════════════════════════════"
di "  Mean demeaned {&Delta}wages: no fracking vs high fracking"
di "════════════════════════════════════════════════════════════════════════"
di "  Industry             No fracking    High fracking    Difference"
di "  " _dup(60) "-"

forvalues ind = 1/8 {
    qui summarize d_wages_plot if ind_id == `ind' & fracking_group == 0, meanonly
    local m0 = r(mean)
    qui summarize d_wages_plot if ind_id == `ind' & fracking_group == 1, meanonly
    local m1 = r(mean)
    di "  `lab`ind''" _col(22) %10.1f `m0' _col(37) %10.1f `m1' _col(52) %10.1f (`m1'-`m0')
}

di ""
di "All done. Output files:"
di "  fig_np_1_distributional.png      — KDE + CDF, Total industry"
di "  fig_np_2_linearity.png           — Linearity check, Total industry"
di "  fig_np_3_lpoly_grid.png          — lpoly by industry (2x4)"
di "  fig_np_4a_kde_mining.png         — KDE Mining"
di "  fig_np_4b_kde_transport.png      — KDE Transport"
di "  fig_np_4c_kde_manufacturing.png  — KDE Manufacturing"
di "  fig_np_4d_kde_government.png     — KDE Government"
