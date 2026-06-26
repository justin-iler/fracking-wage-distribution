/*
Nonparametric Replication of Figure 3:
"Wage Income Effects Including Neighbors within a Given Distance"
Feyrer, Mansur & Sacerdote (AER 2017)

Strategy:
  - Paper: plots OLS/IV β at each distance (0,20,...,200 miles)
  - Our version: plots QR β at Q0.1, Q0.5, Q0.9 at each distance
  - Both use two-way FE (county + year) via within-group demeaning
  - Final figure overlays paper's OLS line with our three quantile lines

The key insight: the paper's β collapses heterogeneity into a single mean.
Quantile regression at each distance reveals whether spillovers benefit
all counties equally or are concentrated in high-wage counties.

Run from working directory containing BLS_Distance_working.dta
*/

clear all
set more off

* ── 0. Install required packages ─────────────────────────────────────────────
* Uncomment if not already installed:
* ssc install estout, replace

* ── 1. Load data ──────────────────────────────────────────────────────────────

use "BLS_Distance_working.dta", clear

keep if industry == "1 Total"
keep if sample   == 1

egen county_id = group(fips)
egen year_id   = group(year)

* ── 2. Define distance bands and variable names ───────────────────────────────

* Distances matching the paper's Figure 3
local distances "0 20 40 60 80 100 120 140 160 180 200"

* ── 3. Two-way demeaning helper program ──────────────────────────────────────
* We demean y and x for each distance separately since the y and x variables
* change at each distance (different geographic aggregations)

capture program drop twoway_demean
program define twoway_demean
    * args: varname county_id year_id -> creates varname_dm
    args var g1 g2
    capture drop `var'_dm
    gen `var'_dm = `var'
    local tol     = 1e-8
    local maxiter = 100
    forvalues iter = 1/`maxiter' {
        qui gen _old_`var' = `var'_dm
        qui bysort `g1': egen _mn = mean(`var'_dm)
        qui replace `var'_dm = `var'_dm - _mn
        drop _mn
        qui bysort `g2': egen _mn = mean(`var'_dm)
        qui replace `var'_dm = `var'_dm - _mn
        drop _mn
        qui gen _diff_`var' = abs(`var'_dm - _old_`var')
        qui summarize _diff_`var', meanonly
        local maxd = r(max)
        drop _old_`var' _diff_`var'
        if `maxd' < `tol' {
            continue, break
        }
    }
end

* ── 4. Loop over distances: demean and run QR ─────────────────────────────────

* Store results in a matrix: rows = distances, cols = OLS + Q10 + Q50 + Q90
* Each stores: beta, CI lower, CI upper

local ndist : word count `distances'
matrix RESULTS = J(`ndist', 13, .)
* Cols: distance(1) | OLS_b OLS_l OLS_u(2-4) | Q10_b Q10_l Q10_u(5-7) | Q50_b Q50_l Q50_u(8-10) | Q90_b Q90_l Q90_u(11-13)

local row = 0

foreach d of local distances {

    local row = `row' + 1
    matrix RESULTS[`row', 1] = `d'

    local yvar "d_d`d'wages_cap"
    local xvar "d`d'newfossilval_cap"

    di ""
    di "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    di "  Distance = `d' miles"
    di "  y = `yvar'   x = `xvar'"
    di "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    * Skip if variables missing
    capture confirm variable `yvar'
    if _rc != 0 {
        di "  Variable `yvar' not found, skipping."
        continue
    }

    * Keep complete cases for this distance
    preserve
    keep if !missing(`yvar', `xvar')

    * Demean y and x
    twoway_demean `yvar' county_id year_id
    twoway_demean `xvar' county_id year_id

    * ── OLS ──
    qui reg `yvar'_dm `xvar'_dm, robust
    matrix RESULTS[`row', 2] = _b[`xvar'_dm]
    matrix RESULTS[`row', 3] = _b[`xvar'_dm] - 1.96*_se[`xvar'_dm]
    matrix RESULTS[`row', 4] = _b[`xvar'_dm] + 1.96*_se[`xvar'_dm]
    di "  OLS beta = " %10.1f _b[`xvar'_dm]

    * ── Quantile regressions (bsqreg, 200 reps for speed) ──
    foreach q in 10 50 90 {
        local tau = `q'/100
        qui bsqreg `yvar'_dm `xvar'_dm, quantile(`tau') reps(200)

        * Pre-compute column indices as locals — Stata cannot evaluate
        * arithmetic inside matrix[] subscripts directly
        if `q' == 10 {
            local col1 = 5
            local col2 = 6
            local col3 = 7
        }
        if `q' == 50 {
            local col1 = 8
            local col2 = 9
            local col3 = 10
        }
        if `q' == 90 {
            local col1 = 11
            local col2 = 12
            local col3 = 13
        }

        local beta = _b[`xvar'_dm]
        local se   = _se[`xvar'_dm]

        matrix RESULTS[`row', `col1'] = `beta'
        matrix RESULTS[`row', `col2'] = `beta' - 1.96*`se'
        matrix RESULTS[`row', `col3'] = `beta' + 1.96*`se'
        di "  Q`q' beta = " %10.1f `beta'
    }

    restore
}

* ── 5. Save results matrix to dataset ────────────────────────────────────────

clear
svmat RESULTS, names(col)

rename c1  distance
rename c2  ols_b
rename c3  ols_l
rename c4  ols_u
rename c5  q10_b
rename c6  q10_l
rename c7  q10_u
rename c8  q50_b
rename c9  q50_l
rename c10 q50_u
rename c11 q90_b
rename c12 q90_l
rename c13 q90_u

* Scale to thousands for readability (matching paper's y-axis)
foreach v of varlist ols_b ols_l ols_u q10_b q10_l q10_u ///
                     q50_b q50_l q50_u q90_b q90_l q90_u {
    replace `v' = `v' / 1000
}

save "fig3_qr_distance_results.dta", replace
di "Results saved to fig3_qr_distance_results.dta"

* ── 6. Merge in paper's IV results from distancefigures1.dta ──────────────────

merge 1:1 distance using "distancefigures1.dta", keepusing(wages_cap_iv_b wages_cap_iv_l wages_cap_iv_u)
drop _merge

* Scale paper's IV results to thousands
foreach v of varlist wages_cap_iv_b wages_cap_iv_l wages_cap_iv_u {
    replace `v' = `v' / 1000
}

* ── 7. Figure A: Quantile lines vs paper's OLS + IV ──────────────────────────

twoway ///
    (rcap ols_l ols_u distance, ///
        lcolor(gs12) lwidth(vthin)) ///
    (connected ols_b distance, ///
        lcolor(gs10) lwidth(thin) lpattern(shortdash) msymbol(none)) ///
    (connected wages_cap_iv_b distance, ///
        lcolor(black) lwidth(medthick) lpattern(solid) msymbol(square) msize(small)) ///
    (connected q10_b distance, ///
        lcolor(navy) lwidth(medium) lpattern(dash) msymbol(circle) msize(small)) ///
    (connected q50_b distance, ///
        lcolor(cranberry) lwidth(medium) lpattern(solid) msymbol(diamond) msize(small)) ///
    (connected q90_b distance, ///
        lcolor(dkorange) lwidth(medium) lpattern(longdash) msymbol(triangle) msize(small)), ///
    legend(order( ///
        2 "OLS (this paper)" ///
        3 "IV (paper, Fig. 3)" ///
        4 "QR: Q0.1" ///
        5 "QR: Q0.5" ///
        6 "QR: Q0.9") ///
        position(11) ring(0) cols(1) size(small)) ///
    xtitle("Distance from new production (miles)", size(small)) ///
    ytitle("ΔWage income per capita ($ thousands)", size(small)) ///
    xlabel(0(20)200) ///
    yline(0, lcolor(gs14) lpattern(dot)) ///
    title("Geographic Spillovers: OLS vs Quantile Regression", size(medsmall)) ///
    subtitle("Two-way FE removed. QR: 200 bootstrap reps.", size(small)) ///
    note("Dependent variable: ΔBLSwages per capita. Regressor: NewFossilValue per capita." ///
         "County + year FE via iterative demeaning. Source: Feyrer et al. (2017).", size(vsmall)) ///
    scheme(s2color) graphregion(color(white))

graph export "fig3a_qr_vs_paper.png", replace width(1400)
di "Saved: fig3a_qr_vs_paper.png"

* ── 8. Figure B: Quantile spread (Q90 - Q10) by distance ─────────────────────
* Shows whether distributional heterogeneity grows or shrinks with distance

gen qspread = q90_b - q10_b

twoway ///
    (bar qspread distance, ///
        barwidth(15) fcolor(navy%40) lcolor(navy) lwidth(thin)) ///
    (connected wages_cap_iv_b distance, ///
        yaxis(2) lcolor(cranberry) lwidth(medium) msymbol(square) msize(small)), ///
    legend(order(1 "Q0.9 - Q0.1 spread (left axis)" ///
                 2 "IV coefficient (right axis)") ///
           position(11) ring(0) cols(1) size(small)) ///
    xtitle("Distance from new production (miles)", size(small)) ///
    ytitle("QR spread: Q0.9 − Q0.1 ($ thousands)", size(small)) ///
    ytitle("IV β ($ thousands)", axis(2) size(small)) ///
    xlabel(0(20)200) ///
    title("Distributional Heterogeneity of Wage Spillovers by Distance", size(medsmall)) ///
    subtitle("Wider spread = more inequality in who benefits from fracking spillovers", size(small)) ///
    note("Source: BLS_Distance_working.dta. Feyrer et al. (2017).", size(vsmall)) ///
    scheme(s2color) graphregion(color(white))

graph export "fig3b_qr_spread.png", replace width(1400)
di "Saved: fig3b_qr_spread.png"

* ── 9. Summary table ──────────────────────────────────────────────────────────

di ""
di "══════════════════════════════════════════════════════════════════════════"
di "  Coefficients by distance ($ thousands per capita)"
di "  OLS = this paper; IV = paper's Figure 3; Q10/Q50/Q90 = quantile reg."
di "══════════════════════════════════════════════════════════════════════════"
di %5s "Dist" %10s "OLS" %10s "IV(paper)" %10s "Q0.1" %10s "Q0.5" %10s "Q0.9"
di "──────────────────────────────────────────────────────────────────────────"

forvalues i = 1/`=_N' {
    di %5.0f distance[`i'] ///
       %10.1f ols_b[`i'] ///
       %10.1f wages_cap_iv_b[`i'] ///
       %10.1f q10_b[`i'] ///
       %10.1f q50_b[`i'] ///
       %10.1f q90_b[`i']
}
di "══════════════════════════════════════════════════════════════════════════"

di ""
di "All output saved:"
di "  fig3_qr_distance_results.dta  — results dataset"
di "  fig3a_qr_vs_paper.png         — quantile lines vs OLS/IV"
di "  fig3b_qr_spread.png           — Q90-Q10 spread by distance"
