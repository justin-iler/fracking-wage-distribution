/*
================================================================================
Quantile Regression with Lag Structure — By Industry
Feyrer, Mansur & Sacerdote (AER 2015)

Specification: ΔWages_it = β0·NewValue_it + β1·NewValue_{i,t-1} + α_i + ω_t + ε_it
  County + year FE absorbed via Gauss-Seidel demeaning (Canay 2011).
  Quantiles: 0.10, 0.50, 0.90.  Bootstrap SE: 100 reps.
  NOTE: No IV. Expect slightly higher magnitudes than paper's IV results.

Produces:
  quantile_lagged_results_byindustry.csv  — full tables (8 industries)
  fig_qr_coefficients.png                 — dot plot of L0+L1 sums by industry
================================================================================
*/

clear all
set more off
set seed 12345

local lab1 "Total"
local lab2 "Mining"
local lab3 "Transport"
local lab4 "Construction"
local lab5 "Manufacturing"
local lab6 "Educ/Health"
local lab7 "Other Services"
local lab8 "Government"

local tol     = 1e-8
local maxiter = 100
local breps   = 100
local quantiles "10 50 90"

* ── Clean up previous outputs ────────────────────────────────────────────────
capture erase "quantile_lagged_results_byindustry.csv"
capture erase "quantile_lagged_results_total.csv"
capture erase "fig_qr_coefficients.png"

* ── Load data, create lag ─────────────────────────────────────────────────────
use "/Users/jiler20/Desktop/Nonparametric/Data/BLS_IRS_fossil_working.dta", clear
keep if sample == 1

sort fips id_industry year
by fips id_industry (year): gen newvalue_capita_L1 = newvalue_capita[_n-1]

drop if year == 2004
drop if missing(d_wages_capita, newvalue_capita, newvalue_capita_L1)

tempfile base
save `base'

* ── Matrix to collect L0+L1 sums for the coefficient plot ────────────────────
* Rows = 8 industries, cols = Q10 / Q50 / Q90 / OLS
matrix SUMS = J(8, 4, .)

* ── Main loop ────────────────────────────────────────────────────────────────
forvalues ind = 1/8 {

    di ""
    di "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    di "  Industry `ind': `lab`ind''"
    di "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    use `base', clear
    keep if id_industry == `ind'

    egen county_id = group(fips)
    egen year_id   = group(year)

    * Gauss-Seidel demeaning
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

    * Quantile regressions
    foreach q of local quantiles {
        local tau = `q' / 100
        qui bsqreg d_wages_capita_dm newvalue_capita_dm newvalue_capita_L1_dm, ///
            quantile(`tau') reps(`breps')
        est store qr_i`ind'_q`q'
        matrix SUMS[`ind', cond("`q'"=="10",1,cond("`q'"=="50",2,3))] = ///
            _b[newvalue_capita_dm] + _b[newvalue_capita_L1_dm]
    }

    * OLS
    qui reg d_wages_capita_dm newvalue_capita_dm newvalue_capita_L1_dm, robust
    est store ols_i`ind'
    matrix SUMS[`ind', 4] = _b[newvalue_capita_dm] + _b[newvalue_capita_L1_dm]

    * Print L0+L1 sums to log
    di "  L0+L1 sums:"
    di "    Q0.10 = " %10.1f SUMS[`ind',1]
    di "    Q0.50 = " %10.1f SUMS[`ind',2]
    di "    Q0.90 = " %10.1f SUMS[`ind',3]
    di "    OLS   = " %10.1f SUMS[`ind',4]

    * Export table (append industries)
    local append_mode = cond(`ind'==1, "replace", "append")
    esttab qr_i`ind'_q10 qr_i`ind'_q50 qr_i`ind'_q90 ols_i`ind' ///
        using "quantile_lagged_results_byindustry.csv", ///
        keep(newvalue_capita_dm newvalue_capita_L1_dm) ///
        mtitles("Q0.10" "Q0.50" "Q0.90" "OLS") ///
        coeflabels(newvalue_capita_dm    "NewValue_t (L0)" ///
                   newvalue_capita_L1_dm "NewValue_{t-1} (L1)") ///
        b(%12.4f) se(%12.4f) star(* 0.10 ** 0.05 *** 0.01) ///
        title("Industry: `lab`ind''") `append_mode'

} // end industry loop

* ── Figure: dot plot of L0+L1 sums ──────────────────────────────────────────
* Convert matrix to dataset for plotting

di ""
di "━━ Figure: Coefficient plot ━━"

drop _all
svmat SUMS

* svmat names columns SUMS1, SUMS2, ... when matrix is named SUMS
* rename robustly regardless of naming convention
ds
local allvars `r(varlist)'
local v1 : word 1 of `allvars'
local v2 : word 2 of `allvars'
local v3 : word 3 of `allvars'
local v4 : word 4 of `allvars'
rename `v1' q10
rename `v2' q50
rename `v3' q90
rename `v4' ols

* Industry labels as numeric
gen ind = _n
label define indlbl 1 "Total" 2 "Mining" 3 "Transport" 4 "Construct" ///
                    5 "Manufact" 6 "Educ/Health" 7 "Other Svcs" 8 "Govt"
label values ind indlbl

* Add zero reference
gen zero = 0

* Plot
twoway ///
    (scatter ind q10, msymbol(O) mcolor(navy) msize(medium)) ///
    (scatter ind q50, msymbol(D) mcolor(maroon) msize(medium)) ///
    (scatter ind q90, msymbol(T) mcolor(forest_green) msize(medium)) ///
    (scatter ind ols, msymbol(S) mcolor(dkorange) msize(medium)) ///
    (pcspike ind zero ind q90, lcolor(gs12) lwidth(thin)) ///
    (pcspike ind zero ind q10, lcolor(gs12) lwidth(thin)), ///
    xline(0, lcolor(gs10) lpattern(dash)) ///
    ylabel(1 "Total" 2 "Mining" 3 "Transport" 4 "Construct" ///
           5 "Manufact" 6 "Educ/Health" 7 "Other Svcs" 8 "Govt", ///
           angle(0) labsize(small)) ///
    ytitle("") ///
    xtitle("L0 + L1 coefficient sum ($/person per M$/person)", size(small)) ///
    legend(order(1 "Q0.10" 2 "Q0.50" 3 "Q0.90" 4 "OLS") ///
           position(3) ring(0) cols(1) size(small)) ///
    title("Total Wage Effect of Fracking by Industry and Quantile", size(medsmall)) ///
    subtitle("β_L0 + β_L1 sum. No IV; OLS-QR with Canay (2011) FE demeaning.", size(small)) ///
    note("Paper's IV benchmark (Total, OLS): ~34,000. Horizontal line at zero." ///
         "Source: Feyrer et al. (2015).", size(vsmall)) ///
    scheme(s2color) graphregion(color(white))

graph export "fig_qr_coefficients.png", replace width(1400)
di "  Saved: fig_qr_coefficients.png"

* ══════════════════════════════════════════════════════════════════════════════
* STEP 1: Comparison table — Our OLS (L0) vs Paper's Table 2 Panel A (county)
* ══════════════════════════════════════════════════════════════════════════════
* The paper's Panel A OLS county-level estimates (β, SE) from Table 2:
*   Total       33,957  (9,655)
*   Mining      16,932  (2,220)
*   Transport    7,980  (2,691)
*   Construct    4,304  (1,937)
*   Manufact      -758  (1,369)
*   Educ/Health -1,255  (1,524)
*   Govt         3,227  (1,004)
*   Oth. serv.   3,227  (1,004)
*
* Our L0 OLS (from stored estimates) should replicate these closely.
* We add Q0.10, Q0.50, Q0.90 to show distributional heterogeneity.
* ──────────────────────────────────────────────────────────────────────────────

di ""
di "━━ Step 1: Comparison Table — Our L0 QR+OLS vs Paper Table 2 Panel A ━━"

* Build a single esttab pulling only the L0 (newvalue_capita_dm) coefficient
* from all stored estimates, plus hardcoded paper benchmarks

* Print formatted table to log
di ""
di "═══════════════════════════════════════════════════════════════════════════════════"
di "  Table: L0 Contemporaneous Effect on Wages — OLS Replication + QR Extension"
di "  Spec: ΔWages_it = β0·NV_it + β1·NV_{t-1} + α_i + ω_t  (county level, OLS)"
di "═══════════════════════════════════════════════════════════════════════════════════"
di %15s "Industry" %12s "Paper OLS" %12s "Our OLS" %12s "Q0.10" %12s "Q0.50" %12s "Q0.90"
di "───────────────────────────────────────────────────────────────────────────────────"

* Paper's Panel A OLS values (hardcoded from Table 2)
local paper_b  33957  16932  7980  4304  -758  -1255  3227  3227
local paper_se  9655   2220  2691  1937  1369   1524  1004  1004

* Retrieve our estimates from stored results
local ind_labels `" "Total" "Mining" "Transport" "Construct" "Manufact" "Educ/Health" "Oth.Serv" "Govt" "'

forvalues ind = 1/8 {
    local pval : word `ind' of `paper_b'
    local indname : word `ind' of `ind_labels'

    * Retrieve our stored OLS and QR L0 coefficients
    est restore ols_i`ind'
    local our_ols = _b[newvalue_capita_dm]

    est restore qr_i`ind'_q10
    local our_q10 = _b[newvalue_capita_dm]

    est restore qr_i`ind'_q50
    local our_q50 = _b[newvalue_capita_dm]

    est restore qr_i`ind'_q90
    local our_q90 = _b[newvalue_capita_dm]

    di %15s "`indname'" ///
       %12.0f `pval'    ///
       %12.0f `our_ols' ///
       %12.0f `our_q10' ///
       %12.0f `our_q50' ///
       %12.0f `our_q90'
}

di "───────────────────────────────────────────────────────────────────────────────────"
di "  Paper OLS: FMS Table 2 Panel A, county level, IV instrument not used here."
di "  Our OLS/QR: Two-way FE demeaning (Canay 2011), L0 coefficient only, 100 boot reps."
di "  Units: $/person per $M/person new fossil fuel production value."
di "═══════════════════════════════════════════════════════════════════════════════════"

* ── Export to CSV ─────────────────────────────────────────────────────────────
* Build a tidy dataset for export

drop _all

* Industry labels
local ilabs `" "Total" "Mining" "Transport" "Construction" "Manufacturing" "Educ/Health" "Other Services" "Government" "'
local paper_bvals  33957  16932  7980  4304  -758  -1255  3227  3227
local paper_sevals  9655   2220  2691  1937  1369   1524  1004  1004

* Create dataset
set obs 8
gen ind = _n
gen str20 industry = ""
gen paper_ols    = .
gen paper_se     = .
gen our_ols      = .
gen our_q10      = .
gen our_q50      = .
gen our_q90      = .

forvalues ind = 1/8 {
    local iname : word `ind' of `ilabs'
    local pb    : word `ind' of `paper_bvals'
    local pse   : word `ind' of `paper_sevals'

    replace industry  = "`iname'" in `ind'
    replace paper_ols = `pb'      in `ind'
    replace paper_se  = `pse'     in `ind'

    est restore ols_i`ind'
    replace our_ols = _b[newvalue_capita_dm] in `ind'

    est restore qr_i`ind'_q10
    replace our_q10 = _b[newvalue_capita_dm] in `ind'

    est restore qr_i`ind'_q50
    replace our_q50 = _b[newvalue_capita_dm] in `ind'

    est restore qr_i`ind'_q90
    replace our_q90 = _b[newvalue_capita_dm] in `ind'
}

* Replication check column: pct diff between our OLS and paper OLS
gen pct_diff = round((our_ols - paper_ols) / abs(paper_ols) * 100, 0.1)

export delimited using "table2_comparison.csv", replace
di "  Saved: table2_comparison.csv"

* ══════════════════════════════════════════════════════════════════════════════
* STEP 2: Quantile gradient figure — L0 coefficient by industry and quantile
*         (directly comparable to Table 2 Panel A, one dot per estimate)
* ══════════════════════════════════════════════════════════════════════════════

di ""
di "━━ Step 2: Quantile gradient figure (L0 only, vs Paper OLS) ━━"

* Add paper OLS values for reference line/scatter
gen paper_ols_k = paper_ols / 1000
gen our_ols_k   = our_ols   / 1000
gen our_q10_k   = our_q10   / 1000
gen our_q50_k   = our_q50   / 1000
gen our_q90_k   = our_q90   / 1000
gen zero        = 0

* Reverse industry order so Total is at top (ind=1 → top of y-axis)
gen y = 9 - ind

label define ylabel ///
    8 "Total"          ///
    7 "Mining"         ///
    6 "Transport"      ///
    5 "Construction"   ///
    4 "Manufacturing"  ///
    3 "Educ/Health"    ///
    2 "Other Services" ///
    1 "Government",    replace
label values y ylabel

* Horizontal CI-style spikes spanning Q0.10 to Q0.90
twoway ///
    (pcspike y our_q10_k y our_q90_k, lcolor(gs10) lwidth(medthick)) ///
    (scatter y paper_ols_k, msymbol(S) mcolor(cranberry) msize(medlarge) ///
        mlabel(paper_ols) mlabformat(%6.0f) mlabsize(vsmall) mlabposition(12)) ///
    (scatter y our_ols_k,   msymbol(D) mcolor(dkorange)  msize(medium)) ///
    (scatter y our_q10_k,   msymbol(O) mcolor(navy)      msize(medium) mfcolor(white)) ///
    (scatter y our_q50_k,   msymbol(O) mcolor(navy)      msize(medium)) ///
    (scatter y our_q90_k,   msymbol(T) mcolor(navy%70)   msize(medium)), ///
    xline(0, lcolor(gs8) lpattern(dash) lwidth(thin)) ///
    ylabel(1/8, valuelabel angle(0) labsize(small)) ///
    ytitle("") ///
    xtitle("L0 coefficient: ΔWages per capita ($ thousands)", size(small)) ///
    legend(order( ///
        2 "Paper OLS (Table 2, Panel A)" ///
        3 "Our OLS" ///
        4 "Our Q0.10" ///
        5 "Our Q0.50" ///
        6 "Our Q0.90") ///
        position(3) ring(0) cols(1) size(small)) ///
    title("Contemporaneous Wage Effect by Industry and Quantile", size(medsmall)) ///
    subtitle("L0 coefficient only. Horizontal bar spans Q0.10–Q0.90.", size(small)) ///
    note("Paper OLS = FMS Table 2 Panel A (county, IV not used). Our estimates: OLS-QR," ///
         "Canay (2011) two-way FE demeaning, 100 bootstrap reps, no IV." ///
         "Source: BLS_IRS_fossil_working.dta.", size(vsmall)) ///
    scheme(s2color) graphregion(color(white))

graph export "fig_table2_comparison.png", replace width(1400)
di "  Saved: fig_table2_comparison.png"

di ""
di "All done. Output files:"
di "  quantile_lagged_results_byindustry.csv  — full QR tables by industry"
di "  fig_qr_coefficients.png                 — L0+L1 sum dot plot"
di "  table2_comparison.csv                   — replication comparison table"
di "  fig_table2_comparison.png               — quantile gradient figure vs paper"
