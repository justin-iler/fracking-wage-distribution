/*
Quantile Regression: Geographic Dispersion of Economic Shocks
Feyrer, Mansur & Sacerdote (AER 2015)

Specification: ΔWages_it = β0·NewValue_it + β1·L.NewValue_it + α_i + ω_t + ε_it
  - α_i : county fixed effects
  - ω_t : year fixed effects
  - Two-way FE absorbed via iterative within-group demeaning (Canay 2011)
  - Quantiles: 0.1, 0.5, 0.9
  - Bootstrap SE: 500 replications
*/

clear all
set more off

* ── 1. Load & filter data ────────────────────────────────────────────────────

use "BLS_IRS_fossil_working.dta", clear

* Keep total industry, in-sample observations
keep if industry == "1 Total"
keep if sample   == 1

* Encode identifiers as numeric
egen county_id = group(fips)
egen year_id   = group(year)

* ── 2. Generate lag of NewValue within county ────────────────────────────────
* Sort by county and year before creating the lag

sort county_id year
by county_id: gen L_newvalue_capita = newvalue_capita[_n-1]

* Drop obs missing any regression variable (including the new lag)
keep if !missing(d_wages_capita, newvalue_capita, L_newvalue_capita)

di "Observations (after lag): `=_N'"
di "Counties                : `=county_id[_N]'"
tab year

* ── 3. Two-way within-group demeaning ────────────────────────────────────────
* Demean y, current NewValue, AND lagged NewValue

foreach var of varlist d_wages_capita newvalue_capita L_newvalue_capita {
    gen `var'_dm = `var'
}

local tol     = 1e-8
local maxiter = 100

forvalues iter = 1/`maxiter' {

    gen _y_old  = d_wages_capita_dm
    gen _x0_old = newvalue_capita_dm
    gen _x1_old = L_newvalue_capita_dm

    * Demean by county
    foreach var of varlist d_wages_capita_dm newvalue_capita_dm L_newvalue_capita_dm {
        qui bysort county_id: egen _mean = mean(`var')
        qui replace `var' = `var' - _mean
        drop _mean
    }

    * Demean by year
    foreach var of varlist d_wages_capita_dm newvalue_capita_dm L_newvalue_capita_dm {
        qui bysort year_id: egen _mean = mean(`var')
        qui replace `var' = `var' - _mean
        drop _mean
    }

    * Check convergence
    qui gen _diff_y  = abs(d_wages_capita_dm  - _y_old)
    qui gen _diff_x0 = abs(newvalue_capita_dm  - _x0_old)
    qui gen _diff_x1 = abs(L_newvalue_capita_dm - _x1_old)
    qui summarize _diff_y,  meanonly
    local dy = r(max)
    qui summarize _diff_x0, meanonly
    local dx0 = r(max)
    qui summarize _diff_x1, meanonly
    local dx1 = r(max)
    drop _y_old _x0_old _x1_old _diff_y _diff_x0 _diff_x1

    if max(`dy', `dx0', `dx1') < `tol' {
        di "Demeaning converged at iteration `iter'"
        continue, break
    }
}

* ── 4. Quantile regressions (q = 0.1, 0.5, 0.9) ─────────────────────────────

foreach q in 10 50 90 {

    local tau = `q' / 100

    di ""
    di "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    di "Quantile q = `tau'"
    di "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    bsqreg d_wages_capita_dm newvalue_capita_dm L_newvalue_capita_dm, ///
        quantile(`tau') reps(500)

    est store qreg_q`q'
}

* ── 5. OLS for comparison ─────────────────────────────────────────────────────

di ""
di "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
di "OLS (comparison)"
di "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

reg d_wages_capita_dm newvalue_capita_dm L_newvalue_capita_dm, robust
est store ols

* ── 6. Summary table ──────────────────────────────────────────────────────────

di ""
di "══════════════════════════════════════════════════════════════════════"
di "  Results: ΔWages_it = β0·NewValue_it + β1·L.NewValue_it + FE"
di "══════════════════════════════════════════════════════════════════════"

esttab qreg_q10 qreg_q50 qreg_q90 ols, ///
    keep(newvalue_capita_dm L_newvalue_capita_dm) ///
    varlabels(newvalue_capita_dm   "NewValue (t)" ///
              L_newvalue_capita_dm "NewValue (t-1)") ///
    mtitles("Q0.1" "Q0.5" "Q0.9" "OLS") ///
    b(%12.1f) se(%12.1f) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    scalars("N Observations") ///
    title("Effect of new fossil fuel value on wages per capita (with lag)") ///
    note("Two-way FE (county + year) via iterative demeaning." ///
         "Bootstrap SE: 500 reps (bsqreg). OLS: HC1-robust SE." ///
         "Units: $/capita wages per $/capita new fossil fuel value.")

* ── 7. Save results ───────────────────────────────────────────────────────────

esttab qreg_q10 qreg_q50 qreg_q90 ols ///
    using "quantile_regression_results_lag.csv", ///
    keep(newvalue_capita_dm L_newvalue_capita_dm) ///
    varlabels(newvalue_capita_dm   "NewValue (t)" ///
              L_newvalue_capita_dm "NewValue (t-1)") ///
    mtitles("Q0.1" "Q0.5" "Q0.9" "OLS") ///
    b(%12.4f) se(%12.4f) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    replace

di ""
di "Results saved to: quantile_regression_results_lag.csv"
