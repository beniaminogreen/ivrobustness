//! Rust-only robustness values, using the existing generic covariance core.

use argmin::core::{CostFunction, Error, Executor, Gradient, State, TerminationReason};
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::quasinewton::LBFGS;
use num_dual::{Dual64, DualNum};

use crate::{CovarianceMatrix, ExtendedCovMatCore};

const TOL: f64 = 1.0e-10; // Same default tolerance as the R calc_rv function.
                          //
use extendr_api::prelude::*;

#[derive(Clone, Copy, Debug)]
pub(crate) enum Method {
    Soo,
    Iv,
    ExogenousIV,
    ExcludableIV,
}

impl Method {
    fn estimate(self, s: &ExtendedCovMatCore<f64>) -> f64 {
        match self {
            Self::Soo => s.soo_estimate(),
            Self::Iv | Self::ExogenousIV | Self::ExcludableIV => s.iv_estimate(),
        }
    }

    fn arrows<D: DualNum<Primitive = f64> + Copy>(self, s: &ExtendedCovMatCore<D>) -> [D; 2] {
        match self {
            Self::Soo => [s.soo_z(), s.soo_y()],
            Self::Iv => [s.iv_z(), s.iv_y()],
            Self::ExogenousIV => [s.iv_y(), D::from(0.0)],
            Self::ExcludableIV => [s.iv_z(), D::from(0.0)],
        }
    }

    fn assumption_gap<D: DualNum<Primitive = f64> + Copy>(self, s: &ExtendedCovMatCore<D>) -> D {
        match self {
            Self::Soo | Self::Iv => D::from(0.0),
            Self::ExogenousIV => s.iv_z(),
            Self::ExcludableIV => s.iv_y(),
        }
    }
}

#[derive(Debug, IntoList)]
pub(crate) struct RvResult {
    method: String,
    pub rv: f64,
    pub rho_1: f64,
    pub rho_2: f64,
    pub rho_3: f64,
    pub bias: f64,
    pub assumption_gap: f64,
}

const RHO_LIMITS: [f64; 3] = [0.99999, 0.9746794344808963, 0.99999];

fn to_rho<D: DualNum<Primitive = f64> + Copy>(x: &[D]) -> [D; 3] {
    std::array::from_fn(|i| {
        let one = D::from(1.0);
        let two = D::from(2.0);

        // Keep exponential arguments nonpositive to avoid overflow.
        let t = if x[i].re() >= 0.0 {
            let e = (-two * x[i]).exp();
            (one - e) / (one + e)
        } else {
            let e = (two * x[i]).exp();
            (e - one) / (e + one)
        };

        D::from(RHO_LIMITS[i]) * t
    })
}

#[derive(Copy, Debug, Clone)]
struct Objective<'a> {
    covariance: &'a CovarianceMatrix,
    method: Method,
    estimate: f64,
    target_bias: f64,
    penalty: f64,
}

impl Objective<'_> {
    fn evaluate<D: DualNum<Primitive = f64> + Copy>(&self, x: &[D]) -> D {
        let rho = to_rho(x);
        let s = self.covariance.extend_core(rho);

        let [a, b] = self.method.arrows(&s).map(|v| v.powi(2));
        let rv = if a.re() >= b.re() { a } else { b };

        let gap = D::from(self.estimate - self.target_bias) - s.tau();
        let assumption_gap = self.method.assumption_gap(&s);

        let ridge = rho.iter().fold(D::from(0.0), |sum, p| sum + p.powi(2));
        rv + D::from(self.penalty) * (gap.powi(2) + assumption_gap.powi(2)) + D::from(0.01) * ridge
    }
}

impl CostFunction for Objective<'_> {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, x: &Self::Param) -> Result<f64, Error> {
        let value = self.evaluate(x);
        assert!(value.is_finite(), "Non-finite RV objective");
        Ok(value)
    }
}

impl Gradient for Objective<'_> {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;

    fn gradient(&self, x: &Self::Param) -> Result<Vec<f64>, Error> {
        Ok((0..3)
            .map(|active| {
                let dual_x: Vec<Dual64> = x
                    .iter()
                    .enumerate()
                    .map(|(i, &v)| Dual64::new(v, if i == active { 1.0 } else { 0.0 }))
                    .collect();
                let value = self.evaluate(&dual_x);
                assert!(
                    value.re.is_finite() && value.eps.is_finite(),
                    "Non-finite RV automatic derivative"
                );
                value.eps
            })
            .collect())
    }
}

fn minimize(objective: Objective<'_>, start: Vec<f64>) -> Vec<f64> {
    // Warm start plus all eight sign combinations of these magnitudes.
    let starts = std::iter::once(start).chain((0..8).map(|mask| {
        [0.2, 0.4, 0.6]
            .iter()
            .enumerate()
            .map(|(j, &v)| if mask & (1 << j) == 0 { v } else { -v })
            .collect::<Vec<f64>>()
    }));

    let mut best = None;
    let mut best_cost = f64::INFINITY;
    let mut last_failure = String::from("No finite solution");

    for initial in starts {
        let linesearch = MoreThuenteLineSearch::new()
            .with_bounds(1.0e-20, f64::INFINITY)
            .unwrap();

        let solver: LBFGS<_, Vec<f64>, Vec<f64>, f64> = LBFGS::new(linesearch, 15)
            .with_tolerance_grad(1.0e-9)
            .unwrap()
            .with_tolerance_cost(1.0e-12)
            .unwrap();

        let result = match Executor::new(objective, solver)
            .configure(|state| state.param(initial).max_iters(100_000))
            .run()
        {
            Ok(result) => result,
            Err(error) => {
                last_failure = error.to_string();
                continue;
            }
        };

        let state = result.state();
        if !matches!(
            state.get_termination_reason(),
            Some(TerminationReason::SolverConverged | TerminationReason::TargetCostReached)
        ) {
            last_failure = format!("{:?}", state.get_termination_reason());
            continue;
        }

        let Some(x) = state.get_best_param() else {
            continue;
        };
        let cost = state.get_best_cost();

        // Prefer the earlier start when costs agree within solver tolerance.
        if cost.is_finite() && x.iter().all(|v| v.is_finite()) && cost < best_cost - 1.0e-12 {
            best_cost = cost;
            best = Some(x.clone());
        }
    }

    best.unwrap_or_else(|| {
        panic!("All RV optimization starts failed; last failure: {last_failure}")
    })
}

/// Target bias is signed: selected_estimate - tau = target_bias.
/// Panics on invalid input, solver failure, or failure to meet the tolerance.
pub(crate) fn calc_rv(covariance: &CovarianceMatrix, target_bias: f64, method: Method) -> RvResult {
    assert!(target_bias.is_finite(), "Target bias must be finite");
    assert!(
        covariance.mat.iter().all(|v| v.is_finite()) && covariance.mat.cholesky().is_some(),
        "Covariance matrix must be finite and positive definite"
    );
    let baseline = covariance.extend_core([0.0_f64; 3]);
    if matches!(method, Method::Iv) {
        assert!(
            baseline.mat[(4, 4)].is_finite() && baseline.mat[(4, 4)] > 0.0,
            "IV requires nonzero covariance between W and Z"
        );
    }
    let estimate = method.estimate(&baseline);
    assert!(estimate.is_finite(), "Undefined regression estimate");

    // Zero SOO bias has the exact solution rho=0 and RV=0.
    if matches!(method, Method::Soo) && target_bias == 0.0 {
        return RvResult {
            rv: 0.0,
            rho_1: 0.0,
            rho_2: 0.0,
            rho_3: 0.0,
            bias: 0.0,
            assumption_gap: 0.0,
            method: format!("{method:?}"),
        };
    }

    // Asymmetric, nonzero start avoids SOO's stationary point at the origin.
    let mut x = vec![0.05, 0.1, 0.15];
    let mut gap = f64::INFINITY;
    for stage in 0..14 {
        let previous_rho = to_rho(&x);
        let objective = Objective {
            covariance,
            method,
            estimate,
            target_bias,
            penalty: 100.0 * 10.0_f64.powi(stage),
        };

        x = minimize(objective, x);
        let rho = to_rho(&x);

        let s = covariance.extend_core(rho);
        let bias = estimate - s.tau();

        let assumption_gap = method.assumption_gap(&s).abs();

        gap = (bias - target_bias).abs();
        let movement = (0..3)
            .map(|i| (rho[i] - previous_rho[i]).abs())
            .fold(0.0, f64::max);
        if gap <= TOL && movement <= TOL && assumption_gap <= TOL {
            let [a, b] = method.arrows(&s);
            let rv = a.powi(2).max(b.powi(2));
            assert!(rv.is_finite(), "Non-finite robustness value");
            return RvResult {
                rv,
                rho_1: rho[0],
                rho_2: rho[1],
                rho_3: rho[2],
                bias,
                assumption_gap,
                method: format!("{method:?}"),
            };
        }
    }

    panic!("RV did not converge within tolerance {TOL}; final bias gap = {gap}");
}
