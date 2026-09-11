use extendr_api::prelude::*;

mod calc_rv;
mod cov_mat;
mod regr;

use crate::cov_mat::{covariance_matrix_wzy, weighted_covariance_matrix_wzy};
use crate::regr::{part_r, regr};

#[allow(dead_code)]
const U: usize = 0;

const W: usize = 1;
const Z: usize = 2;
const Y: usize = 3;
const ZHAT: usize = 4;

use nalgebra::{Matrix3, Matrix4, Matrix5, Matrix5x4};

use num_dual::DualNum;

use crate::calc_rv::{calc_rv, Method, RvResult};

use rand::{rngs::SmallRng, Rng, SeedableRng};
use rand_distr::{Distribution, Exp};

fn new_weights_vec<R: Rng + ?Sized>(n: usize, rng: &mut R) -> Vec<f64> {
    assert!(n > 0, "Must have at least one observation");

    let exp = Exp::new(1.0).unwrap();
    let mut weights: Vec<f64> = (0..n).map(|_| exp.sample(rng)).collect();

    let total: f64 = weights.iter().sum();
    for weight in &mut weights {
        *weight /= total;
    }

    weights
}

#[allow(dead_code)]
#[extendr]
#[derive(Clone, Debug, Copy)]
struct CovarianceMatrix {
    mat: Matrix3<f64>,
}

impl CovarianceMatrix {
    fn extend_core<D>(&self, params: [D; 3]) -> ExtendedCovMatCore<D>
    where
        D: DualNum<Primitive = f64> + Copy,
    {
        assert!(
            params.iter().all(|&p| p >= -1.0 && p <= 1.0),
            "params must be finite and lie in [-1, 1]"
        );

        // S0 is fixed: factor it using ordinary f64 arithmetic.
        let lower_f64 = self
            .mat
            .cholesky()
            .expect("Could not factor covariance matrix")
            .l();

        // Lift constants into D with zero derivatives.
        let lower = lower_f64.map(|value| D::from(value));
        let one = D::from(1.0);

        let mut l = Matrix4::<D>::identity();
        l.fixed_view_mut::<3, 3>(1, 1).copy_from(&lower);

        let r0 = (one - params[0].powi(2)).sqrt();
        let r1 = (one - params[1].powi(2)).sqrt();
        let r2 = (one - params[2].powi(2)).sqrt();

        l[(0, 1)] = params[0];
        l[(0, 2)] = params[1] * r0;
        l[(0, 3)] = params[2] * r0 * r1;

        // Algebraically equivalent to sqrt(1 - sum(l[0, 1:4]^2)),
        // without cancellation from subtracting nearly equal values.
        l[(0, 0)] = r0 * r1 * r2;

        let su = l * l.transpose();
        let coef_zhat = regr(&su, &[W], Z)[0];

        let mut extended_l = Matrix5x4::<D>::zeros();
        extended_l.fixed_view_mut::<4, 4>(0, 0).copy_from(&l);

        for column in 0..4 {
            extended_l[(ZHAT, column)] = coef_zhat * l[(W, column)];
        }

        ExtendedCovMatCore {
            mat: extended_l * extended_l.transpose(),
        }
    }
}

#[allow(dead_code)]
#[extendr]
impl CovarianceMatrix {
    fn new(w: Vec<f64>, z: Vec<f64>, y: Vec<f64>) -> Self {
        let mat = covariance_matrix_wzy(&w, &z, &y).expect("Could Not Create Covaraiance Matrix");

        Self { mat }
    }

    fn iv_estimate(&self) -> f64 {
        // Cov(W, Y) / Cov(W, Z)
        self.mat[(0, 2)] / self.mat[(0, 1)]
    }

    fn soo_estimate(&self) -> f64 {
        // Coefficient on Z in the regression Y ~ Z + W.
        regr(&self.mat, &[1, 0], 2)[0]
    }

    fn to_r_mat(&self) -> RMatrix<f64> {
        RMatrix::new_matrix(3, 3, |row, col| self.mat[(row, col)])
    }

    fn extend(&self, params: [f64; 3]) -> ExtendedCovarianceMatrix {
        ExtendedCovarianceMatrix {
            core: self.extend_core(params),
        }
    }

    fn exogenous_iv_rv(&self, b_iv: f64) -> RvResult {
        calc_rv(&self, b_iv, Method::ExogenousIV)
    }

    fn excludable_iv_rv(&self, b_iv: f64) -> RvResult {
        calc_rv(&self, b_iv, Method::ExcludableIV)
    }

    fn iv_rv(&self, b_iv: f64) -> RvResult {
        calc_rv(&self, b_iv, Method::Iv)
    }

    fn soo_rv(&self, b_soo: f64) -> RvResult {
        calc_rv(&self, b_soo, Method::Soo)
    }

    fn calc_iv_rv_with_se(w: Vec<f64>, z: Vec<f64>, y: Vec<f64>, b_iv: f64) -> RvResult {
        let mat = covariance_matrix_wzy(&w, &z, &y).expect("Could Not Create Covaraiance Matrix");

        let cov_mat = Self { mat };

        let rv = calc_rv(&cov_mat, b_iv, Method::Iv);

        let mut rng = SmallRng::from_rng(&mut rand::rng());

        let mut bootstrap_rvs = Vec::new();
        for i in 0..1000 {
            let weights = new_weights_vec(y.len(), &mut rng);

            let mat = weighted_covariance_matrix_wzy(&w, &y, &z, &weights).unwrap();
            dbg!(i);
            dbg!(mat);
            let inner_cov_mat = Self { mat };

            let inner_rv = calc_rv(&inner_cov_mat, b_iv, Method::Iv).rv;

            bootstrap_rvs.push(inner_rv)
        }

        println!("bootstrap_rvs {:?}", bootstrap_rvs);

        rv
    }
}

struct ExtendedCovMatCore<D> {
    mat: Matrix5<D>,
}

#[allow(dead_code)]
impl<D> ExtendedCovMatCore<D>
where
    D: DualNum + Copy,
{
    fn soo_estimate(&self) -> D {
        regr(&self.mat, &[Z, W], Y)[0]
    }

    fn iv_estimate(&self) -> D {
        regr(&self.mat, &[ZHAT], Y)[0]
    }

    fn tau(&self) -> D {
        regr(&self.mat, &[Z, W, U], Y)[0]
    }

    fn soo_z(&self) -> D {
        part_r(&self.mat, &[U, Z, W])
    }

    fn soo_y(&self) -> D {
        part_r(&self.mat, &[U, Y, W, Z])
    }

    fn iv_z(&self) -> D {
        part_r(&self.mat, &[U, W])
    }

    fn iv_y(&self) -> D {
        part_r(&self.mat, &[Y, W, U, Z])
    }

    fn arrows(&self) -> [D; 4] {
        [self.soo_z(), self.soo_y(), self.iv_z(), self.iv_y()]
    }

    fn biases(&self) -> (D, D) {
        let tau = self.tau();

        (self.soo_estimate() - tau, self.iv_estimate() - tau)
    }
}

#[allow(dead_code)]
#[extendr]
struct ExtendedCovarianceMatrix {
    core: ExtendedCovMatCore<f64>,
}

#[extendr]
impl ExtendedCovarianceMatrix {
    fn to_r_mat(&self) -> RMatrix<f64> {
        RMatrix::new_matrix(5, 5, |row, col| self.core.mat[(row, col)])
    }

    fn soo_estimate(&self) -> f64 {
        self.core.soo_estimate()
    }

    fn iv_estimate(&self) -> f64 {
        self.core.iv_estimate()
    }

    fn tau(&self) -> f64 {
        self.core.tau()
    }

    fn soo_z(&self) -> f64 {
        self.core.soo_z()
    }

    fn soo_y(&self) -> f64 {
        self.core.soo_y()
    }

    fn iv_z(&self) -> f64 {
        self.core.iv_z()
    }

    fn iv_y(&self) -> f64 {
        self.core.iv_y()
    }
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod ivrobustness;
    impl CovarianceMatrix;
    impl ExtendedCovarianceMatrix;
}
