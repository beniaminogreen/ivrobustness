use nalgebra::{DMatrix, DVector, SMatrix};
use num_dual::{linalg::LU, DualNum};

/// Equivalent to:
/// solve(S[predictors, predictors], S[predictors, response])
pub(crate) fn regr<const N: usize, D>(
    covariance: &SMatrix<D, N, N>,
    predictors: &[usize],
    response: usize,
) -> DVector<D>
where
    D: DualNum + Copy,
{
    assert!(!predictors.is_empty(), "At least one predictor is required");
    assert!(response < N, "Response index is out of bounds");
    assert!(
        predictors.iter().all(|&index| index < N),
        "Predictor index is out of bounds"
    );

    let k = predictors.len();

    let predictor_covariance = DMatrix::from_fn(k, k, |row, column| {
        covariance[(predictors[row], predictors[column])]
    });

    let predictor_response = DVector::from_fn(k, |row, _| covariance[(predictors[row], response)]);

    let decomposition =
        LU::new(predictor_covariance).expect("Singular regression covariance matrix");

    decomposition.solve(&predictor_response)
}

pub(crate) fn part_r<const N: usize, D>(covariance: &SMatrix<D, N, N>, indices: &[usize]) -> D
where
    D: DualNum + Copy,
{
    assert!(indices.len() >= 2, "part_r requires at least two variables");
    assert!(
        indices.iter().all(|&index| index < N),
        "Variable index is out of bounds"
    );

    let size = indices.len();

    let selected_covariance = DMatrix::from_fn(size, size, |row, column| {
        covariance[(indices[row], indices[column])]
    });

    let precision = LU::new(selected_covariance)
        .expect("Could not factor selected covariance matrix")
        .inverse();

    -precision[(0, 1)] / (precision[(0, 0)] * precision[(1, 1)]).sqrt()
}
