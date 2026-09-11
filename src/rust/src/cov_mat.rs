use nalgebra::{Matrix3, Vector3};

fn validate_wzy(w: &[f64], z: &[f64], y: &[f64]) -> std::result::Result<usize, String> {
    let n = w.len();

    if n < 2 {
        return Err("at least two observations are required".into());
    }

    if z.len() != n || y.len() != n {
        return Err("w, z, and y must have equal lengths".into());
    }

    let contains_invalid_value = w
        .iter()
        .chain(z.iter())
        .chain(y.iter())
        .any(|value| !value.is_finite());

    if contains_invalid_value {
        return Err("w, z, and y cannot contain NA, NaN, or infinity".into());
    }

    Ok(n)
}

/// Sample covariance matrix for variables ordered W, Z, Y.
///
/// Uses the conventional n - 1 denominator.
pub fn covariance_matrix_wzy(
    w: &[f64],
    z: &[f64],
    y: &[f64],
) -> std::result::Result<Matrix3<f64>, String> {
    let n = validate_wzy(w, z, y)?;

    let mut mean = Vector3::zeros();

    for i in 0..n {
        mean += Vector3::new(w[i], z[i], y[i]);
    }

    mean /= n as f64;

    let mut covariance = Matrix3::zeros();

    for i in 0..n {
        let observation = Vector3::new(w[i], z[i], y[i]);
        let deviation = observation - mean;

        for row in 0..3 {
            for col in 0..3 {
                covariance[(row, col)] += deviation[row] * deviation[col];
            }
        }
    }

    Ok(covariance / (n as f64 - 1.0))
}

/// Unbiased weighted covariance matrix for variables ordered W, Z, Y.
///
/// Uses reliability weights with denominator:
///
///     sum(weights) - sum(weights^2) / sum(weights)
///
/// This matches R's:
///
///     cov.wt(..., method = "unbiased")
pub fn weighted_covariance_matrix_wzy(
    w: &[f64],
    z: &[f64],
    y: &[f64],
    weights: &[f64],
) -> std::result::Result<Matrix3<f64>, String> {
    let n = validate_wzy(w, z, y)?;

    if weights.len() != n {
        return Err("weights must have the same length as w, z, and y".into());
    }

    if weights
        .iter()
        .any(|weight| !weight.is_finite() || *weight < 0.0)
    {
        return Err("weights must be finite and non-negative".into());
    }

    let weight_sum: f64 = weights.iter().sum();

    if weight_sum <= 0.0 {
        return Err("at least one weight must be positive".into());
    }

    let weight_square_sum: f64 = weights.iter().map(|weight| weight * weight).sum();

    let denominator = weight_sum - weight_square_sum / weight_sum;

    if denominator <= 0.0 {
        return Err("at least two observations must have positive weight".into());
    }

    let mut mean = Vector3::zeros();

    for i in 0..n {
        let observation = Vector3::new(w[i], z[i], y[i]);
        mean += observation * weights[i];
    }

    mean /= weight_sum;

    let mut covariance = Matrix3::zeros();

    for i in 0..n {
        let observation = Vector3::new(w[i], z[i], y[i]);
        let deviation = observation - mean;

        for row in 0..3 {
            for col in 0..3 {
                covariance[(row, col)] += weights[i] * deviation[row] * deviation[col];
            }
        }
    }

    Ok(covariance / denominator)
}
