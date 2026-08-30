//! gromos++ `gmath::Stat` log-exponential statistics (B.A. Berg, Comput. Phys. Comm. 153 (2003)
//! 397): averages of exponentials of large numbers kept in log space — used by `dfmult`.

/// ln( (e^a + e^b) ) without overflow.
fn ln_add(a: f64, b: f64) -> f64 {
    a.max(b) + (1.0 + (a.min(b) - a.max(b)).exp()).ln()
}

/// ln( e^a − e^b ) for a ≥ b (the sign is the caller's business).
fn ln_sub(a: f64, b: f64) -> f64 {
    a.max(b) + (1.0 - (a.min(b) - a.max(b)).exp()).ln()
}

/// gromos++ `lnexpave`: ln < exp(x) > over the values.
pub fn lnexp_ave(vals: &[f64]) -> f64 {
    let mut ave = vals[0];
    for &v in &vals[1..] {
        ave = ln_add(ave, v);
    }
    ave - (vals.len() as f64).ln()
}

/// gromos++ `lnexpcovariance`: ln |cov(exp X, exp Y)| with its sign.
pub fn lnexp_covariance(x: &[f64], y: &[f64]) -> Result<(f64, i32), String> {
    if x.len() != y.len() {
        return Err("can't calculate the covariance: unequal number of elements".into());
    }
    let mut cov = x[0] + y[0];
    for i in 1..x.len() {
        cov = ln_add(cov, x[i] + y[i]);
    }
    cov -= (x.len() as f64).ln();
    let ave_prod = lnexp_ave(x) + lnexp_ave(y);
    let sign = if cov < ave_prod { -1 } else { 1 };
    Ok((ln_sub(cov, ave_prod), sign))
}

/// gromos++ `lnexp_stat_ineff`: ln of the statistical inefficiency from the integrated
/// autocorrelation of exp X · exp Y (stops where the correlation goes negative).
pub fn lnexp_stat_ineff(x: &[f64], y: &[f64]) -> Result<f64, String> {
    if x.len() != y.len() {
        return Err(
            "can't calculate the statistical inefficiency: unequal number of elements".into(),
        );
    }
    let n = x.len();
    let mut si = 0.0f64;
    let mu_x = lnexp_ave(x);
    let mu_y = lnexp_ave(y);
    let (sigma2_xy, mut sign) = lnexp_covariance(x, y)?;
    if sigma2_xy.exp() == 0.0 {
        eprintln!("Variance zero! Unable to compute the statistical inefficiency");
    }
    let mut t = 1usize;
    let mut increment = 1usize;
    loop {
        let mut c = ln_add(x[0] + y[t], y[0] + x[t]);
        for k in 2..=(n - t) {
            let term = ln_add(x[k - 1] + y[k + t - 1], y[k - 1] + x[k + t - 1]);
            c = ln_add(c, term);
        }
        c -= (2.0 * (n - t) as f64).ln();
        let aveprod = mu_x + mu_y;
        if c < aveprod {
            sign *= -1;
        }
        c = ln_sub(c, aveprod);
        c -= sigma2_xy;
        if sign <= 0 {
            break;
        }
        let term = 2f64.ln() + c + (1.0 - t as f64 / n as f64).ln() + (increment as f64).ln();
        si = ln_add(si, term);
        t += increment;
        increment += 1;
        if t >= n - 1 {
            break;
        }
    }
    Ok(si)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_average_of_exponentials() {
        let v = [0.0, (2.0f64).ln(), (3.0f64).ln()]; // exp: 1, 2, 3 → mean 2
        assert!((lnexp_ave(&v) - 2f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn covariance_of_identical_series_is_the_variance() {
        let v: Vec<f64> = [1.0f64, 2.0, 4.0].iter().map(|x| x.ln()).collect();
        let (c, sign) = lnexp_covariance(&v, &v).unwrap();
        // <x²> − <x>² = 7 − 49/9
        assert_eq!(sign, 1);
        assert!((c.exp() - (7.0 - 49.0 / 9.0)).abs() < 1e-10);
    }
}
