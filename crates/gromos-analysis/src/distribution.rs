//! gromos++ `gmath::Distribution`: `nsteps` bins over `[begin, end)`, written as bin centre and
//! count (or count normalised to unit area).

#[derive(Debug, Clone)]
pub struct Distribution {
    begin: f64,
    end: f64,
    step: f64,
    counts: Vec<usize>,
    n_values: usize,
}

impl Distribution {
    pub fn new(begin: f64, end: f64, nsteps: usize) -> Result<Self, String> {
        if begin >= end {
            return Err("distribution: upper boundary should be higher than lower".into());
        }
        if nsteps < 1 {
            return Err("distribution: you need at least one step".into());
        }
        Ok(Distribution {
            begin,
            end,
            step: (end - begin) / nsteps as f64,
            counts: vec![0; nsteps],
            n_values: 0,
        })
    }

    /// gromos++ `Stat::dist_init(nsteps)`: the range of the values themselves.
    pub fn over_values(values: &[f64], nsteps: usize) -> Result<Self, String> {
        let min = values.iter().copied().fold(f64::INFINITY, f64::min);
        let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let mut d = if max == min {
            let dd = max * 0.01;
            Self::new(min - dd, max + dd, nsteps)?
        } else {
            Self::new(min, max + 1e-10, nsteps)?
        };
        for &v in values {
            d.add(v);
        }
        Ok(d)
    }

    /// gromos++ `Stat::dist_init(lower, upper, nsteps, periodic)`.
    pub fn over_values_in(
        values: &[f64],
        lower: f64,
        upper: f64,
        nsteps: usize,
        periodic: bool,
    ) -> Result<Self, String> {
        let mut d = Self::new(lower, upper, nsteps)?;
        let period = upper - lower;
        for &v in values {
            let mut x = v;
            if periodic {
                while x >= upper {
                    x -= period;
                }
                while x < lower {
                    x += period;
                }
            }
            d.add(x);
        }
        Ok(d)
    }

    pub fn add(&mut self, value: f64) {
        if value >= self.begin && value < self.end {
            let q = ((value - self.begin) / self.step) as usize;
            if q < self.counts.len() {
                self.counts[q] += 1;
                self.n_values += 1;
            }
        }
    }

    pub fn n_values(&self) -> usize {
        self.n_values
    }

    /// Centre of the most populated bin.
    pub fn max_val_at(&self) -> f64 {
        let mut x_max = 0;
        for (i, &c) in self.counts.iter().enumerate().skip(1) {
            if c > self.counts[x_max] {
                x_max = i;
            }
        }
        self.begin + (x_max as f64 + 0.5) * self.step
    }

    /// `setw(8) centre \t setw(5) count` per bin, numbers in the caller's stream format
    /// (`fmt`: fixed with n decimals, as `tser` sets, or the default `%g`).
    pub fn write(&self, fmt: &dyn Fn(f64) -> String) -> String {
        let mut o = String::new();
        for (i, c) in self.counts.iter().enumerate() {
            o.push_str(&format!(
                "{:>8}\t{:>5}\n",
                fmt(self.begin + (i as f64 + 0.5) * self.step),
                c
            ));
        }
        o
    }

    /// Counts divided by (number of values × bin width).
    pub fn write_normalized(&self, fmt: &dyn Fn(f64) -> String) -> String {
        let nval = self.n_values.max(1) as f64;
        let mut o = String::new();
        for (i, c) in self.counts.iter().enumerate() {
            o.push_str(&format!(
                "{:>8}\t{:>5}\n",
                fmt(self.begin + (i as f64 + 0.5) * self.step),
                fmt(*c as f64 / (nval * self.step))
            ));
        }
        o
    }
}

/// C++ default stream formatting of a double with `precision` significant digits (`%g`):
/// fixed when the exponent is in [−5, precision), scientific otherwise, trailing zeros dropped.
pub fn cpp_g(x: f64, precision: usize) -> String {
    if x == 0.0 {
        return "0".to_string();
    }
    if !x.is_finite() {
        return if x.is_nan() {
            "nan".into()
        } else if x > 0.0 {
            "inf".into()
        } else {
            "-inf".into()
        };
    }
    let p = precision.max(1);
    // exponent after rounding to p significant digits
    let sci = format!("{:.*e}", p - 1, x);
    let (mant, ex) = sci.split_once('e').unwrap();
    let exp: i32 = ex.parse().unwrap();
    if exp < -4 || exp >= p as i32 {
        let m = if mant.contains('.') {
            mant.trim_end_matches('0').trim_end_matches('.')
        } else {
            mant
        };
        return format!("{m}e{}{:02}", if exp < 0 { '-' } else { '+' }, exp.abs());
    }
    let decimals = (p as i32 - 1 - exp).max(0) as usize;
    let t = format!("{:.*}", decimals, x);
    if t.contains('.') {
        t.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        t
    }
}

/// Six significant digits — the stream default.
pub fn trim_float(x: f64) -> String {
    cpp_g(x, 6)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cpp_default_formatting() {
        assert_eq!(cpp_g(0.153, 6), "0.153");
        assert_eq!(cpp_g(1234567.0, 6), "1.23457e+06");
        assert_eq!(cpp_g(0.00001234, 6), "1.234e-05");
        assert_eq!(cpp_g(100.0, 6), "100");
        assert_eq!(cpp_g(-2.5, 5), "-2.5");
        assert_eq!(cpp_g(0.0001, 6), "0.0001");
    }
}
