//! gromos++ `Arguments`: `@key value value …` on the command line (or in an `@f` file, expanded
//! by `gromos_io::gromos_args`), every key multi-valued, absent / flag / valued distinguished.

use std::collections::BTreeMap;

#[derive(Debug, Clone, Default)]
pub struct Arguments {
    program: String,
    values: BTreeMap<String, Vec<String>>,
    order: Vec<String>,
}

impl Arguments {
    /// Parse the process arguments; `known` are the accepted keys (an unknown key is an error,
    /// as in gromos++). `--help`/`-h` returns `Err(usage)`.
    pub fn parse(known: &[&str], usage: &str) -> Result<Self, String> {
        Self::from_vec(crate::gromos_args(), known, usage)
    }

    pub fn from_vec(raw: Vec<String>, known: &[&str], usage: &str) -> Result<Self, String> {
        let program = raw.first().cloned().unwrap_or_default();
        let mut a = Arguments {
            program,
            ..Default::default()
        };
        let mut key: Option<String> = None;
        for tok in raw.iter().skip(1) {
            if tok == "--help" || tok == "-h" {
                return Err(usage.to_string());
            }
            if let Some(k) = tok.strip_prefix("--").or_else(|| tok.strip_prefix('@')) {
                if !known.contains(&k) {
                    return Err(format!("unknown argument @{k}\n{usage}"));
                }
                a.values.entry(k.to_string()).or_default();
                if !a.order.contains(&k.to_string()) {
                    a.order.push(k.to_string());
                }
                key = Some(k.to_string());
            } else if let Some(k) = &key {
                a.values.get_mut(k).unwrap().push(tok.clone());
            } else {
                return Err(format!("stray argument '{tok}'\n{usage}"));
            }
        }
        Ok(a)
    }

    pub fn program(&self) -> &str {
        &self.program
    }

    /// gromos++ `count`: −1 when the key is absent, else the number of values (0 for a flag).
    pub fn count(&self, key: &str) -> i64 {
        self.values.get(key).map_or(-1, |v| v.len() as i64)
    }

    pub fn has(&self, key: &str) -> bool {
        self.values.contains_key(key)
    }

    /// All values of `key` (empty when absent).
    pub fn values(&self, key: &str) -> &[String] {
        self.values.get(key).map_or(&[], |v| v.as_slice())
    }

    /// The first value of `key`, or an error naming it (gromos++ `args[key]` / `check`).
    pub fn value(&self, key: &str) -> Result<&str, String> {
        self.values(key)
            .first()
            .map(|s| s.as_str())
            .ok_or_else(|| format!("@{key} is required"))
    }

    /// The first value parsed as `T`, or `default` when absent.
    pub fn get<T: std::str::FromStr>(&self, key: &str, default: T) -> Result<T, String> {
        match self.values(key).first() {
            None => Ok(default),
            Some(v) => v
                .parse::<T>()
                .map_err(|_| format!("@{key}: cannot parse '{v}'")),
        }
    }

    /// The first value parsed as `T`, required.
    pub fn require<T: std::str::FromStr>(&self, key: &str) -> Result<T, String> {
        let v = self.value(key)?;
        v.parse::<T>()
            .map_err(|_| format!("@{key}: cannot parse '{v}'"))
    }
}

/// Print `err` the gromos++ way and exit 1 (usage goes out unchanged).
pub fn fail(err: &str) -> ! {
    eprintln!("{err}");
    std::process::exit(1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keys_are_multi_valued_and_flags_count_zero() {
        let a = Arguments::from_vec(
            [
                "p", "--topo", "t.top", "--traj", "a.trc", "b.trc", "--nots", "@time", "0", "0.5",
            ]
            .iter()
            .map(|s| s.to_string())
            .collect(),
            &["topo", "traj", "nots", "time"],
            "usage",
        )
        .unwrap();
        assert_eq!(a.count("traj"), 2);
        assert_eq!(a.count("nots"), 0);
        assert_eq!(a.count("absent"), -1);
        assert_eq!(a.value("topo").unwrap(), "t.top");
        assert_eq!(a.get::<f64>("time", 1.0).unwrap(), 0.0);
        assert_eq!(a.values("time"), &["0", "0.5"]);
    }

    #[test]
    fn unknown_keys_are_rejected() {
        assert!(Arguments::from_vec(vec!["p".into(), "--nope".into()], &["topo"], "u").is_err());
    }
}
