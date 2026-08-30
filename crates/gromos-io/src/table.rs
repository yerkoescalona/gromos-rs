//! Numeric column files — the time series `ene_ana`, `tser` and friends write and the free-energy
//! programs read back: whitespace-separated numbers, `#` starts a comment, blank lines skipped.
//! gromos++ reads these with `istringstream` line by line; a line that does not parse is an error.

use std::path::Path;

use crate::IoError;

/// Every data row of `path` as numbers. Rows may differ in length; the caller checks what it needs.
pub fn read_columns<P: AsRef<Path>>(path: P) -> Result<Vec<Vec<f64>>, IoError> {
    let text = std::fs::read_to_string(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;
    parse_columns(&text, &path.as_ref().display().to_string())
}

/// [`read_columns`] on text already in memory.
pub fn parse_columns(text: &str, what: &str) -> Result<Vec<Vec<f64>>, IoError> {
    let mut rows = Vec::new();
    for (n, line) in text.lines().enumerate() {
        let data = line.split('#').next().unwrap_or("").trim();
        if data.is_empty() {
            continue;
        }
        let row: Result<Vec<f64>, _> = data.split_whitespace().map(|t| t.parse::<f64>()).collect();
        match row {
            Ok(r) => rows.push(r),
            Err(_) => {
                return Err(IoError::ParseError(format!(
                    "{what}, line {}: not a numeric row: {data:?}",
                    n + 1
                )))
            },
        }
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn comments_and_blank_lines_are_skipped() {
        let rows = parse_columns("# t  x\n\n 0.0 1.5 # trailing\n0.002 2.5\n", "test").unwrap();
        assert_eq!(rows, vec![vec![0.0, 1.5], vec![0.002, 2.5]]);
    }

    #[test]
    fn a_word_is_an_error() {
        assert!(parse_columns("0.0 abc\n", "test").is_err());
    }
}
