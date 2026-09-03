//! An energy trajectory read according to a library file's declaration of it — gromos++'s
//! `utils::EnergyTraj` plus the `gmath::Expression` evaluator it leans on.
//!
//! `ene_ana`, `ext_ti_ana` and friends do not know the layout of a `.tre`/`.trg` at compile
//! time. A library file declares it:
//!
//! ```text
//! ENERTRJ
//!   block ENERGY03
//!     subblock ENER 50 1
//!     size NUM_BATHS
//!     subblock KINENER NUM_BATHS 3
//!     subblock NONBONDED matrix_NUM_ENERGY_GROUPS 6
//! END
//! VARIABLES
//!   totpot = ENER[3]
//!   solute_intra = BONDED[1][1] + BONDED[1][2]
//! END
//! ```
//!
//! `block` names a GROMOS block in the file; the `subblock`s that follow consume its values in
//! order. `size` reads one integer from the stream and binds it as a dimension for later
//! subblocks — and simultaneously binds `matrix_<NAME>` to N(N+1)/2, which is how the
//! energy-group pair table is declared. `VARIABLES` then names arithmetic over those tables.
//! That indirection is what lets one library file span GROMOS versions, and why the trajectory
//! carries an `ENEVERSION` block to check against the library's.
//!
//! **Quirk reproduced deliberately:** gromos++ splits an operator chain at the *rightmost*
//! `+`/`-`, but at the *leftmost* `*`/`/` (`Expression.cc:310`), so `*` and `/` come out
//! right-associative: `a / b / c` evaluates as `a / (b / c)`, not `(a / b) / c`. The shipped
//! library's only multi-operator entry is `densit = MASS[1] * 1.66056 / VOLUME[1]`, where the
//! two groupings differ only in the last bits, but reproducing the rule keeps every digit of
//! our output equal to gromos++'s.

use crate::gz::TextReader;
use std::collections::HashMap;

/// A subblock dimension: a literal count, or a `size` bound while reading the frame.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Dim {
    Fixed(usize),
    Named(String),
}

/// One line of an `ENERTRJ`/`FRENERTRJ` declaration.
#[derive(Debug, Clone)]
enum Entry {
    /// A GROMOS block in the trajectory file; the entries after it consume its values.
    Block(String),
    /// Read one integer and bind it (and `matrix_<name>`) as a dimension.
    Size(String),
    /// A named `rows × cols` table of values.
    Sub { name: String, rows: Dim, cols: Dim },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Op {
    Add,
    Sub,
    Mul,
    Div,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Fun {
    Sin,
    Cos,
    Exp,
    Log,
}

#[derive(Debug, Clone)]
enum Expr {
    Const(f64),
    /// `SUBBLOCK[i]` / `SUBBLOCK[i][j]`, or a bare name resolved at evaluation time.
    Name {
        name: String,
        i: usize,
        j: usize,
    },
    Fun(Fun, Box<Expr>),
    Bin(Op, Box<Expr>, Box<Expr>),
}

/// An energy trajectory plus the library that describes it.
#[derive(Debug, Default)]
pub struct EnergyTraj {
    /// File type (`ENERTRJ`, `FRENERTRJ`) → its declaration.
    layouts: HashMap<String, Vec<Entry>>,
    /// Subblock name → `[rows][cols]` of the frame currently held.
    data: HashMap<String, Vec<Vec<f64>>>,
    /// `size` bindings of the frame currently held, including the `matrix_` forms.
    sizes: HashMap<String, usize>,
    constants: HashMap<String, f64>,
    variables: HashMap<String, Expr>,
    version: Option<String>,
}

impl EnergyTraj {
    /// Read a library file: the `ENERTRJ`/`FRENERTRJ` declarations, `VARIABLES` and
    /// `ENEVERSION`.
    pub fn from_library(text: &str) -> Result<Self, String> {
        let mut traj = EnergyTraj::default();
        // gromos++ `set_standards`
        traj.add_constant("BOLTZ", gromos_core::units::kB);

        for (name, body) in blocks(text) {
            match name.as_str() {
                "ENERTRJ" | "FRENERTRJ" => {
                    let mut entries = Vec::new();
                    for line in body {
                        traj.parse_declaration(&line, &mut entries)?;
                    }
                    traj.layouts.insert(name, entries);
                },
                "ENEVERSION" => {
                    let v: String = body.concat().split_whitespace().collect();
                    if traj.version.is_some() {
                        return Err("library has two ENEVERSION blocks".into());
                    }
                    traj.version = Some(v);
                },
                "VARIABLES" => {
                    for (var, expr) in variable_definitions(&body) {
                        let parsed = parse_expression(&expr)
                            .map_err(|e| format!("library variable '{var} = {expr}': {e}"))?;
                        traj.variables.insert(var, parsed);
                    }
                },
                _ => {},
            }
        }
        if traj.layouts.is_empty() {
            return Err("library declares neither ENERTRJ nor FRENERTRJ".into());
        }
        Ok(traj)
    }

    /// A named constant available to every expression (`ene_ana` adds `MASS` and `NUMMOL`
    /// from the topology).
    pub fn add_constant(&mut self, name: &str, value: f64) {
        self.constants.insert(name.to_string(), value);
    }

    pub fn version(&self) -> Option<&str> {
        self.version.as_deref()
    }

    /// The names a library declared, for `@library print`.
    pub fn variable_names(&self) -> Vec<&str> {
        let mut v: Vec<&str> = self.variables.keys().map(|s| s.as_str()).collect();
        v.sort_unstable();
        v
    }

    fn parse_declaration(&mut self, line: &str, entries: &mut Vec<Entry>) -> Result<(), String> {
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.is_empty() {
            return Ok(());
        }
        match t[0].to_ascii_lowercase().as_str() {
            "block" => {
                if t.len() != 2 {
                    return Err(format!("'block' takes one name: {line}"));
                }
                entries.push(Entry::Block(t[1].to_string()));
            },
            "size" => {
                if t.len() != 2 {
                    return Err(format!("'size' takes one name: {line}"));
                }
                // gromos++ registers the plain and the matrix form together
                self.sizes.entry(t[1].to_string()).or_insert(0);
                self.sizes.entry(format!("matrix_{}", t[1])).or_insert(0);
                entries.push(Entry::Size(t[1].to_string()));
            },
            "subblock" => {
                if t.len() != 4 {
                    return Err(format!("'subblock' takes a name and two sizes: {line}"));
                }
                let dim = |s: &str| -> Result<Dim, String> {
                    if let Ok(n) = s.parse::<usize>() {
                        return Ok(Dim::Fixed(n));
                    }
                    if self.sizes.contains_key(s) {
                        Ok(Dim::Named(s.to_string()))
                    } else {
                        Err(format!("variable name {s} in block size not defined"))
                    }
                };
                entries.push(Entry::Sub {
                    name: t[1].to_string(),
                    rows: dim(t[2])?,
                    cols: dim(t[3])?,
                });
                self.data.entry(t[1].to_string()).or_default();
            },
            other => return Err(format!("unknown library keyword '{other}': {line}")),
        }
        Ok(())
    }

    fn dim(&self, d: &Dim) -> usize {
        match d {
            Dim::Fixed(n) => *n,
            Dim::Named(n) => self.sizes.get(n).copied().unwrap_or(0),
        }
    }

    /// Read the next frame of `file_type` (`ENERTRJ`/`FRENERTRJ`); `Ok(false)` at end of file.
    pub fn read_frame(&mut self, reader: &mut TextReader, file_type: &str) -> Result<bool, String> {
        let layout = self
            .layouts
            .get(file_type)
            .ok_or_else(|| format!("don't know how to read from file type {file_type}"))?
            .clone();

        let mut values: std::vec::IntoIter<f64> = Vec::new().into_iter();
        let mut block = String::new();
        // gromos++ guards against a library that does not describe the file it is reading
        // (`EnergyTraj::read_block`); without it a misaligned `size` silently yields empty
        // subblocks and the mismatch only surfaces as a missing value much later.
        let mut check_leftover: Option<String> = None;
        for entry in &layout {
            match entry {
                Entry::Block(name) => {
                    let Some((got, body)) = read_gromos_block(reader)? else {
                        return Ok(false);
                    };
                    if &got != name {
                        return Err(format!("Block {name} expected. Got {got}"));
                    }
                    let mut nums = Vec::new();
                    for tok in body.iter().flat_map(|l| l.split_whitespace()) {
                        nums.push(
                            tok.parse::<f64>()
                                .map_err(|_| format!("{name}: '{tok}' is not a number"))?,
                        );
                    }
                    if let Some(leftover) = check_leftover.take() {
                        return Err(leftover);
                    }
                    values = nums.into_iter();
                    block = name.clone();
                },
                Entry::Size(name) => {
                    let x = values
                        .next()
                        .ok_or_else(|| format!("Tried to read an integer for {name}"))?
                        as usize;
                    self.sizes.insert(name.clone(), x);
                    self.sizes.insert(format!("matrix_{name}"), x * (x + 1) / 2);
                },
                Entry::Sub { name, rows, cols } => {
                    let (si, sj) = (self.dim(rows), self.dim(cols));
                    let mut table = vec![vec![0.0; sj]; si];
                    for row in table.iter_mut() {
                        for cell in row.iter_mut() {
                            *cell = values
                                .next()
                                .ok_or_else(|| format!("Not enough values in {name}"))?;
                        }
                    }
                    self.data.insert(name.clone(), table);
                },
            }
            if !block.is_empty() {
                let rest = values.len();
                check_leftover = (rest > 0).then(|| {
                    format!(
                        "Block definition does not agree with trajectory data for {block}: {rest} leftover values"
                    )
                });
            }
        }
        if let Some(leftover) = check_leftover {
            return Err(leftover);
        }
        Ok(true)
    }

    /// The value of a property: `SUBBLOCK[i]`, `SUBBLOCK[i][j]`, a library variable, or a
    /// constant. Indices are 1-based, as in the library.
    pub fn value(&self, prop: &str) -> Result<f64, String> {
        let expr = parse_expression(prop)
            .map_err(|e| format!("Trying to access an unknown variable {prop}: {e}"))?;
        self.eval(&expr, &mut Vec::new())
    }

    fn eval(&self, e: &Expr, stack: &mut Vec<String>) -> Result<f64, String> {
        match e {
            Expr::Const(v) => Ok(*v),
            Expr::Name { name, i, j } => {
                if let Some(table) = self.data.get(name) {
                    return table
                        .get(*i)
                        .and_then(|r| r.get(*j))
                        .copied()
                        .ok_or_else(|| {
                            format!("{name}[{}][{}] is outside the block read", i + 1, j + 1)
                        });
                }
                if let Some(v) = self.constants.get(name) {
                    return Ok(*v);
                }
                if let Some(v) = self.sizes.get(name) {
                    return Ok(*v as f64);
                }
                if let Some(sub) = self.variables.get(name) {
                    if stack.iter().any(|s| s == name) {
                        return Err(format!("variable {name} is defined in terms of itself"));
                    }
                    stack.push(name.clone());
                    let v = self.eval(sub, stack);
                    stack.pop();
                    return v;
                }
                Err(format!("Trying to access an unknown variable {name}"))
            },
            Expr::Fun(f, a) => {
                let v = self.eval(a, stack)?;
                Ok(match f {
                    Fun::Sin => v.sin(),
                    Fun::Cos => v.cos(),
                    Fun::Exp => v.exp(),
                    Fun::Log => v.ln(),
                })
            },
            Expr::Bin(op, a, b) => {
                let (a, b) = (self.eval(a, stack)?, self.eval(b, stack)?);
                Ok(match op {
                    Op::Add => a + b,
                    Op::Sub => a - b,
                    Op::Mul => a * b,
                    Op::Div => a / b,
                })
            },
        }
    }
}

/// Split a GROMOS text file into `(block name, body lines)`; comments and blank lines dropped.
fn blocks(text: &str) -> Vec<(String, Vec<String>)> {
    let mut out = Vec::new();
    let mut name: Option<String> = None;
    let mut body: Vec<String> = Vec::new();
    for line in text.lines() {
        let line = line.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        if line == "END" {
            if let Some(n) = name.take() {
                out.push((n, std::mem::take(&mut body)));
            }
            continue;
        }
        match &name {
            None => name = Some(line.to_string()),
            Some(_) => body.push(line.to_string()),
        }
    }
    out
}

/// Read the next GROMOS block from a stream: `Ok(None)` at end of file.
fn read_gromos_block(reader: &mut TextReader) -> Result<Option<(String, Vec<String>)>, String> {
    let mut buf = String::new();
    let name = loop {
        if reader.read_line(&mut buf).map_err(|e| e.to_string())? == 0 {
            return Ok(None);
        }
        let line = buf.split('#').next().unwrap_or("").trim();
        if !line.is_empty() {
            break line.to_string();
        }
    };
    let mut body = Vec::new();
    loop {
        if reader.read_line(&mut buf).map_err(|e| e.to_string())? == 0 {
            return Err(format!("no END in {name} block"));
        }
        let line = buf.split('#').next().unwrap_or("").trim();
        if line == "END" {
            break;
        }
        if !line.is_empty() {
            body.push(line.to_string());
        }
    }
    Ok(Some((name, body)))
}

/// `name = expression` pairs from a `VARIABLES` body. gromos++ reads the block as one token
/// stream and cuts each definition at the token before the next `=`.
fn variable_definitions(body: &[String]) -> Vec<(String, String)> {
    let tokens: Vec<&str> = body.iter().flat_map(|l| l.split_whitespace()).collect();
    let eq: Vec<usize> = (0..tokens.len()).filter(|&i| tokens[i] == "=").collect();
    let mut out = Vec::new();
    for (k, &i) in eq.iter().enumerate() {
        if i == 0 {
            continue;
        }
        // up to the name of the next definition, or the end
        let end = eq.get(k + 1).map_or(tokens.len(), |&next| next - 1);
        let rhs = tokens[i + 1..end].join(" ");
        if !rhs.is_empty() {
            out.push((tokens[i - 1].to_string(), rhs));
        }
    }
    out
}

#[derive(Debug, Clone, PartialEq)]
enum Tok {
    Op(Op),
    Fun(Fun),
    Open,
    Close,
    Atom(String),
}

fn tokenize(s: &str) -> Vec<Tok> {
    s.split_whitespace()
        .map(|t| match t {
            "+" => Tok::Op(Op::Add),
            "-" => Tok::Op(Op::Sub),
            "*" => Tok::Op(Op::Mul),
            "/" => Tok::Op(Op::Div),
            "(" => Tok::Open,
            ")" => Tok::Close,
            "sin" => Tok::Fun(Fun::Sin),
            "cos" => Tok::Fun(Fun::Cos),
            "exp" => Tok::Fun(Fun::Exp),
            "log" => Tok::Fun(Fun::Log),
            other => Tok::Atom(other.to_string()),
        })
        .collect()
}

fn parse_expression(s: &str) -> Result<Expr, String> {
    let toks = tokenize(s);
    if toks.is_empty() {
        return Err("empty expression".into());
    }
    parse_tokens(&toks)
}

/// gromos++ `Expression::calc`: split at the rightmost top-level `+`/`-`, else at the
/// *leftmost* top-level `*`/`/` — see the module note on associativity.
fn parse_tokens(toks: &[Tok]) -> Result<Expr, String> {
    if toks.is_empty() {
        return Err("empty expression".into());
    }
    let mut depth = 0i32;
    let mut add_sub: Option<usize> = None;
    let mut mul_div: Option<usize> = None;
    for (i, t) in toks.iter().enumerate() {
        match t {
            Tok::Open => depth += 1,
            Tok::Close => depth -= 1,
            Tok::Op(op) if depth == 0 && i > 0 => match op {
                Op::Add | Op::Sub => add_sub = Some(i), // keep the rightmost
                Op::Mul | Op::Div => {
                    if mul_div.is_none() {
                        mul_div = Some(i) // keep the leftmost
                    }
                },
            },
            _ => {},
        }
    }
    if depth != 0 {
        return Err("unbalanced brackets".into());
    }
    if let Some(k) = add_sub.or(mul_div) {
        let Tok::Op(op) = toks[k] else { unreachable!() };
        return Ok(Expr::Bin(
            op,
            Box::new(parse_tokens(&toks[..k])?),
            Box::new(parse_tokens(&toks[k + 1..])?),
        ));
    }
    parse_atom(toks)
}

fn parse_atom(toks: &[Tok]) -> Result<Expr, String> {
    match toks.first() {
        None => Err("empty expression".into()),
        Some(Tok::Op(_)) => Err("an expression cannot start with an operator".into()),
        Some(Tok::Fun(f)) => Ok(Expr::Fun(*f, Box::new(parse_atom(&toks[1..])?))),
        Some(Tok::Open) => {
            if toks.last() != Some(&Tok::Close) {
                return Err("unbalanced brackets".into());
            }
            parse_tokens(&toks[1..toks.len() - 1])
        },
        Some(Tok::Close) => Err("unbalanced brackets".into()),
        Some(Tok::Atom(a)) => {
            if toks.len() > 1 {
                return Err(format!(
                    "two values in a row ({a}) — an operator is missing"
                ));
            }
            if let Ok(v) = a.parse::<f64>() {
                return Ok(Expr::Const(v));
            }
            parse_reference(a)
        },
    }
}

/// `NAME`, `NAME[i]` or `NAME[i][j]`, with 1-based indices as the library writes them.
fn parse_reference(s: &str) -> Result<Expr, String> {
    let Some(open) = s.find('[') else {
        return Ok(Expr::Name {
            name: s.to_string(),
            i: 0,
            j: 0,
        });
    };
    let name = s[..open].to_string();
    let mut idx = Vec::new();
    for part in s[open..].split('[').skip(1) {
        let close = part
            .find(']')
            .ok_or_else(|| format!("missing ']' in {s}"))?;
        let n: usize = part[..close]
            .trim()
            .parse()
            .map_err(|_| format!("index in {s} is not an integer"))?;
        if n == 0 {
            return Err(format!("Indexes start at 1, invalid variable {s}"));
        }
        idx.push(n - 1);
    }
    match idx.len() {
        1 => Ok(Expr::Name {
            name,
            i: idx[0],
            j: 0,
        }),
        2 => Ok(Expr::Name {
            name,
            i: idx[0],
            j: idx[1],
        }),
        _ => Err(format!("{s} has more than two indices")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const LIB: &str = "\
TITLE
a library
END
ENEVERSION
  2023-04-15
END
ENERTRJ
  block TIMESTEP
    subblock TIME 2 1
  block ENERGY03
    subblock ENER 4 1
    size NUM_BATHS
    subblock KINENER NUM_BATHS 3
    size NUM_ENERGY_GROUPS
    subblock NONBONDED matrix_NUM_ENERGY_GROUPS 2
END
VARIABLES
time = TIME[2]
totene = ENER[1]
totpot = ENER[3]
sum = ENER[1] + ENER[3]
scaled = ENER[1] * 2 / 4
nested = totene + totpot
END
";

    /// A unique scratch path: the tests run in parallel and each deletes its own file.
    fn scratch(tag: &str) -> std::path::PathBuf {
        static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        std::env::temp_dir().join(format!(
            "gromos_etrj_{tag}_{}_{}.tre",
            std::process::id(),
            N.fetch_add(1, std::sync::atomic::Ordering::SeqCst)
        ))
    }

    fn frame() -> (EnergyTraj, TextReader) {
        // two baths (3 values each), two energy groups → 3 pairs of 2 values
        let text = "\
TIMESTEP
   1  0.002
END
ENERGY03
  10.0 20.0 30.0 40.0
  2
  1 2 3   4 5 6
  2
  7 8   9 10   11 12
END
";
        let path = scratch("frame");
        std::fs::write(&path, text).unwrap();
        let reader = TextReader::open(&path).unwrap();
        std::fs::remove_file(&path).ok();
        (EnergyTraj::from_library(LIB).unwrap(), reader)
    }

    #[test]
    fn library_declares_the_layout_and_the_names() {
        let traj = EnergyTraj::from_library(LIB).unwrap();
        assert_eq!(traj.version(), Some("2023-04-15"));
        assert_eq!(
            traj.variable_names(),
            ["nested", "scaled", "sum", "time", "totene", "totpot"]
        );
    }

    #[test]
    fn a_frame_fills_the_declared_subblocks() {
        let (mut traj, mut reader) = frame();
        assert!(traj.read_frame(&mut reader, "ENERTRJ").unwrap());

        assert_eq!(traj.value("TIME[2]").unwrap(), 0.002);
        assert_eq!(traj.value("totene").unwrap(), 10.0);
        assert_eq!(traj.value("totpot").unwrap(), 30.0);
        // the size read from the stream sets the following subblock's row count
        assert_eq!(traj.value("KINENER[2][3]").unwrap(), 6.0);
        // matrix_N expands 2 energy groups to 3 pairs
        assert_eq!(traj.value("NONBONDED[3][2]").unwrap(), 12.0);
        assert_eq!(traj.value("sum").unwrap(), 40.0);
        assert_eq!(traj.value("nested").unwrap(), 40.0);

        // one frame only
        assert!(!traj.read_frame(&mut reader, "ENERTRJ").unwrap());
    }

    #[test]
    fn indices_are_one_based_and_bounds_checked() {
        let (mut traj, mut reader) = frame();
        traj.read_frame(&mut reader, "ENERTRJ").unwrap();
        assert!(traj.value("ENER[0]").unwrap_err().contains("start at 1"));
        assert!(traj.value("ENER[5]").unwrap_err().contains("outside"));
        assert!(traj.value("NOSUCH").unwrap_err().contains("unknown"));
    }

    #[test]
    fn a_block_out_of_order_is_an_error() {
        let text = "ENERGY03\n 1 2 3 4\n 1\n 1 2 3\n 1\n 1 2\nEND\n";
        let path = scratch("order");
        std::fs::write(&path, text).unwrap();
        let mut reader = TextReader::open(&path).unwrap();
        std::fs::remove_file(&path).ok();
        let mut traj = EnergyTraj::from_library(LIB).unwrap();
        let err = traj.read_frame(&mut reader, "ENERTRJ").unwrap_err();
        assert!(err.contains("TIMESTEP expected"), "{err}");
    }

    /// gromos++ splits at the rightmost `+`/`-` (left-associative) but the leftmost `*`/`/`
    /// (right-associative). Both are reproduced.
    #[test]
    fn operator_associativity_matches_gromospp() {
        let ev = |s: &str| {
            EnergyTraj::default()
                .eval(&parse_expression(s).unwrap(), &mut Vec::new())
                .unwrap()
        };
        assert_eq!(ev("10 - 4 - 3"), 3.0); // (10-4)-3, not 10-(4-3)
        assert_eq!(ev("2 / 4 / 8"), 4.0); // 2/(4/8), the gromos++ grouping
        assert_eq!(ev("1 + 2 * 3"), 7.0); // + binds loosest
        assert_eq!(ev("( 1 + 2 ) * 3"), 9.0);
        assert_eq!(ev("exp ( 0 ) + 4"), 5.0);
    }

    #[test]
    fn a_malformed_expression_is_rejected() {
        assert!(parse_expression("+ 1").is_err());
        assert!(parse_expression("1 2").is_err());
        assert!(parse_expression("( 1 + 2").is_err());
    }
}
