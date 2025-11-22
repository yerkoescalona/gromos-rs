# Developer Tools Summary

Quick reference for all developer tools in gromosXX, inspired by [Polars](https://github.com/pola-rs/polars).

## Quick Start

```bash
# Setup everything
./scripts/setup-dev.sh

# Install pre-commit hooks
pip install pre-commit
pre-commit install

# Run quality checks
make pre-commit
```

## Configuration Files

| File | Purpose | Inspired By |
|------|---------|-------------|
| `.pre-commit-config.yaml` | Pre-commit hook configuration | Polars workflow |
| `rustfmt.toml` | Rust code formatting rules | Polars `rustfmt.toml` |
| `.typos.toml` | Spell checker configuration | Polars `.typos.toml` |
| `.editorconfig` | Editor settings (cross-IDE) | Standard practice |
| `.gitattributes` | Git file handling | Standard practice |
| `.markdownlint.yaml` | Markdown linting rules | Common practice |
| `py-gromos/pyproject.toml` | Python tool configurations | Polars `pyproject.toml` |

## Tool Overview

### Rust Tools

#### rustfmt
**Purpose:** Code formatting
**Config:** `rustfmt.toml`
**Run:** `cargo fmt --all`

**Key Settings:**
- Line width: 100
- Import grouping: StdExternalCrate
- Field init shorthand: enabled

#### clippy
**Purpose:** Linting and code quality
**Run:** `cargo clippy --workspace --all-features -- -D warnings`

**What it checks:**
- Common mistakes
- Performance issues
- Idiomatic patterns
- Type safety

### Python Tools

#### black
**Purpose:** Code formatting
**Config:** `py-gromos/pyproject.toml`
**Run:** `black python/ tests/ examples/`

**Key Settings:**
- Line length: 100
- Target: Python 3.8+

#### ruff
**Purpose:** Fast linting (replaces flake8, isort, etc.)
**Config:** `py-gromos/pyproject.toml`
**Run:** `ruff check --fix python/ tests/`

**Checks:**
- E/W: pycodestyle
- F: pyflakes
- I: isort
- B: bugbear
- C4: comprehensions
- UP: pyupgrade
- ARG: unused arguments
- SIM: simplifications

#### mypy
**Purpose:** Static type checking
**Config:** `py-gromos/pyproject.toml`
**Run:** `mypy python/`

**Settings:**
- Check untyped defs
- Warn on redundant casts
- Strict equality

#### pytest
**Purpose:** Testing
**Config:** `py-gromos/pyproject.toml`
**Run:** `pytest tests/ -v`

**Features:**
- Markers: slow, benchmark, integration
- Coverage integration
- Strict warnings

### Cross-Language Tools

#### typos
**Purpose:** Spell checking
**Config:** `.typos.toml`
**Run:** `typos`

**Features:**
- Custom dictionary for GROMOS terms
- Scientific terminology support
- Exclude patterns

#### markdownlint
**Purpose:** Markdown linting
**Config:** `.markdownlint.yaml`
**Run:** `markdownlint '**/*.md'`

**Settings:**
- Line length: 120
- Allow inline HTML
- Code block style: fenced

#### pre-commit
**Purpose:** Git hook automation
**Config:** `.pre-commit-config.yaml`
**Run:** `pre-commit run --all-files`

**Hooks:**
- File checks (trailing whitespace, EOF, etc.)
- Rust: rustfmt, clippy
- Python: black, ruff
- Documentation: typos, markdownlint
- Config: YAML, TOML validation

## Comparison with Polars

| Feature | Polars | gromosXX | Status |
|---------|--------|----------|--------|
| rustfmt | ✓ | ✓ | Identical style |
| clippy | ✓ | ✓ | Same strictness |
| black | ✓ | ✓ | Same config |
| ruff | ✓ | ✓ | Similar rules |
| typos | ✓ | ✓ | Custom dictionary |
| dprint | ✓ | - | Using other tools |
| pre-commit | - | ✓ | Added feature |
| mypy | ✓ | ✓ | Similar config |
| pytest | ✓ | ✓ | Similar config |

## Makefile Commands

All tools are integrated into the Makefile:

### Code Quality
```bash
make fmt          # Format Rust and Python
make fmt-check    # Check formatting without changes
make lint         # Run all linters
make check        # Run cargo check
```

### Testing
```bash
make test         # Run all tests
make test-coverage # Run with coverage
make test-rust    # Rust tests only
make test-python  # Python tests only
```

### Complete Workflow
```bash
make pre-commit   # Format + Lint + Test
make ci           # Complete CI pipeline
make all          # Everything
```

## IDE Integration

### VSCode

**Extensions:**
- rust-analyzer
- Python
- Black Formatter
- Ruff
- EditorConfig

**Settings:**
```json
{
  "editor.formatOnSave": true,
  "rust-analyzer.check.command": "clippy",
  "[python]": {
    "editor.defaultFormatter": "ms-python.black-formatter"
  },
  "python.linting.ruffEnabled": true
}
```

### PyCharm

**Settings:**
- Enable black formatter
- Enable ruff linter
- Use rustfmt on save
- Enable EditorConfig support

### Neovim/Vim

**Plugins:**
- ALE or nvim-lspconfig
- rust.vim or rust-tools.nvim
- python-mode or vim-python-pep8-indent

## CI/CD Integration

All tools run in GitHub Actions:

**test-rust.yml:**
- rustfmt check
- clippy
- cargo test

**test-python.yml:**
- black check
- ruff check
- pytest with coverage

**ci.yml:**
- Complete integration
- All quality checks
- Coverage reporting

## Common Workflows

### Before Committing
```bash
make fmt          # Format code
make lint         # Check for issues
make test         # Run tests
# or
make pre-commit   # All of the above
```

### Daily Development
```bash
# Auto-format on save (IDE)
# Pre-commit hooks run automatically
git commit -m "feat: new feature"
```

### Before PR
```bash
make ci           # Full CI pipeline locally
./scripts/run-ci-locally.sh  # Comprehensive check
```

### Fixing Issues
```bash
# Auto-fix formatting
cargo fmt --all
black python/ tests/ examples/

# Auto-fix some lint issues
cargo clippy --fix --allow-dirty
ruff check --fix python/

# Check spelling
typos --write-changes
```

## Tool Update Strategy

**Weekly (automated):**
- Dependency updates via dependabot
- Pre-commit hook updates

**Monthly (manual review):**
- Tool configuration refinement
- New rule adoption

**Per Release:**
- Tool version pinning
- Configuration documentation

## Resources

### Documentation
- [CI/CD Guide](ci-cd.md) - Complete CI/CD guide
- [Pre-commit Guide](pre-commit.md) - Pre-commit setup
- [Contributing Guide](contributing.md) - Contribution guidelines

### External Resources
- [Polars Repository](https://github.com/pola-rs/polars) - Our inspiration
- [pre-commit.com](https://pre-commit.com/) - Pre-commit documentation
- [Ruff Documentation](https://docs.astral.sh/ruff/) - Ruff guide
- [Clippy Lints](https://rust-lang.github.io/rust-clippy/master/) - All clippy lints

### Tool Repositories
- [rustfmt](https://github.com/rust-lang/rustfmt)
- [clippy](https://github.com/rust-lang/rust-clippy)
- [black](https://github.com/psf/black)
- [ruff](https://github.com/astral-sh/ruff)
- [typos](https://github.com/crate-ci/typos)
- [pre-commit](https://github.com/pre-commit/pre-commit)

## Quick Reference Card

```
┌─────────────────────────────────────────────────┐
│ gromosXX Developer Tools Quick Reference        │
├─────────────────────────────────────────────────┤
│ Format:     make fmt                            │
│ Lint:       make lint                           │
│ Test:       make test                           │
│ All:        make pre-commit                     │
│                                                 │
│ Rust fmt:   cargo fmt --all                     │
│ Rust lint:  cargo clippy --workspace           │
│ Python fmt: black python/ tests/               │
│ Python lint: ruff check python/                │
│ Spell:      typos                               │
│                                                 │
│ Pre-commit: pre-commit run --all-files         │
│ CI local:   ./scripts/run-ci-locally.sh        │
└─────────────────────────────────────────────────┘
```

## Summary

gromosXX uses industry-standard developer tools inspired by Polars to ensure:

✓ **Consistent code style** across Rust and Python
✓ **Early error detection** via pre-commit hooks
✓ **Automated quality checks** in CI/CD
✓ **Easy onboarding** for new contributors
✓ **Maintainable codebase** over time

Install once with `./scripts/setup-dev.sh` and enjoy automated code quality!
