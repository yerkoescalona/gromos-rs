# Pre-commit Hooks Guide

This guide explains the pre-commit hook setup for gromosXX, inspired by Polars' developer tooling.

## Overview

Pre-commit hooks automatically check your code before you commit, catching issues early and ensuring consistent code quality across the project.

## Installation

### 1. Install pre-commit

```bash
pip install pre-commit
```

Or include it in your development environment:

```bash
cd py-gromos
make install-dev  # Includes pre-commit
```

### 2. Install Git Hooks

```bash
# From the repository root
pre-commit install
```

This installs the pre-commit hooks into your `.git/hooks/` directory. Now hooks will run automatically on `git commit`.

## What Gets Checked

Our pre-commit configuration includes:

### General Checks
- ✓ Trailing whitespace removal
- ✓ End-of-file fixer (ensures files end with newline)
- ✓ YAML/TOML/JSON syntax validation
- ✓ Large file detection (> 1MB)
- ✓ Merge conflict detection
- ✓ Case conflict detection
- ✓ Line ending normalization (LF)

### Rust Checks
- ✓ **rustfmt**: Code formatting
- ✓ **clippy**: Linting and code quality

### Python Checks
- ✓ **black**: Code formatting
- ✓ **ruff**: Linting, import sorting, code quality

### Documentation Checks
- ✓ **typos**: Spell checking
- ✓ **markdownlint**: Markdown linting

### Configuration Checks
- ✓ YAML formatting
- ✓ TOML validation

## Usage

### Automatic (Recommended)

Hooks run automatically on `git commit`:

```bash
git add .
git commit -m "feat: add new feature"
# Hooks run automatically and fix issues
```

If hooks make changes:
```bash
# Review the changes made by hooks
git diff

# Add the fixed files
git add .

# Commit again
git commit -m "feat: add new feature"
```

### Manual Run

Run hooks on all files without committing:

```bash
# Run all hooks
pre-commit run --all-files

# Run specific hook
pre-commit run black --all-files
pre-commit run cargo-fmt --all-files
pre-commit run typos --all-files
```

### Skip Hooks (Not Recommended)

In rare cases where you need to skip hooks:

```bash
git commit --no-verify -m "WIP: temporary commit"
```

**Note:** CI will still check your code, so skipping hooks locally just delays feedback.

## Hook Details

### rustfmt (Rust Formatting)

Automatically formats Rust code according to `rustfmt.toml`.

**Configuration:** `/rustfmt.toml`

**Key settings:**
- Line width: 100 characters
- Import grouping: StdExternalCrate
- Trailing commas in match blocks
- Field init shorthand

**Manual run:**
```bash
cargo fmt --all
```

### clippy (Rust Linting)

Catches common mistakes and suggests improvements.

**Manual run:**
```bash
cargo clippy --workspace --all-features --all-targets -- -D warnings
```

### black (Python Formatting)

Formats Python code with consistent style.

**Configuration:** `py-gromos/pyproject.toml`

**Key settings:**
- Line length: 100 characters
- Target Python 3.8+

**Manual run:**
```bash
cd py-gromos
black python/ tests/ examples/
```

### ruff (Python Linting)

Fast Python linter that replaces multiple tools.

**Configuration:** `py-gromos/pyproject.toml`

**Checks:**
- pycodestyle (E, W)
- pyflakes (F)
- isort (I)
- flake8-bugbear (B)
- flake8-comprehensions (C4)
- pyupgrade (UP)
- flake8-unused-arguments (ARG)
- flake8-simplify (SIM)

**Manual run:**
```bash
cd py-gromos
ruff check python/ tests/ examples/
ruff check --fix python/ tests/ examples/  # Auto-fix
```

### typos (Spell Checking)

Catches typos in code, comments, and documentation.

**Configuration:** `.typos.toml`

**Custom dictionary:**
- GROMOS-specific terms
- Scientific terminology
- Author names

**Manual run:**
```bash
typos
typos --write-changes  # Auto-fix
```

### markdownlint (Markdown Linting)

Ensures consistent Markdown formatting.

**Configuration:** `.markdownlint.yaml`

**Key settings:**
- Line length: 120 characters
- Allow inline HTML (for tables)

**Manual run:**
```bash
markdownlint '**/*.md' --fix
```

## Configuration Files

### `.pre-commit-config.yaml`
Main configuration file defining all hooks and their versions.

### `rustfmt.toml`
Rust formatting rules (inspired by Polars).

### `.typos.toml`
Spell checker dictionary and exclusions.

### `.editorconfig`
Editor settings for consistent formatting across IDEs.

### `.gitattributes`
Git attributes for line endings and diff behavior.

### `.markdownlint.yaml`
Markdown linting rules.

### `py-gromos/pyproject.toml`
Python tool configurations (black, ruff, pytest, mypy, coverage).

## Troubleshooting

### Hook is failing unexpectedly

```bash
# Update hooks to latest versions
pre-commit autoupdate

# Clear cache and retry
pre-commit clean
pre-commit run --all-files
```

### Hook takes too long

```bash
# Run specific file types only
pre-commit run --files *.rs
pre-commit run --files *.py
```

### Disable specific hook temporarily

Edit `.pre-commit-config.yaml` and comment out the hook:

```yaml
# - repo: https://github.com/pre-commit/pre-commit-hooks
#   rev: v4.5.0
#   hooks:
#     - id: trailing-whitespace
```

### False positive in typos

Add to `.typos.toml`:

```toml
[default.extend-words]
myword = "myword"
```

## CI Integration

Pre-commit hooks are also run in CI via GitHub Actions. If you forget to run them locally, CI will catch issues.

**CI Pre-commit workflow:**
1. Install pre-commit
2. Run all hooks
3. Fail if any hook fails

This ensures all code merged to main passes quality checks.

## Best Practices

### 1. Run Before Committing
Always let pre-commit run before committing. It's faster to fix issues locally than wait for CI.

### 2. Fix, Don't Skip
If a hook fails, fix the issue rather than skipping with `--no-verify`.

### 3. Update Regularly
Keep hooks updated:
```bash
pre-commit autoupdate
```

### 4. Understand Failures
Read hook output to understand what needs fixing:
```bash
pre-commit run --all-files --verbose
```

### 5. IDE Integration
Configure your IDE to use the same formatters:
- **VSCode**: Install Rust Analyzer, Python, Black Formatter extensions
- **PyCharm**: Enable black and ruff in settings
- **Vim/Neovim**: Use ALE or similar linting plugins

## Advanced Usage

### Custom Hooks

Add project-specific hooks to `.pre-commit-config.yaml`:

```yaml
- repo: local
  hooks:
    - id: my-custom-check
      name: My Custom Check
      entry: ./scripts/my-check.sh
      language: system
      pass_filenames: false
```

### Selective Commits

Run specific hooks only:

```bash
SKIP=cargo-clippy git commit -m "message"  # Skip clippy
SKIP=black,ruff git commit -m "message"    # Skip Python checks
```

### Pre-push Hooks

Install pre-push hooks to run tests before pushing:

```bash
pre-commit install --hook-type pre-push
```

Add to `.pre-commit-config.yaml`:

```yaml
- repo: local
  hooks:
    - id: tests
      name: run tests
      entry: make test
      language: system
      pass_filenames: false
      stages: [push]
```

## Summary

Pre-commit hooks help maintain code quality by:
- ✓ Catching issues before they reach CI
- ✓ Enforcing consistent style
- ✓ Reducing review iterations
- ✓ Making the codebase more maintainable

Install them once, and they'll help you write better code automatically!

## Additional Resources

- [Pre-commit Documentation](https://pre-commit.com/)
- [Contributing Guide](contributing.md)
- [CI/CD Guide](ci-cd.md)
- [Polars Project](https://github.com/pola-rs/polars) (Our inspiration)
