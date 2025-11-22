# Naming Conventions for Scientific Code

This project uses flexible naming conventions that respect established scientific notation rather than strictly enforcing programming language conventions.

## Philosophy

In scientific computing, certain variable names have universally recognized meanings. Forcing these into snake_case can make code **less** readable to domain experts. For example:

- `kT` (Boltzmann constant × Temperature) is clearer than `k_t`
- `pH` (hydrogen ion concentration) is clearer than `p_h`
- `deltaG` (Gibbs free energy change) is clearer than `delta_g`

## Allowed Naming Patterns

### Rust

Clippy warnings for naming conventions are relaxed:
- ✓ `kT` instead of `k_t`
- ✓ `pH` instead of `p_h`
- ✓ `pKa` instead of `p_ka`
- ✓ `deltaG` instead of `delta_g`
- ✓ Acronyms in any case: `MD`, `PME`, `REMD`

**Configuration:** `clippy.toml`

```rust
// All of these are acceptable:
let kT = 2.479;           // Boltzmann constant * Temperature
let pH = 7.4;             // Hydrogen ion concentration
let deltaG = -50.0;       // Free energy change
let rmsd = 0.15;          // Root mean square deviation
let REMD_steps = 1000;    // Replica Exchange MD steps
```

### Python

Ruff naming checks (N series) are disabled for scientific variables:
- ✓ Variables: `kT`, `pH`, `deltaG`
- ✓ Function arguments: `def calculate_energy(kT, pH)`
- ✓ Class attributes: `self.kT = kT`

**Configuration:** `py-gromos/pyproject.toml`

```python
# All of these are acceptable:
kT = 2.479              # Boltzmann constant * Temperature
pH = 7.4                # Hydrogen ion concentration
deltaG = -50.0          # Free energy change

def calculate_boltzmann_factor(kT: float, energy: float) -> float:
    """Calculate exp(-E/kT)."""
    return np.exp(-energy / kT)

class ThermodynamicState:
    def __init__(self, kT: float, pH: float):
        self.kT = kT
        self.pH = pH
```

## When to Use Snake_Case

Still use snake_case for:
- Non-scientific variable names: `atom_count`, `file_path`
- Function names: `calculate_energy()`, `read_topology()`
- Module names: `interaction.py`, `topology.py`

## Common Scientific Variables

Here are common scientific variables that should **not** be converted to snake_case:

### Thermodynamic Variables
- `kT` - Thermal energy (Boltzmann constant × Temperature)
- `T` - Temperature
- `P` - Pressure
- `V` - Volume
- `E` - Energy
- `H` - Enthalpy
- `S` - Entropy
- `G` - Gibbs free energy
- `deltaG`, `deltaH`, `deltaS` - Changes in thermodynamic quantities

### Physical Chemistry
- `pH` - Hydrogen ion concentration
- `pKa` - Acid dissociation constant
- `Kd` - Dissociation constant
- `Ka` - Association constant
- `Kb` - Base dissociation constant

### Molecular Dynamics
- `MD` - Molecular Dynamics
- `MC` - Monte Carlo
- `REMD` - Replica Exchange MD
- `PME` - Particle Mesh Ewald
- `LJ` - Lennard-Jones
- `rmsd` - Root Mean Square Deviation
- `rmsf` - Root Mean Square Fluctuation

### Mathematical/Physical Constants
- `kB` - Boltzmann constant
- `NA` - Avogadro's number
- `hbar` - Reduced Planck constant

### Computational
- `MPI` - Message Passing Interface
- `CUDA` - Compute Unified Device Architecture
- `SIMD` - Single Instruction Multiple Data
- `FFT` - Fast Fourier Transform
- `FFTW` - Fastest Fourier Transform in the West

## Examples in Context

### Good (Scientific Notation)

```rust
// Rust
pub fn boltzmann_factor(energy: f64, kT: f64) -> f64 {
    (-energy / kT).exp()
}

pub struct MDConfig {
    pub temperature: f64,  // K
    pub kT: f64,           // kJ/mol (computed)
    pub pH: f64,           // Unitless
}
```

```python
# Python
def calculate_pKa(Ka: float) -> float:
    """Calculate pKa from Ka."""
    return -np.log10(Ka)

class Thermostat:
    def __init__(self, temperature: float, kT: float):
        self.temperature = temperature
        self.kT = kT
```

### Avoid (Over-enforcement of snake_case)

```rust
// Don't do this - less readable to domain experts
pub fn boltzmann_factor(energy: f64, k_t: f64) -> f64 {  // ❌
    (-energy / k_t).exp()
}

pub struct MDConfig {
    pub temperature: f64,
    pub k_t: f64,      // ❌ Not standard notation
    pub p_h: f64,      // ❌ Not standard notation
}
```

## Linter Configuration

### Rust (clippy.toml)
```toml
# Allow mixed-case names for scientific variables
allow-mixed-case-names = true
```

### Python (pyproject.toml)
```toml
[tool.ruff.lint]
ignore = [
    "N801",  # class names
    "N802",  # function names
    "N803",  # argument names (allow kT, pH, deltaG)
    "N806",  # variable names in functions
    "N815",  # mixed case variable in class scope
    "N816",  # mixed case variable in global scope
]
```

### Pre-commit Hooks
Clippy is configured to allow scientific naming:
```bash
cargo clippy -- -D warnings -A clippy::upper_case_acronyms -A clippy::non_snake_case
```

## Style Guide Summary

| Type | Convention | Example |
|------|-----------|---------|
| Scientific variables | Domain standard | `kT`, `pH`, `deltaG` |
| Regular variables | snake_case | `atom_count`, `step_size` |
| Functions | snake_case | `calculate_energy()` |
| Structs/Classes | PascalCase | `SystemState`, `Topology` |
| Modules | snake_case | `topology.py`, `interaction.rs` |
| Constants | SCREAMING_SNAKE | `AVOGADRO_NUMBER` |
| Acronyms | Flexible | `MD`, `PME`, `REMD` acceptable |

## References

This approach is common in scientific computing projects:
- **OpenMM**: Uses `kT` throughout
- **GROMACS**: Uses scientific notation in code and papers
- **MDAnalysis**: Allows flexible naming for scientific variables
- **ASE**: Uses `kB`, `eV`, `Angstrom`

## Questions?

See [Contributing Guide](contributing.md) for general guidelines or open an issue if you're unsure about a specific naming case.
