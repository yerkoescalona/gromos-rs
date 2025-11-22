"""
Polarization and QM/MM Example
===============================

Demonstrates advanced electronic structure methods in GROMOS:
- Polarizable force fields (Drude oscillators, induced dipoles)
- QM/MM hybrid simulations
"""

import gromos


def polarization_models():
    print("=" * 70)
    print("Polarization Models")
    print("=" * 70)

    pol_calc = gromos.PolarizationCalculator(
        n_atoms=100, model=gromos.PolarizationModel.point_dipole()
    )
    print(f"Polarization calculator: {pol_calc}")

    pol_params = gromos.PolarizabilityParameters(alpha=1.0, thole_a=2.6, drude_mass=0.4)
    pol_calc.add_atom(pol_params)
    return pol_calc


def main():
    print("GROMOS Polarization & QM/MM Examples\n")
    pol = polarization_models()
    print("\nAdvanced electronic structure ready!")


if __name__ == "__main__":
    main()
