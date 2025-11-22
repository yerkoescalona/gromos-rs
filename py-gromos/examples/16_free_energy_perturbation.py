"""
Free Energy Perturbation (FEP) Example
=======================================

Demonstrates FEP calculations in GROMOS for alchemical transformations.
"""

import gromos


def basic_fep_setup():
    print("=" * 70)
    print("Free Energy Perturbation Setup")
    print("=" * 70)

    lambda_ctrl = gromos.LambdaController(lambda_value=0.0, dlambda_dt=0.0001, lambda_exponent=1.0)
    print(f"Lambda controller: {lambda_ctrl}")

    pert_atom = gromos.PerturbedAtom(
        atom_index=10,
        a_charge=0.0,
        b_charge=0.0,
        a_iac=12,
        b_iac=15,
        lj_softcore=0.5,
        crf_softcore=0.5,
    )
    print(f"Perturbed atom: {pert_atom}")
    return lambda_ctrl, pert_atom


def main():
    print("GROMOS Free Energy Perturbation Examples\n")
    lambda_ctrl, pert = basic_fep_setup()
    print("\nFEP simulation ready!")


if __name__ == "__main__":
    main()
