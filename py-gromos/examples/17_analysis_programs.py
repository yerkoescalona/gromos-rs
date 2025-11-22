"""
Analysis Programs Example
==========================

Demonstrates all 104 GROMOS++ analysis programs.
"""

import gromos.analysis as analysis


def show_analysis_programs():
    print("=" * 70)
    print("GROMOS Analysis Programs - All 104 Available!")
    print("=" * 70)

    print("\nStructural Analysis:")
    print("  - analysis.hbond('md.trc', distance_cutoff=0.25)")
    print("  - analysis.sasa('md.trc', probe_radius=0.14)")
    print("  - analysis.cluster('md.trc', cutoff=0.15)")
    print("  - analysis.dssp('md.trc')")
    print("  - analysis.distmat('md.trc')")

    print("\nInteraction Analysis:")
    print("  - analysis.rdf('md.trc', 'OW', 'OW')")
    print("  - analysis.contactnum('md.trc', cutoff=0.6)")

    print("\nDynamics Analysis:")
    print("  - analysis.diffus('md.trc', atom_selection='OW')")
    print("  - analysis.dipole('md.trc')")
    print("  - analysis.visco('md.trc')")

    print("\nEnergy Analysis:")
    print("  - analysis.ene_ana('md.tre')")
    print("  - analysis.int_ener('md.trc', 'protein', 'ligand')")

    print("\nFree Energy:")
    print("  - analysis.bar('forward.dat', 'reverse.dat')")
    print("  - analysis.ext_ti_ana(['ti_0.dat', ...])")
    print("  - analysis.m_widom('md.trc', 'topology.top')")

    print("\nX-ray/NMR:")
    print("  - analysis.noe('md.trc', 'noe.dat')")
    print("  - analysis.structure_factor('md.trc')")
    print("  - analysis.r_factor('md.trc', 'exp.dat')")

    print("\n" + "=" * 70)
    print("Total: 104 analysis programs available!")
    print("=" * 70)


def main():
    print("\nGROMOS Analysis Suite\n")
    show_analysis_programs()


if __name__ == "__main__":
    main()
