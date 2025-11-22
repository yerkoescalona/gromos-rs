"""
NMR Restraints Example
======================

This example demonstrates the use of NMR restraints in GROMOS simulations:
- J-value restraints (Karplus relation for dihedral angles)
- RDC restraints (Residual Dipolar Couplings)

These restraints are essential for NMR structure refinement and validation.
"""

import gromos
import numpy as np


def jvalue_restraints_example():
    """
    J-value restraints using Karplus relation.

    The Karplus equation relates the scalar coupling constant J to the
    dihedral angle φ:
        J(φ) = A cos²(φ) + B cos(φ) + C

    Standard parameters for ³J(HN-Hα):
    - A = 6.4 Hz
    - B = -1.4 Hz
    - C = 1.9 Hz
    """
    print("=" * 70)
    print("J-Value Restraints (Karplus Relation)")
    print("=" * 70)

    # Create J-value restraints for backbone dihedral angles
    # Typical values for protein backbone ³J(HN-Hα) couplings

    # Restraint 1: φ angle in alpha-helix region (~-60°)
    # Expected J-coupling: ~4 Hz
    jval1 = gromos.JValueRestraint(
        atoms=[4, 6, 8, 14],  # N-CA-C-N atoms defining φ
        target_j=4.2,  # Experimental J-value (Hz)
        force_constant=10.0,  # Force constant (kJ/mol/Hz²)
        karplus_a=6.4,  # Karplus parameter A
        karplus_b=-1.4,  # Karplus parameter B
        karplus_c=1.9,  # Karplus parameter C
        karplus_delta=0.0,  # Phase shift
    )

    print(f"\nJ-value restraint 1: {jval1}")
    print(f"  Atoms: {jval1.atoms}")
    print(f"  Target J: {jval1.target_j} Hz")
    print(f"  Karplus parameters: A={jval1.karplus_a}, B={jval1.karplus_b}, C={jval1.karplus_c}")

    # Restraint 2: φ angle in beta-sheet region (~-120°)
    # Expected J-coupling: ~9 Hz
    jval2 = gromos.JValueRestraint(
        atoms=[24, 26, 28, 34],
        target_j=8.8,
        force_constant=10.0,
        karplus_a=6.4,
        karplus_b=-1.4,
        karplus_c=1.9,
    )

    print(f"\nJ-value restraint 2: {jval2}")
    print(f"  Target J: {jval2.target_j} Hz (beta-sheet region)")

    # Calculate J-value from dihedral angle
    phi_helix = -60.0 * np.pi / 180.0  # Convert to radians
    j_calculated = (
        jval1.karplus_a * np.cos(phi_helix) ** 2
        + jval1.karplus_b * np.cos(phi_helix)
        + jval1.karplus_c
    )

    print(f"\n  Predicted J for φ = -60°: {j_calculated:.2f} Hz")

    # Different Karplus parameters for different couplings
    # ³J(Hα-Hβ) has different parameters
    jval_side = gromos.JValueRestraint(
        atoms=[6, 8, 10, 12],
        target_j=6.5,
        force_constant=10.0,
        karplus_a=9.5,  # Different A for side chain
        karplus_b=-1.6,  # Different B
        karplus_c=1.8,  # Different C
    )

    print(f"\nSide chain J-value restraint: {jval_side}")
    print(f"  Different Karplus parameters for ³J(Hα-Hβ)")

    return [jval1, jval2, jval_side]


def rdc_restraints_example():
    """
    RDC (Residual Dipolar Coupling) restraints.

    RDCs provide information about the orientation of bond vectors
    relative to the alignment tensor. They are extremely valuable
    for determining the global fold of proteins.

    RDC depends on:
    - Maximum dipolar coupling D_max (bond type specific)
    - Alignment tensor (Saupe matrix)
    - Internuclear vector orientation
    """
    print("\n" + "=" * 70)
    print("RDC Restraints (Residual Dipolar Couplings)")
    print("=" * 70)

    # Create RDC restraints for N-H bonds
    # Typical D_max for N-H: ~24 kHz = -10,000 Hz in common convention

    # Saupe matrix components (5 independent: Sxx, Syy, Szz, Sxy, Sxz)
    # These describe the alignment tensor
    # Typically from SVD fit to experimental RDCs
    saupe_matrix = [
        0.00015,  # Sxx
        -0.00008,  # Syy
        -0.00007,  # Szz (Sxx + Syy + Szz = 0)
        0.00003,  # Sxy
        0.00001,  # Sxz
    ]

    # RDC 1: Helix N-H bond
    rdc1 = gromos.RDCRestraint(
        i=5,  # N atom
        j=6,  # H atom
        d_max=-10000.0,  # Maximum dipolar coupling for N-H (Hz)
        target_rdc=15.2,  # Experimental RDC (Hz)
        force_constant=5.0,  # Force constant (kJ/mol/Hz²)
        saupe_matrix=saupe_matrix,
        tau=0.0,  # No time averaging
        flat_bottom_width=0.0,  # No flat-bottom potential
        rdc_type=1,  # Type 1 = N-H bond
    )

    print(f"\nRDC restraint 1 (N-H): {rdc1}")
    print(f"  Atoms: {rdc1.i}-{rdc1.j}")
    print(f"  Target RDC: {rdc1.target_rdc} Hz")
    print(f"  D_max: {rdc1.d_max} Hz")

    # RDC 2: Different alignment
    rdc2 = gromos.RDCRestraint(
        i=25,
        j=26,
        d_max=-10000.0,
        target_rdc=-8.5,  # Negative RDC
        force_constant=5.0,
        saupe_matrix=saupe_matrix,
        tau=0.0,
        flat_bottom_width=2.0,  # Allow 2 Hz tolerance
        rdc_type=1,
    )

    print(f"\nRDC restraint 2 (N-H): {rdc2}")
    print(f"  Target RDC: {rdc2.target_rdc} Hz")
    print(f"  Flat-bottom width: {rdc2.flat_bottom_width} Hz")

    # RDC for C-H bonds (different D_max)
    # Typical D_max for C-H: ~60 kHz = -25,000 Hz
    rdc_ch = gromos.RDCRestraint(
        i=8,  # C atom
        j=9,  # H atom
        d_max=-25000.0,  # Larger for C-H
        target_rdc=12.3,
        force_constant=5.0,
        saupe_matrix=saupe_matrix,
        tau=10.0,  # Time averaging over 10 ps
        flat_bottom_width=0.0,
        rdc_type=0,  # Type 0 = C-H bond
    )

    print(f"\nRDC restraint (C-H): {rdc_ch}")
    print(f"  D_max: {rdc_ch.d_max} Hz (larger for C-H)")
    print(f"  Time averaging: τ = {rdc_ch.tau} ps")

    # Print Saupe matrix
    print(f"\nAlignment tensor (Saupe matrix):")
    print(f"  Sxx = {saupe_matrix[0]:.6f}")
    print(f"  Syy = {saupe_matrix[1]:.6f}")
    print(f"  Szz = {saupe_matrix[2]:.6f}")
    print(f"  Sxy = {saupe_matrix[3]:.6f}")
    print(f"  Sxz = {saupe_matrix[4]:.6f}")
    print(f"  Trace (should be ~0): {sum(saupe_matrix[:3]):.8f}")

    return [rdc1, rdc2, rdc_ch]


def combined_nmr_refinement():
    """
    Combined NMR refinement with both J-values and RDCs.

    This is a typical scenario for NMR structure refinement where
    both local (J-values) and global (RDCs) restraints are used.
    """
    print("\n" + "=" * 70)
    print("Combined NMR Refinement")
    print("=" * 70)

    # Collect all restraints
    j_restraints = jvalue_restraints_example()
    rdc_restraints = rdc_restraints_example()

    print(f"\nTotal restraints:")
    print(f"  J-value restraints: {len(j_restraints)}")
    print(f"  RDC restraints: {len(rdc_restraints)}")
    print(f"  Total: {len(j_restraints) + len(rdc_restraints)}")

    print("\nThese restraints would be added to the simulation:")
    print("  1. During equilibration: Gradually increase force constants")
    print("  2. During production: Full force constants for refinement")
    print("  3. Analysis: Calculate violations and NOE-like scores")

    # Example: Calculate restraint energy contribution
    print("\nTypical force constants:")
    print("  J-values: 5-20 kJ/mol/Hz²")
    print("  RDCs: 1-10 kJ/mol/Hz²")
    print("  (Start with lower, gradually increase)")


def advanced_nmr_features():
    """
    Advanced NMR restraint features.
    """
    print("\n" + "=" * 70)
    print("Advanced NMR Features")
    print("=" * 70)

    # Time-averaged restraints
    print("\n1. Time-Averaged Restraints:")
    print("   - Use τ parameter for averaging over timescale")
    print("   - Important for flexible regions")
    print("   - Typical τ: 1-100 ps")

    rdc_avg = gromos.RDCRestraint(
        i=15,
        j=16,
        d_max=-10000.0,
        target_rdc=7.5,
        force_constant=5.0,
        saupe_matrix=[0.00015, -0.00008, -0.00007, 0.00003, 0.00001],
        tau=50.0,  # 50 ps averaging
        flat_bottom_width=0.0,
        rdc_type=1,
    )
    print(f"   Example: {rdc_avg}")

    # Flat-bottom restraints
    print("\n2. Flat-Bottom Restraints:")
    print("   - Allow tolerance around experimental value")
    print("   - Useful for uncertain measurements")
    print("   - No penalty within flat-bottom width")

    rdc_flat = gromos.RDCRestraint(
        i=35,
        j=36,
        d_max=-10000.0,
        target_rdc=10.0,
        force_constant=5.0,
        saupe_matrix=[0.00015, -0.00008, -0.00007, 0.00003, 0.00001],
        tau=0.0,
        flat_bottom_width=3.0,  # ±3 Hz tolerance
        rdc_type=1,
    )
    print(f"   Example: {rdc_flat}")
    print(
        f"   Tolerance: {rdc_flat.target_rdc - rdc_flat.flat_bottom_width} to "
        f"{rdc_flat.target_rdc + rdc_flat.flat_bottom_width} Hz"
    )

    # Different bond types
    print("\n3. Different RDC Types:")
    print("   Type 0 (C-H):  D_max ≈ -25,000 Hz")
    print("   Type 1 (N-H):  D_max ≈ -10,000 Hz")
    print("   Type 2 (C-C):  D_max ≈ -5,000 Hz")
    print("   Type 3 (custom): User-defined D_max")

    # Multiple alignment media
    print("\n4. Multiple Alignment Media:")
    print("   - Use different Saupe matrices for different conditions")
    print("   - Improves structural accuracy")
    print("   - Reduces degeneracy in orientation determination")

    saupe_media1 = [0.00015, -0.00008, -0.00007, 0.00003, 0.00001]
    saupe_media2 = [-0.00010, 0.00012, -0.00002, -0.00005, 0.00004]

    print(f"   Media 1 Saupe: {saupe_media1}")
    print(f"   Media 2 Saupe: {saupe_media2}")


def main():
    """Run all NMR restraint examples."""
    print("\n" + "=" * 70)
    print("GROMOS NMR Restraints - Complete Example")
    print("=" * 70)

    # J-value restraints
    j_restraints = jvalue_restraints_example()

    # RDC restraints
    rdc_restraints = rdc_restraints_example()

    # Combined refinement
    combined_nmr_refinement()

    # Advanced features
    advanced_nmr_features()

    print("\n" + "=" * 70)
    print("Example complete!")
    print("=" * 70)
    print("\nKey takeaways:")
    print("  1. J-values constrain local dihedral angles via Karplus relation")
    print("  2. RDCs constrain global orientation via alignment tensor")
    print("  3. Time averaging helps with flexible regions")
    print("  4. Flat-bottom restraints account for experimental uncertainty")
    print("  5. Multiple alignment media improve structural accuracy")
    print("\nFor NMR structure refinement:")
    print("  - Start with loose restraints during equilibration")
    print("  - Gradually increase force constants")
    print("  - Monitor violations throughout simulation")
    print("  - Use ensemble of structures for final analysis")


if __name__ == "__main__":
    main()
