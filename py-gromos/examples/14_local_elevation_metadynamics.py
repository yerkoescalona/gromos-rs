"""
Local Elevation and Metadynamics Example
=========================================

This example demonstrates advanced sampling with local elevation (LE) and
metadynamics in GROMOS:

- Reaction coordinate definition
- Gaussian hill deposition
- Multi-dimensional umbrellas
- Free energy surface reconstruction
- Adaptive biasing force methods

Local elevation fills energy basins with Gaussian hills to encourage
exploration of rare events and calculation of free energy surfaces.
"""

import gromos
import numpy as np
import matplotlib.pyplot as plt


def simple_1d_metadynamics():
    """
    1D metadynamics along a dihedral angle.

    This is the simplest case: biasing along a single reaction coordinate
    (a dihedral angle) to explore conformational space.
    """
    print("=" * 70)
    print("1D Metadynamics - Dihedral Angle")
    print("=" * 70)

    # Define reaction coordinate: φ dihedral angle
    # Typical for protein backbone sampling
    coord = gromos.LECoordinate(
        umbrella_id=1,
        coord_type=gromos.CoordinateType.dihedral(),
        atoms=[4, 6, 8, 14],  # N-CA-C-N atoms
        reference_positions=None,  # No reference needed for dihedral
    )

    print(f"\nReaction coordinate: {coord}")
    print(f"  Type: Dihedral angle")
    print(f"  Atoms: {coord.atoms}")

    # Create umbrella for metadynamics
    umbrella = gromos.Umbrella(
        umbrella_id=1,
        coordinates=[coord],
        grid_sizes=[360],  # 360 bins for 360 degrees
        grid_mins=[-180.0],  # -180 degrees
        grid_maxs=[180.0],  # +180 degrees
        grid_spacings=[1.0],  # 1 degree per bin
        hill_height=0.5,  # 0.5 kJ/mol per hill
        gaussian_widths=[5.0],  # 5 degree width (σ)
        deposition_frequency=100,  # Add hill every 100 steps
        enabled=True,
        building=True,  # Build potential (not frozen)
        periodic=[True],  # Dihedral is periodic
    )

    print(f"\nUmbrella configuration: {umbrella}")
    print(f"  Dimensions: {umbrella.dimensionality}")
    print(f"  Grid size: {umbrella.grid_sizes}")
    print(f"  Hill height: {umbrella.hill_height} kJ/mol")
    print(f"  Gaussian σ: {umbrella.gaussian_widths[0]}°")
    print(f"  Deposition: Every {umbrella.deposition_frequency} steps")

    # Simulation parameters
    print("\nSimulation setup:")
    print("  Timestep: 2 fs")
    print("  Total time: 100 ns")
    print("  Hills deposited: ~500,000")
    print("  Expected barrier crossings: Multiple")

    print("\nFree energy surface:")
    print("  After simulation, integrate hill potential")
    print("  F(φ) = -V_bias(φ)")
    print("  Expect minima at α-helix (-60°) and β-sheet (-120°)")

    return umbrella


def distance_metadynamics():
    """
    Metadynamics along a distance coordinate.

    Useful for:
    - Protein-ligand unbinding
    - Ion permeation through channels
    - Conformational transitions
    """
    print("\n" + "=" * 70)
    print("1D Metadynamics - Distance Coordinate")
    print("=" * 70)

    # Distance between centers of mass of two groups
    # E.g., ligand COM to protein binding site COM
    coord = gromos.LECoordinate(
        umbrella_id=2,
        coord_type=gromos.CoordinateType.distance(),
        atoms=[10, 25],  # Atom indices for distance
        reference_positions=None,
    )

    print(f"\nReaction coordinate: {coord}")
    print(f"  Type: Distance")
    print(f"  Between atoms: {coord.atoms}")

    # Umbrella for unbinding pathway
    umbrella = gromos.Umbrella(
        umbrella_id=2,
        coordinates=[coord],
        grid_sizes=[200],  # 200 bins
        grid_mins=[0.0],  # 0 nm (bound)
        grid_maxs=[2.0],  # 2 nm (unbound)
        grid_spacings=[0.01],  # 0.01 nm spacing
        hill_height=0.3,  # Smaller hills for smoother FES
        gaussian_widths=[0.05],  # 0.05 nm = 0.5 Å width
        deposition_frequency=500,  # Less frequent for distance
        enabled=True,
        building=True,
        periodic=[False],  # Distance is NOT periodic
    )

    print(f"\nUmbrella: {umbrella}")
    print(f"  Range: {umbrella.grid_mins[0]:.1f} - {umbrella.grid_maxs[0]:.1f} nm")
    print(f"  Resolution: {umbrella.grid_spacings[0]:.3f} nm")

    print("\nApplication: Protein-ligand unbinding")
    print("  Bound state: ~0.0-0.5 nm")
    print("  Transition: ~0.5-1.2 nm")
    print("  Unbound: >1.5 nm")
    print("  Binding free energy: ΔG = F(unbound) - F(bound)")

    return umbrella


def metadynamics_2d():
    """
    2D metadynamics with two reaction coordinates.

    Allows exploration of complex free energy landscapes with
    multiple transition pathways.
    """
    print("\n" + "=" * 70)
    print("2D Metadynamics - φ/ψ Ramachandran")
    print("=" * 70)

    # Coordinate 1: φ angle
    coord_phi = gromos.LECoordinate(
        umbrella_id=3,
        coord_type=gromos.CoordinateType.dihedral(),
        atoms=[4, 6, 8, 14],  # φ dihedral
        reference_positions=None,
    )

    # Coordinate 2: ψ angle
    coord_psi = gromos.LECoordinate(
        umbrella_id=3,
        coord_type=gromos.CoordinateType.dihedral(),
        atoms=[6, 8, 14, 16],  # ψ dihedral
        reference_positions=None,
    )

    print(f"\nReaction coordinates:")
    print(f"  1. φ (phi): {coord_phi.atoms}")
    print(f"  2. ψ (psi): {coord_psi.atoms}")

    # 2D umbrella
    umbrella_2d = gromos.Umbrella(
        umbrella_id=3,
        coordinates=[coord_phi, coord_psi],
        grid_sizes=[360, 360],  # 360x360 grid
        grid_mins=[-180.0, -180.0],  # Both -180°
        grid_maxs=[180.0, 180.0],  # Both +180°
        grid_spacings=[1.0, 1.0],  # 1° spacing
        hill_height=0.2,  # Smaller for 2D
        gaussian_widths=[5.0, 5.0],  # 5° width in each dimension
        deposition_frequency=200,  # More frequent for 2D
        enabled=True,
        building=True,
        periodic=[True, True],  # Both periodic
    )

    print(f"\nUmbrella: {umbrella_2d}")
    print(f"  Dimensions: {umbrella_2d.dimensionality}D")
    print(f"  Grid: {umbrella_2d.grid_sizes[0]} × {umbrella_2d.grid_sizes[1]}")
    print(f"  Total bins: {umbrella_2d.grid_sizes[0] * umbrella_2d.grid_sizes[1]:,}")

    print("\nRamachandran plot exploration:")
    print("  α-helix region: φ ≈ -60°, ψ ≈ -45°")
    print("  β-sheet region: φ ≈ -120°, ψ ≈ +120°")
    print("  Left-handed helix: φ ≈ +60°, ψ ≈ +45°")

    print("\nComputational cost:")
    grid_points = umbrella_2d.grid_sizes[0] * umbrella_2d.grid_sizes[1]
    print(f"  Grid points: {grid_points:,}")
    print(f"  Memory: ~{grid_points * 8 / 1024 / 1024:.1f} MB (double precision)")
    print(f"  Convergence: Requires longer simulation than 1D")

    return umbrella_2d


def rmsd_metadynamics():
    """
    Metadynamics along RMSD coordinate.

    Useful for:
    - Exploring protein conformational changes
    - Sampling different folded states
    - Transition path sampling
    """
    print("\n" + "=" * 70)
    print("Metadynamics - RMSD Coordinate")
    print("=" * 70)

    # RMSD requires reference positions
    # E.g., crystal structure or native state
    ref_positions = [
        [0.0, 0.0, 0.0],  # Atom 0 position in reference
        [0.15, 0.0, 0.0],  # Atom 1
        [0.15, 0.15, 0.0],  # Atom 2
        [0.0, 0.15, 0.0],  # Atom 3
        # ... etc for all atoms
    ]

    coord_rmsd = gromos.LECoordinate(
        umbrella_id=4,
        coord_type=gromos.CoordinateType.rmsd(),
        atoms=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],  # Atoms for RMSD
        reference_positions=ref_positions[:10],  # Reference structure
    )

    print(f"\nReaction coordinate: {coord_rmsd}")
    print(f"  Type: RMSD")
    print(f"  Atoms: {len(coord_rmsd.atoms)}")

    # Umbrella for conformational sampling
    umbrella = gromos.Umbrella(
        umbrella_id=4,
        coordinates=[coord_rmsd],
        grid_sizes=[100],  # 100 bins
        grid_mins=[0.0],  # 0 nm (native)
        grid_maxs=[1.0],  # 1 nm (unfolded)
        grid_spacings=[0.01],  # 0.01 nm = 0.1 Å
        hill_height=0.5,
        gaussian_widths=[0.02],  # 0.02 nm = 0.2 Å
        deposition_frequency=100,
        enabled=True,
        building=True,
        periodic=[False],  # RMSD not periodic
    )

    print(f"\nUmbrella: {umbrella}")
    print(f"  Range: {umbrella.grid_mins[0]:.2f} - {umbrella.grid_maxs[0]:.2f} nm")

    print("\nApplications:")
    print("  - Protein folding/unfolding pathways")
    print("  - Conformational transitions (open ↔ closed)")
    print("  - Domain movements")
    print("  - Oligomerization interfaces")

    return umbrella


def advanced_metadynamics_features():
    """
    Advanced metadynamics features and variants.
    """
    print("\n" + "=" * 70)
    print("Advanced Metadynamics Features")
    print("=" * 70)

    print("\n1. Well-Tempered Metadynamics:")
    print("   - Hill height decreases over time")
    print("   - W(t) = W₀ * exp(-V(s,t) / (k_B * ΔT))")
    print("   - ΔT: bias temperature")
    print("   - Improves convergence")
    print("   - Self-limiting bias")

    print("\n2. Multiple Walkers:")
    print("   - Run multiple replicas simultaneously")
    print("   - Share bias potential in real-time")
    print("   - Parallel exploration of FES")
    print("   - Faster convergence")

    print("\n3. Bias Exchange:")
    print("   - Multiple replicas with different bias")
    print("   - Periodic exchanges between replicas")
    print("   - Improved sampling efficiency")

    print("\n4. Funnel Metadynamics:")
    print("   - External restraint confines sampling")
    print("   - Gradually widened funnel")
    print("   - Prevents unphysical conformations")

    print("\n5. Reconnaissance Metadynamics:")
    print("   - Short exploratory runs")
    print("   - Identify important regions")
    print("   - Focused longer simulations")


def free_energy_reconstruction():
    """
    Reconstructing free energy surfaces from metadynamics.
    """
    print("\n" + "=" * 70)
    print("Free Energy Surface Reconstruction")
    print("=" * 70)

    print("\n1. From Bias Potential:")
    print("   F(s) = -V_bias(s) + C")
    print("   Where:")
    print("     F(s): Free energy")
    print("     V_bias(s): Sum of all deposited hills")
    print("     C: Constant (arbitrary zero)")

    print("\n2. Hill Summation:")
    print("   V_bias(s,t) = Σ W_i * exp(-|s - s_i|² / (2σ²))")
    print("   Where:")
    print("     W_i: Height of hill i")
    print("     s_i: Position of hill i")
    print("     σ: Gaussian width")

    print("\n3. Error Estimation:")
    print("   - Block averaging")
    print("   - Multiple independent runs")
    print("   - Well-tempered convergence check")
    print("   - Free energy change monitoring")

    print("\n4. Reweighting:")
    print("   - Can calculate properties at any s")
    print("   - Ensemble averages: <A> = ∫ A(s) e^(-F(s)/kT) ds")
    print("   - Unbiased distribution: ρ(s) = e^(-F(s)/kT)")


def practical_metadynamics_workflow():
    """
    Practical workflow for running metadynamics simulations.
    """
    print("\n" + "=" * 70)
    print("Practical Metadynamics Workflow")
    print("=" * 70)

    print("\n1. Coordinate Selection:")
    print("   ✓ Choose collective variables wisely")
    print("   ✓ Should capture slow degrees of freedom")
    print("   ✓ Avoid too many dimensions (curse of dimensionality)")
    print("   ✓ Test with short exploratory runs")

    print("\n2. Parameter Tuning:")
    print("   Hill height W:")
    print("     - Start small (~kT)")
    print("     - Monitor barrier crossing")
    print("     - Adjust based on barrier height")

    print("\n   Gaussian width σ:")
    print("     - Related to local fluctuations")
    print("     - σ ≈ √(2 k_B T / k_coordinate)")
    print("     - Too small: rough FES")
    print("     - Too large: loss of resolution")

    print("\n   Deposition frequency:")
    print("     - Balance speed vs. accuracy")
    print("     - Too frequent: hysteresis effects")
    print("     - Too rare: slow convergence")
    print("     - Typical: 100-1000 steps")

    print("\n3. Convergence Checks:")
    print("   ✓ Monitor free energy evolution")
    print("   ✓ Check barrier heights stabilize")
    print("   ✓ Observe multiple transitions")
    print("   ✓ Compare forward/backward paths")

    print("\n4. Production Run:")
    print("   ✓ Long enough for convergence")
    print("   ✓ Save bias potential periodically")
    print("   ✓ Monitor visited grid points")
    print("   ✓ Typical: 50-500 ns depending on system")

    print("\n5. Analysis:")
    print("   ✓ Reconstruct free energy surface")
    print("   ✓ Locate minima and transition states")
    print("   ✓ Calculate barrier heights")
    print("   ✓ Extract minimum energy paths")
    print("   ✓ Estimate errors")


def coordinate_types_reference():
    """
    Reference for all available coordinate types.
    """
    print("\n" + "=" * 70)
    print("Available Coordinate Types")
    print("=" * 70)

    print("\n1. Distance (2 atoms):")
    print("   coord = gromos.CoordinateType.distance()")
    print("   s = |r_j - r_i|")

    print("\n2. Angle (3 atoms):")
    print("   coord = gromos.CoordinateType.angle()")
    print("   s = ∠(i-j-k)")

    print("\n3. Dihedral (4 atoms):")
    print("   coord = gromos.CoordinateType.dihedral()")
    print("   s = dihedral angle around j-k")

    print("\n4. RMSD (N atoms + reference):")
    print("   coord = gromos.CoordinateType.rmsd()")
    print("   s = RMSD from reference structure")

    print("\n5. Coordinate (1 atom, 1 dimension):")
    print("   coord = gromos.CoordinateType.coordinate()")
    print("   s = x, y, or z coordinate")

    print("\n6. Distance Difference (4 atoms):")
    print("   coord = gromos.CoordinateType.distance_difference()")
    print("   s = |r_k - r_j| - |r_j - r_i|")

    print("\n7. Radius of Gyration (N atoms):")
    print("   coord = gromos.CoordinateType.radius_of_gyration()")
    print("   s = R_g of selected atoms")


def example_simulation_output():
    """
    Example output and what to expect from metadynamics.
    """
    print("\n" + "=" * 70)
    print("Expected Simulation Output")
    print("=" * 70)

    print("\nOutput files:")
    print("  1. bias_potential.dat")
    print("     - Grid points and bias values")
    print("     - Updated every N steps")
    print("     - Use for FES reconstruction")

    print("\n  2. hills.dat")
    print("     - Time, position, height of each hill")
    print("     - Can rebuild V_bias from scratch")
    print("     - Useful for analysis")

    print("\n  3. colvar.dat")
    print("     - Time series of CV values")
    print("     - Monitor sampling")
    print("     - Check diffusion in CV space")

    print("\n  4. fes.dat")
    print("     - Free energy surface")
    print("     - F(s) vs. s")
    print("     - Final result")

    print("\nTypical FES for dihedral:")
    print("  φ (deg)    F(φ) (kJ/mol)")
    print("  -180       10.5")
    print("  -120       2.3  ← β-sheet minimum")
    print("   -60       0.0  ← α-helix minimum (reference)")
    print("     0       8.7")
    print("   +60       12.4")
    print("  +120       15.1")
    print("  +180       10.5")


def main():
    """Run all metadynamics examples."""
    print("\n" + "=" * 70)
    print("GROMOS Local Elevation/Metadynamics - Complete Example")
    print("=" * 70)

    # 1D examples
    umb_dihedral = simple_1d_metadynamics()
    umb_distance = distance_metadynamics()

    # 2D example
    umb_2d = metadynamics_2d()

    # RMSD example
    umb_rmsd = rmsd_metadynamics()

    # Advanced features
    advanced_metadynamics_features()

    # Reconstruction
    free_energy_reconstruction()

    # Workflow
    practical_metadynamics_workflow()

    # Coordinate types
    coordinate_types_reference()

    # Expected output
    example_simulation_output()

    print("\n" + "=" * 70)
    print("Example complete!")
    print("=" * 70)
    print("\nKey takeaways:")
    print("  1. Metadynamics fills energy basins with Gaussian hills")
    print("  2. Choose collective variables carefully")
    print("  3. Tune hill height, width, and frequency")
    print("  4. Monitor convergence continuously")
    print("  5. 2D/3D metadynamics scales as O(N^d)")
    print("\nApplications:")
    print("  • Protein folding/unfolding")
    print("  • Conformational transitions")
    print("  • Ligand binding/unbinding")
    print("  • Chemical reactions")
    print("  • Crystal nucleation")
    print("  • Membrane permeation")
    print("\nReferences:")
    print("  - Laio & Parrinello, PNAS 2002")
    print("  - Barducci et al., Phys Rev Lett 2008 (well-tempered)")
    print("  - Valsson & Parrinello, Phys Rev Lett 2014 (variationally optimized)")


if __name__ == "__main__":
    main()
