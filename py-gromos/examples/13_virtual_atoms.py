"""
Virtual Atoms Example
=====================

This example demonstrates the use of virtual atoms (dummy atoms) in GROMOS:
- Coarse-grained water models (TIP4P, TIP5P)
- Martini coarse-grained force fields
- Virtual sites for constraints
- Lone pair representations

Virtual atoms have no mass and their positions are constructed geometrically
from real atoms, enabling accurate representations of molecular interactions.
"""

import gromos
import numpy as np


def tip4p_water_example():
    """
    TIP4P water model with virtual M-site.

    TIP4P places a virtual site (M) on the bisector of the H-O-H angle,
    displaced from the oxygen towards the hydrogens. The negative charge
    resides on M instead of O, improving the electrostatic representation.

    Geometry:
        H1
         \
          O --- M (virtual site)
         /
        H2

    M is located at: r_OM = 0.15 Å along the H-H bisector
    """
    print("=" * 70)
    print("TIP4P Water Model - Virtual M-Site")
    print("=" * 70)

    # Create TIP4P water virtual atom manager
    virt_manager = gromos.VirtualAtomManager()

    # TIP4P geometry parameters
    # O-H distance: 0.9572 Å
    # H-O-H angle: 104.52°
    # O-M distance: 0.15 Å

    # Create virtual M-site for first water molecule
    # Atom indices: O=0, H1=1, H2=2, M=3 (virtual)
    tip4p_m1 = gromos.VirtualAtom(
        virtual_index=3,  # Index of virtual M-site
        construct_atoms=[0, 1, 2],  # O, H1, H2
        construct_type=1,  # Type 1: Out-of-plane (TIP4P geometry)
        distance_om=0.015,  # 0.15 Å = 0.015 nm from O
    )

    print(f"\nTIP4P M-site 1: {tip4p_m1}")
    print(f"  Virtual atom index: {tip4p_m1.virtual_index}")
    print(f"  Constructed from atoms: {tip4p_m1.construct_atoms} (O, H1, H2)")
    print(f"  O-M distance: {tip4p_m1.distance_om * 10:.2f} Å")

    # Add to manager
    virt_manager.add_virtual_atom(tip4p_m1)

    # Create virtual sites for multiple water molecules
    # Water 2: atoms 4-7 (O, H1, H2, M)
    tip4p_m2 = gromos.VirtualAtom(
        virtual_index=7, construct_atoms=[4, 5, 6], construct_type=1, distance_om=0.015
    )
    virt_manager.add_virtual_atom(tip4p_m2)

    # Water 3: atoms 8-11
    tip4p_m3 = gromos.VirtualAtom(
        virtual_index=11, construct_atoms=[8, 9, 10], construct_type=1, distance_om=0.015
    )
    virt_manager.add_virtual_atom(tip4p_m3)

    print(f"\nTotal TIP4P water molecules: 3")
    print(f"Virtual atoms managed: {virt_manager.num_virtual_atoms()}")

    print("\nTIP4P advantages:")
    print("  ✓ Better electrostatics than TIP3P")
    print("  ✓ More accurate water structure")
    print("  ✓ Improved dielectric properties")
    print("  ✓ No additional computational atoms (M is virtual)")

    return virt_manager


def tip5p_water_example():
    """
    TIP5P water model with two virtual L-sites.

    TIP5P places two virtual sites (L1, L2) representing lone pairs
    on oxygen, arranged tetrahedrally. This gives an even better
    representation of water's electronic structure.

    Geometry:
         L1
          |
       H1-O-H2
          |
         L2

    L-sites arranged tetrahedrally, ~0.70 Å from O
    """
    print("\n" + "=" * 70)
    print("TIP5P Water Model - Virtual Lone Pairs")
    print("=" * 70)

    virt_manager = gromos.VirtualAtomManager()

    # TIP5P: 2 virtual sites per water
    # Atom indices: O=0, H1=1, H2=2, L1=3, L2=4

    # First lone pair (L1)
    tip5p_l1 = gromos.VirtualAtom(
        virtual_index=3,
        construct_atoms=[0, 1, 2],  # O, H1, H2
        construct_type=2,  # Type 2: Tetrahedral lone pair
        distance_om=0.070,  # 0.70 Å from O
    )

    # Second lone pair (L2)
    tip5p_l2 = gromos.VirtualAtom(
        virtual_index=4, construct_atoms=[0, 1, 2], construct_type=2, distance_om=0.070
    )

    virt_manager.add_virtual_atom(tip5p_l1)
    virt_manager.add_virtual_atom(tip5p_l2)

    print(f"\nTIP5P lone pair sites: 2 per water")
    print(f"  L1: {tip5p_l1}")
    print(f"  L2: {tip5p_l2}")
    print(f"  Distance from O: {tip5p_l1.distance_om * 10:.2f} Å")

    print("\nTIP5P advantages:")
    print("  ✓ Most accurate TIP model for water structure")
    print("  ✓ Excellent reproduction of liquid water properties")
    print("  ✓ Better temperature dependence")
    print("  ✗ More complex setup (2 virtual sites per water)")

    return virt_manager


def martini_cg_example():
    """
    Martini coarse-grained force field virtual sites.

    Martini uses virtual sites to represent anisotropic particles
    and improve directionality in coarse-grained models.

    Example: Benzene ring represented by 3 beads + virtual sites
    """
    print("\n" + "=" * 70)
    print("Martini Coarse-Grained Virtual Sites")
    print("=" * 70)

    virt_manager = gromos.VirtualAtomManager()

    # Martini benzene: 3 SC4 beads in triangle + 1 virtual site at center
    # This maintains ring planarity and improves stacking interactions

    # Virtual site at center of benzene ring
    # Constructed from 3 ring beads
    benzene_center = gromos.VirtualAtom(
        virtual_index=3,  # Central virtual site
        construct_atoms=[0, 1, 2],  # 3 SC4 beads forming triangle
        construct_type=3,  # Type 3: Center of mass construction
        distance_om=0.0,  # At geometric center
    )

    print(f"\nMartini benzene virtual site: {benzene_center}")
    print(f"  Constructed from 3 SC4 beads")
    print(f"  Location: Geometric center of triangle")

    virt_manager.add_virtual_atom(benzene_center)

    # Martini cholesterol: Multiple virtual sites for ring system
    # Rings represented by beads + virtual sites for rigidity
    chol_vs1 = gromos.VirtualAtom(
        virtual_index=10,
        construct_atoms=[4, 5, 6],  # Ring A beads
        construct_type=3,
        distance_om=0.0,
    )

    chol_vs2 = gromos.VirtualAtom(
        virtual_index=11,
        construct_atoms=[7, 8, 9],  # Ring B beads
        construct_type=3,
        distance_om=0.0,
    )

    virt_manager.add_virtual_atom(chol_vs1)
    virt_manager.add_virtual_atom(chol_vs2)

    print(f"\nMartini cholesterol virtual sites: 2")
    print(f"  Ring A center: {chol_vs1.virtual_index}")
    print(f"  Ring B center: {chol_vs2.virtual_index}")

    print(f"\nTotal Martini virtual sites: {virt_manager.num_virtual_atoms()}")

    print("\nMartini virtual site benefits:")
    print("  ✓ Maintain molecular rigidity")
    print("  ✓ Improve directional interactions")
    print("  ✓ Enable proper π-π stacking")
    print("  ✓ No mass (purely geometric)")

    return virt_manager


def constraint_virtual_sites():
    """
    Virtual sites for constraint algorithms.

    Virtual sites can be used to:
    1. Represent rigid bodies (e.g., rigid water)
    2. Implement holonomic constraints
    3. Define constraint reference points
    """
    print("\n" + "=" * 70)
    print("Virtual Sites for Constraints")
    print("=" * 70)

    virt_manager = gromos.VirtualAtomManager()

    # Rigid water: All atoms positions defined relative to COM virtual site
    # This allows treating water as a rigid body

    # Virtual site at water center of mass
    water_com = gromos.VirtualAtom(
        virtual_index=3,
        construct_atoms=[0, 1, 2],  # O, H1, H2
        construct_type=3,  # Center of mass
        distance_om=0.0,
    )

    print(f"\nRigid water COM virtual site: {water_com}")
    print("  Used for:")
    print("    - SETTLE algorithm (rigid water)")
    print("    - Faster integration (3 DOF instead of 9)")
    print("    - Exact conservation of bond lengths/angles")

    virt_manager.add_virtual_atom(water_com)

    # Constraint reference point for distance constraints
    # E.g., for pulling simulations or umbrella sampling
    constraint_ref = gromos.VirtualAtom(
        virtual_index=100,
        construct_atoms=[10, 11, 12, 13],  # 4 atoms defining reference
        construct_type=3,  # Geometric center
        distance_om=0.0,
    )

    print(f"\nConstraint reference point: {constraint_ref}")
    print("  Applications:")
    print("    - Umbrella sampling reference")
    print("    - Steered MD pulling point")
    print("    - Flexible distance constraints")

    virt_manager.add_virtual_atom(constraint_ref)

    return virt_manager


def advanced_virtual_atom_types():
    """
    Different virtual atom construction types.
    """
    print("\n" + "=" * 70)
    print("Virtual Atom Construction Types")
    print("=" * 70)

    print("\nType 0: Linear (2 atoms)")
    print("  V is on line between atoms A and B")
    print("  Position: r_V = r_A + d * (r_B - r_A) / |r_B - r_A|")
    print("  Example: Bond midpoint, extended bond")

    linear_vs = gromos.VirtualAtom(
        virtual_index=10,
        construct_atoms=[0, 1],
        construct_type=0,
        distance_om=0.05,  # 0.5 Å from atom 0 toward atom 1
    )
    print(f"  {linear_vs}")

    print("\nType 1: Out-of-plane (3 atoms)")
    print("  V is on bisector of angle ABC, displaced along bisector")
    print("  Position: On angle bisector at distance d from B")
    print("  Example: TIP4P M-site, trigonal planar")

    oop_vs = gromos.VirtualAtom(
        virtual_index=11, construct_atoms=[0, 1, 2], construct_type=1, distance_om=0.015
    )
    print(f"  {oop_vs}")

    print("\nType 2: Tetrahedral (3 atoms)")
    print("  V is at tetrahedral position relative to 3 atoms")
    print("  Position: Tetrahedral geometry, distance d from central atom")
    print("  Example: TIP5P lone pairs, sp3 representation")

    tetra_vs = gromos.VirtualAtom(
        virtual_index=12, construct_atoms=[0, 1, 2], construct_type=2, distance_om=0.070
    )
    print(f"  {tetra_vs}")

    print("\nType 3: Center of Mass (N atoms)")
    print("  V is at geometric center of N atoms")
    print("  Position: r_V = (1/N) * Σ r_i")
    print("  Example: Molecular COM, ring centers")

    com_vs = gromos.VirtualAtom(
        virtual_index=13,
        construct_atoms=[0, 1, 2, 3, 4],  # 5 atoms
        construct_type=3,
        distance_om=0.0,
    )
    print(f"  {com_vs}")

    print("\nType 4: Weighted COM (N atoms with masses)")
    print("  V is at mass-weighted center")
    print("  Position: r_V = Σ (m_i * r_i) / Σ m_i")
    print("  Example: True center of mass with masses")


def virtual_atom_workflow():
    """
    Complete workflow for using virtual atoms in simulation.
    """
    print("\n" + "=" * 70)
    print("Virtual Atoms Workflow")
    print("=" * 70)

    print("\n1. Setup Phase:")
    print("   a. Create VirtualAtomManager")
    manager = gromos.VirtualAtomManager()

    print("   b. Define virtual atoms with construction rules")
    va1 = gromos.VirtualAtom(
        virtual_index=3, construct_atoms=[0, 1, 2], construct_type=1, distance_om=0.015
    )
    manager.add_virtual_atom(va1)

    print("   c. Add all virtual atoms to manager")
    print(f"      Total virtual atoms: {manager.num_virtual_atoms()}")

    print("\n2. Integration Phase:")
    print("   a. Before force calculation:")
    print("      - Construct virtual atom positions from parent atoms")
    print("      - Virtual atoms have NO forces initially")
    print("   b. During force calculation:")
    print("      - Calculate forces on virtual atoms (if charged)")
    print("      - These are purely electrostatic/vdW")
    print("   c. After force calculation:")
    print("      - Distribute virtual atom forces to parent atoms")
    print("      - Maintain conservation of momentum")

    print("\n3. Constraints:")
    print("   - Virtual atoms automatically satisfy geometry constraints")
    print("   - No SHAKE/RATTLE needed for virtual-parent distances")
    print("   - Improves numerical stability")

    print("\n4. Output:")
    print("   - Virtual atom positions written to trajectory")
    print("   - Can be visualized like regular atoms")
    print("   - Useful for analysis (e.g., charge distribution)")


def practical_tips():
    """
    Practical tips for using virtual atoms.
    """
    print("\n" + "=" * 70)
    print("Practical Tips & Best Practices")
    print("=" * 70)

    print("\n✓ Water Models:")
    print("  - TIP3P: No virtual atoms (3 atoms)")
    print("  - TIP4P: 1 virtual M-site (4 atoms total)")
    print("  - TIP5P: 2 virtual L-sites (5 atoms total)")
    print("  → Choose based on accuracy vs. speed trade-off")

    print("\n✓ Topology Setup:")
    print("  - Virtual atoms must be LAST in atom list")
    print("  - Number virtual atoms starting after real atoms")
    print("  - Exclude virtual-parent interactions from pair list")

    print("\n✓ Force Distribution:")
    print("  - Virtual atom forces distributed to parents")
    print("  - Preserves force conservation")
    print("  - Automatic in GROMOS virtual atom algorithm")

    print("\n✓ Performance:")
    print("  - Virtual atoms: ~5-10% overhead for construction")
    print("  - TIP4P: ~10% slower than TIP3P")
    print("  - TIP5P: ~15% slower than TIP3P")
    print("  → Worthwhile for improved accuracy")

    print("\n✓ Analysis:")
    print("  - Include virtual atoms in RDF calculations")
    print("  - Useful for charge distribution visualization")
    print("  - Better representation of electrostatics")

    print("\n⚠ Common Pitfalls:")
    print("  - Forgetting to exclude virtual-parent pairs")
    print("  - Incorrect atom ordering in topology")
    print("  - Not updating virtual positions before analysis")
    print("  - Using virtual atoms in bond/angle terms (invalid!)")


def main():
    """Run all virtual atom examples."""
    print("\n" + "=" * 70)
    print("GROMOS Virtual Atoms - Complete Example")
    print("=" * 70)

    # TIP4P water
    tip4p_mgr = tip4p_water_example()

    # TIP5P water
    tip5p_mgr = tip5p_water_example()

    # Martini CG
    martini_mgr = martini_cg_example()

    # Constraints
    constraint_mgr = constraint_virtual_sites()

    # Advanced types
    advanced_virtual_atom_types()

    # Workflow
    virtual_atom_workflow()

    # Tips
    practical_tips()

    print("\n" + "=" * 70)
    print("Example complete!")
    print("=" * 70)
    print("\nKey takeaways:")
    print("  1. Virtual atoms have no mass - purely geometric")
    print("  2. TIP4P/TIP5P improve water electrostatics")
    print("  3. Martini uses virtual sites for rigidity")
    print("  4. Multiple construction types for different geometries")
    print("  5. Forces distributed to parent atoms automatically")
    print("\nApplications:")
    print("  • Accurate water models (TIP4P, TIP5P)")
    print("  • Coarse-grained simulations (Martini)")
    print("  • Rigid body constraints")
    print("  • Improved charge distributions")
    print("  • Lone pair representations")


if __name__ == "__main__":
    main()
