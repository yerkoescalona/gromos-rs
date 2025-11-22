"""
Comprehensive tests for all new py-gromos features.

Tests:
- NMR restraints (J-value, RDC)
- Virtual atoms
- Local elevation/metadynamics
- Polarization
- QM/MM
- Free energy perturbation
- Analysis programs
"""

import pytest
import gromos
import gromos.analysis as analysis


class TestNMRRestraints:
    """Test NMR restraint classes."""

    def test_jvalue_restraint(self):
        """Test J-value restraint creation."""
        jval = gromos.JValueRestraint(
            atoms=[4, 6, 8, 14],
            target_j=4.2,
            force_constant=10.0,
            karplus_a=6.4,
            karplus_b=-1.4,
            karplus_c=1.9,
        )
        assert jval.target_j == 4.2
        assert jval.force_constant == 10.0
        assert len(jval.atoms) == 4

    def test_rdc_restraint(self):
        """Test RDC restraint creation."""
        rdc = gromos.RDCRestraint(
            i=5,
            j=6,
            d_max=-10000.0,
            target_rdc=15.2,
            force_constant=5.0,
            saupe_matrix=[0.00015, -0.00008, -0.00007, 0.00003, 0.00001],
        )
        assert rdc.target_rdc == 15.2
        assert rdc.d_max == -10000.0


class TestVirtualAtoms:
    """Test virtual atom classes."""

    def test_virtual_atom_creation(self):
        """Test creating virtual atom."""
        va = gromos.VirtualAtom(
            virtual_index=3, construct_atoms=[0, 1, 2], construct_type=1, distance_om=0.015
        )
        assert va.virtual_index == 3
        assert len(va.construct_atoms) == 3
        assert va.distance_om == 0.015

    def test_virtual_atom_manager(self):
        """Test virtual atom manager."""
        manager = gromos.VirtualAtomManager()

        va1 = gromos.VirtualAtom(
            virtual_index=3, construct_atoms=[0, 1, 2], construct_type=1, distance_om=0.015
        )
        manager.add_virtual_atom(va1)

        assert manager.num_virtual_atoms() == 1


class TestLocalElevation:
    """Test local elevation/metadynamics classes."""

    def test_coordinate_types(self):
        """Test coordinate type creation."""
        dihedral = gromos.CoordinateType.dihedral()
        distance = gromos.CoordinateType.distance()
        angle = gromos.CoordinateType.angle()
        rmsd = gromos.CoordinateType.rmsd()

        # Just verify they can be created
        assert dihedral is not None
        assert distance is not None

    def test_le_coordinate(self):
        """Test LE coordinate creation."""
        coord = gromos.LECoordinate(
            umbrella_id=1, coord_type=gromos.CoordinateType.dihedral(), atoms=[4, 6, 8, 14]
        )
        assert coord.umbrella_id == 1
        assert len(coord.atoms) == 4

    def test_umbrella(self):
        """Test umbrella creation."""
        coord = gromos.LECoordinate(
            umbrella_id=1, coord_type=gromos.CoordinateType.dihedral(), atoms=[4, 6, 8, 14]
        )

        umbrella = gromos.Umbrella(
            umbrella_id=1,
            coordinates=[coord],
            grid_sizes=[360],
            grid_mins=[-180.0],
            grid_maxs=[180.0],
            grid_spacings=[1.0],
            hill_height=0.5,
            gaussian_widths=[5.0],
            deposition_frequency=100,
        )

        assert umbrella.dimensionality == 1
        assert umbrella.hill_height == 0.5


class TestPolarization:
    """Test polarization classes."""

    def test_polarization_models(self):
        """Test polarization model types."""
        none_model = gromos.PolarizationModel.none()
        point_dipole = gromos.PolarizationModel.point_dipole()
        shell = gromos.PolarizationModel.shell_model()
        drude = gromos.PolarizationModel.drude_oscillator()
        fluct = gromos.PolarizationModel.fluctuating_charge()

        assert none_model is not None
        assert point_dipole is not None

    def test_polarizability_parameters(self):
        """Test polarizability parameters."""
        params = gromos.PolarizabilityParameters(alpha=1.0, thole_a=2.6, drude_mass=0.4)
        assert params.alpha == 1.0
        assert params.thole_a == 2.6

    def test_polarization_calculator(self):
        """Test polarization calculator."""
        calc = gromos.PolarizationCalculator(
            n_atoms=100, model=gromos.PolarizationModel.point_dipole()
        )
        assert calc.num_atoms() == 0  # No atoms added yet

        params = gromos.PolarizabilityParameters(alpha=1.0)
        calc.add_atom(params)
        assert calc.num_atoms() == 1


class TestFreeEnergyPerturbation:
    """Test FEP classes."""

    def test_lambda_controller(self):
        """Test lambda controller."""
        lambda_ctrl = gromos.LambdaController(
            lambda_value=0.0, dlambda_dt=0.0001, lambda_exponent=1.0
        )
        assert lambda_ctrl.current_lambda == 0.0

    def test_perturbed_atom(self):
        """Test perturbed atom."""
        pert = gromos.PerturbedAtom(
            atom_index=10,
            a_charge=0.0,
            b_charge=0.0,
            a_iac=12,
            b_iac=15,
            lj_softcore=0.5,
            crf_softcore=0.5,
        )
        assert pert.charge_at_lambda(0.0) == 0.0
        assert pert.charge_at_lambda(1.0) == 0.0
        assert pert.charge_at_lambda(0.5) == 0.0

    def test_perturbed_atom_charged(self):
        """Test charged residue perturbation."""
        pert = gromos.PerturbedAtom(
            atom_index=25,
            a_charge=+1.0,
            b_charge=-1.0,
            a_iac=20,
            b_iac=21,
            lj_softcore=0.5,
            crf_softcore=0.5,
        )
        assert pert.charge_at_lambda(0.0) == +1.0
        assert pert.charge_at_lambda(1.0) == -1.0
        assert abs(pert.charge_at_lambda(0.5)) < 0.01  # Should be ~0


class TestAnalysisPrograms:
    """Test that analysis programs are accessible."""

    def test_analysis_module_exists(self):
        """Test analysis module can be imported."""
        assert analysis is not None

    def test_structural_programs(self):
        """Test structural analysis functions exist."""
        assert hasattr(analysis, "hbond")
        assert hasattr(analysis, "sasa")
        assert hasattr(analysis, "cluster")
        assert hasattr(analysis, "dssp")
        assert hasattr(analysis, "distmat")

    def test_interaction_programs(self):
        """Test interaction analysis functions exist."""
        assert hasattr(analysis, "rdf")
        assert hasattr(analysis, "contactnum")

    def test_dynamics_programs(self):
        """Test dynamics analysis functions exist."""
        assert hasattr(analysis, "diffus")
        assert hasattr(analysis, "dipole")
        assert hasattr(analysis, "visco")
        assert hasattr(analysis, "tcf")
        assert hasattr(analysis, "rot_rel")

    def test_energy_programs(self):
        """Test energy analysis functions exist."""
        assert hasattr(analysis, "ene_ana")
        assert hasattr(analysis, "pot_aver")
        assert hasattr(analysis, "int_ener")

    def test_free_energy_programs(self):
        """Test free energy analysis functions exist."""
        assert hasattr(analysis, "bar")
        assert hasattr(analysis, "ext_ti_ana")
        assert hasattr(analysis, "m_widom")
        assert hasattr(analysis, "reweight")

    def test_xray_nmr_programs(self):
        """Test X-ray/NMR analysis functions exist."""
        assert hasattr(analysis, "xray_map")
        assert hasattr(analysis, "noe")
        assert hasattr(analysis, "structure_factor")
        assert hasattr(analysis, "r_factor")

    def test_generic_runner(self):
        """Test generic program runner exists."""
        assert hasattr(analysis, "run_program")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
