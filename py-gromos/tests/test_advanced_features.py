"""
Tests for advanced py-gromos features.

These tests are placeholders for when the Rust bindings are extended.
Currently only Vec3, Energy, Frame, rmsd, rdf are exposed.
"""

import pytest


@pytest.mark.skip(reason="NMR restraints not yet exposed in Rust bindings")
class TestNMRRestraints:
    def test_jvalue_restraint(self):
        pass

    def test_rdc_restraint(self):
        pass


@pytest.mark.skip(reason="Virtual atoms not yet exposed in Rust bindings")
class TestVirtualAtoms:
    def test_virtual_atom_creation(self):
        pass


@pytest.mark.skip(reason="Local elevation not yet exposed in Rust bindings")
class TestLocalElevation:
    def test_coordinate_types(self):
        pass


@pytest.mark.skip(reason="Polarization not yet exposed in Rust bindings")
class TestPolarization:
    def test_polarization_model(self):
        pass


@pytest.mark.skip(reason="QM/MM not yet exposed in Rust bindings")
class TestQMMM:
    def test_qmmm_calculator(self):
        pass


@pytest.mark.skip(reason="FEP not yet exposed in Rust bindings")
class TestFreeEnergyPerturbation:
    def test_lambda_controller(self):
        pass

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
        assert umbrella.building == True


@pytest.mark.skip(reason="Polarization not yet exposed in Rust bindings")
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
        # Note: num_atoms() returns the initialized size, not params added
        assert calc.num_atoms() == 100
        # Add atoms
        params = gromos.PolarizabilityParameters(alpha=1.0, thole_a=2.6)
        calc.add_atom(params)
        # After adding, num_atoms still reflects total capacity
        assert calc.num_atoms() == 101  # Added one more


@pytest.mark.skip(reason="FEP not yet exposed in Rust bindings")
class TestFreeEnergyPerturbation:
    """Test FEP classes."""

    def test_lambda_controller(self):
        """Test lambda controller."""
        lambda_ctrl = gromos.LambdaController(
            lambda_value=0.0, lambda_exponent=1, dlambda_dt=0.0001
        )
        assert lambda_ctrl.lambda_value == 0.0
        assert lambda_ctrl.lambda_exponent == 1

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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
