"""
Basic tests for GROMOS-RS Python bindings.

Tests only what the Rust extension actually exports:
  Vec3, Energy, Frame, rmsd, rdf
"""

import numpy as np
import pytest

import gromos


class TestImport:
    def test_import(self):
        assert gromos.__version__ == "0.1.0"

    def test_available_types(self):
        assert hasattr(gromos, "Vec3")
        assert hasattr(gromos, "Energy")
        assert hasattr(gromos, "Frame")
        assert hasattr(gromos, "rmsd")
        assert hasattr(gromos, "rdf")


class TestVec3:
    def test_creation(self):
        v = gromos.Vec3(1.0, 2.0, 3.0)
        assert v.x == 1.0
        assert v.y == 2.0
        assert v.z == 3.0

    def test_repr(self):
        v = gromos.Vec3(1.0, 2.0, 3.0)
        r = repr(v)
        assert "1" in r and "2" in r and "3" in r

    def test_addition(self):
        v1 = gromos.Vec3(1.0, 0.0, 0.0)
        v2 = gromos.Vec3(0.0, 1.0, 0.0)
        v3 = v1 + v2
        assert v3.x == pytest.approx(1.0)
        assert v3.y == pytest.approx(1.0)
        assert v3.z == pytest.approx(0.0)

    def test_subtraction(self):
        v1 = gromos.Vec3(3.0, 2.0, 1.0)
        v2 = gromos.Vec3(1.0, 1.0, 1.0)
        v3 = v1 - v2
        assert v3.x == pytest.approx(2.0)
        assert v3.y == pytest.approx(1.0)
        assert v3.z == pytest.approx(0.0)

    def test_scalar_mul(self):
        v = gromos.Vec3(1.0, 2.0, 3.0)
        v2 = v * 2.0
        assert v2.x == pytest.approx(2.0)
        assert v2.y == pytest.approx(4.0)
        assert v2.z == pytest.approx(6.0)

    def test_dot(self):
        v1 = gromos.Vec3(1.0, 0.0, 0.0)
        v2 = gromos.Vec3(0.0, 1.0, 0.0)
        assert v1.dot(v2) == pytest.approx(0.0)

        v3 = gromos.Vec3(1.0, 2.0, 3.0)
        v4 = gromos.Vec3(4.0, 5.0, 6.0)
        assert v3.dot(v4) == pytest.approx(32.0)

    def test_cross(self):
        v1 = gromos.Vec3(1.0, 0.0, 0.0)
        v2 = gromos.Vec3(0.0, 1.0, 0.0)
        c = v1.cross(v2)
        assert c.x == pytest.approx(0.0)
        assert c.y == pytest.approx(0.0)
        assert c.z == pytest.approx(1.0)

    def test_length(self):
        v = gromos.Vec3(3.0, 4.0, 0.0)
        assert v.length() == pytest.approx(5.0)

    def test_normalize(self):
        v = gromos.Vec3(3.0, 0.0, 0.0)
        n = v.normalize()
        assert n.x == pytest.approx(1.0)
        assert n.y == pytest.approx(0.0)
        assert n.length() == pytest.approx(1.0)


class TestEnergy:
    def test_creation(self):
        e = gromos.Energy()
        assert e.kinetic == pytest.approx(0.0)
        assert e.potential == pytest.approx(0.0)
        assert e.total == pytest.approx(0.0)


class TestFrame:
    def test_creation(self):
        f = gromos.Frame(1.5, 100)
        assert f.time == pytest.approx(1.5)
        assert f.step == 100
        assert f.n_atoms == 0


class TestRmsd:
    def test_identical(self):
        pos = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=np.float32)
        assert gromos.rmsd(pos, pos) == pytest.approx(0.0, abs=1e-6)

    def test_translated(self):
        pos = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=np.float32)
        ref = pos + 1.0
        r = gromos.rmsd(pos, ref)
        assert r == pytest.approx(np.sqrt(3.0), rel=1e-5)

    def test_shape_mismatch(self):
        pos = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.float32)
        ref = np.array([[0, 0, 0]], dtype=np.float32)
        with pytest.raises(ValueError):
            gromos.rmsd(pos, ref)


class TestRdf:
    def test_basic(self):
        # 3 atoms at known distances
        pos = np.array(
            [
                [0, 0, 0],
                [1, 0, 0],
                [0, 2, 0],
            ],
            dtype=np.float32,
        )
        r_vals, g_vals = gromos.rdf(pos, [0], [1, 2], n_bins=10, r_max=3.0)
        assert len(r_vals) == 10
        assert len(g_vals) == 10
        # g(r) should have non-zero values at distance bins containing 1.0 and 2.0
        assert any(g > 0 for g in g_vals)
