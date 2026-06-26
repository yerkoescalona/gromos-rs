//! Improper dihedral and cross-dihedral force calculations.

use gromos_core::configuration::Configuration;
use gromos_core::math::Vec3;
use gromos_core::topology::Topology;

use super::ForceEnergy;

/// Calculate improper dihedral forces (for planarity/chirality).
///
/// Potential: V = (1/2) * K * (ζ - ζ₀)²
pub fn calculate_improper_dihedral_forces(topo: &Topology, conf: &Configuration) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    for improper in topo.all_improper_dihedrals_global() {
        if improper.dihedral_type >= topo.improper_dihedral_parameters.len() {
            continue;
        }

        let params = &topo.improper_dihedral_parameters[improper.dihedral_type];

        // gromosXX convention: rkj = pos(k)-pos(j), rij = pos(i)-pos(j), rkl = pos(k)-pos(l)
        let r_kj = conf.current().pos[improper.k] - conf.current().pos[improper.j];
        let r_ij = conf.current().pos[improper.i] - conf.current().pos[improper.j];
        let r_kl = conf.current().pos[improper.k] - conf.current().pos[improper.l];

        let r_mj = r_ij.cross(r_kj);
        let r_nk = r_kj.cross(r_kl);

        let d_kj2 = r_kj.dot(r_kj);
        let d_mj2 = r_mj.dot(r_mj);
        let d_nk2 = r_nk.dot(r_nk);

        let d_kj = d_kj2.sqrt();
        let d_mj = d_mj2.sqrt();
        let d_nk = d_nk2.sqrt();

        if d_mj < 1e-10 || d_nk < 1e-10 {
            continue;
        }

        let acs = (r_mj.dot(r_nk) / (d_mj * d_nk)).clamp(-1.0, 1.0);
        let mut zeta = acs.acos();

        if r_ij.dot(r_nk) < 0.0 {
            zeta = -zeta;
        }

        // Bring ζ to interval [ζ₀ - π, ζ₀ + π]
        let mut zeta_adj = zeta;
        while zeta_adj < (params.q0 - std::f64::consts::PI) {
            zeta_adj += 2.0 * std::f64::consts::PI;
        }
        while zeta_adj > (params.q0 + std::f64::consts::PI) {
            zeta_adj -= 2.0 * std::f64::consts::PI;
        }

        let d_zeta = zeta_adj - params.q0;
        let energy = 0.5 * params.k * d_zeta * d_zeta;

        let mut k_i = -params.k * d_zeta * d_kj;
        let mut k_l = -k_i;

        if d_mj2 < 1.0e-10 * d_kj2 {
            k_i = 0.0;
        } else {
            k_i /= d_mj2;
        }

        if d_nk2 < 1.0e-10 * d_kj2 {
            k_l = 0.0;
        } else {
            k_l /= d_nk2;
        }

        let k_j1 = r_ij.dot(r_kj) / d_kj2 - 1.0;
        let k_j2 = r_kl.dot(r_kj) / d_kj2;

        let f_i = r_mj * k_i;
        let f_l = r_nk * k_l;
        let f_j = f_i * k_j1 - f_l * k_j2;
        let f_k = -(f_i + f_j + f_l);

        result.energy += energy;
        result.forces[improper.i] += f_i;
        result.forces[improper.j] += f_j;
        result.forces[improper.k] += f_k;
        result.forces[improper.l] += f_l;

        // gromosXX: virial += rij*fi + rkj*fk + rlj*fl
        let r_lj = conf.current().pos[improper.l] - conf.current().pos[improper.j];
        let rij_v = [r_ij.x, r_ij.y, r_ij.z];
        let rkj_v = [r_kj.x, r_kj.y, r_kj.z];
        let rlj_v = [r_lj.x, r_lj.y, r_lj.z];
        let fi_v  = [f_i.x, f_i.y, f_i.z];
        let fk_v  = [f_k.x, f_k.y, f_k.z];
        let fl_v  = [f_l.x, f_l.y, f_l.z];
        for a in 0..3 {
            for b in 0..3 {
                result.virial[a][b] += rij_v[a]*fi_v[b] + rkj_v[a]*fk_v[b] + rlj_v[a]*fl_v[b];
            }
        }
    }

    result
}

/// Calculate cross-dihedral forces (8-atom coupled torsional term).
///
/// Potential: V = K * (1 + cos(m*(φ + ψ) - δ))
pub fn calculate_crossdihedral_forces(topo: &Topology, conf: &Configuration) -> ForceEnergy {
    let mut result = ForceEnergy::new(topo.num_atoms());

    for crossdih in topo.all_cross_dihedrals_global() {
        if crossdih.cross_dihedral_type >= topo.dihedral_parameters.len() {
            continue;
        }

        let params = &topo.dihedral_parameters[crossdih.cross_dihedral_type];

        let r_ab = conf.current().pos[crossdih.b] - conf.current().pos[crossdih.a];
        let r_cb = conf.current().pos[crossdih.b] - conf.current().pos[crossdih.c];
        let r_cd = conf.current().pos[crossdih.d] - conf.current().pos[crossdih.c];

        let r_ef = conf.current().pos[crossdih.f] - conf.current().pos[crossdih.e];
        let r_gf = conf.current().pos[crossdih.f] - conf.current().pos[crossdih.g];
        let r_gh = conf.current().pos[crossdih.h] - conf.current().pos[crossdih.g];

        let r_nc = r_cb.cross(r_cd);
        let d_cb2 = r_cb.dot(r_cb);
        if d_cb2 < 1e-10 { continue; }

        let f_ram = r_ab.dot(r_cb) / d_cb2;
        let f_rdn = r_cd.dot(r_cb) / d_cb2;
        let r_am  = r_ab - r_cb * f_ram;
        let r_dn  = r_cb * f_rdn - r_cd;
        let d_am  = r_am.length();
        let d_dn  = r_dn.length();
        if d_am < 1e-10 || d_dn < 1e-10 { continue; }

        let cos_phi = (r_am.dot(r_dn) / (d_am * d_dn)).clamp(-1.0, 1.0);
        let mut phi  = cos_phi.acos();
        if r_ab.dot(r_nc) < 0.0 { phi = -phi; }

        let r_ng  = r_gf.cross(r_gh);
        let d_gf2 = r_gf.dot(r_gf);
        if d_gf2 < 1e-10 { continue; }

        let f_rem = r_ef.dot(r_gf) / d_gf2;
        let f_rgn = r_gh.dot(r_gf) / d_gf2;
        let r_em  = r_ef - r_gf * f_rem;
        let r_gn  = r_gf * f_rgn - r_gh;
        let d_em  = r_em.length();
        let d_gn  = r_gn.length();
        if d_em < 1e-10 || d_gn < 1e-10 { continue; }

        let cos_psi = (r_em.dot(r_gn) / (d_em * d_gn)).clamp(-1.0, 1.0);
        let mut psi  = cos_psi.acos();
        if r_ef.dot(r_ng) < 0.0 { psi = -psi; }

        let arg    = params.m as f64 * (phi + psi) - params.pd;
        let energy = params.k * (1.0 + arg.cos());
        let k      = params.k * params.m as f64 * arg.sin();

        let r_mb  = r_ab.cross(r_cb);
        let r_mf  = r_ef.cross(r_gf);
        let d_mb2 = r_mb.dot(r_mb);
        let d_mf2 = r_mf.dot(r_mf);
        let d_nc2 = r_nc.dot(r_nc);
        let d_ng2 = r_ng.dot(r_ng);
        let d_cb  = d_cb2.sqrt();
        let d_gf  = d_gf2.sqrt();

        let mut f_a = Vec3::ZERO;
        let mut f_d = Vec3::ZERO;
        if d_mb2 > 1e-10 * d_cb2 { f_a = r_mb * ( k * d_cb / d_mb2); }
        if d_nc2 > 1e-10 * d_cb2 { f_d = r_nc * (-k * d_cb / d_nc2); }
        let f_b = f_a * (f_ram - 1.0) - f_d * f_rdn;
        let f_c = -(f_a + f_b + f_d);

        let mut f_e = Vec3::ZERO;
        let mut f_h = Vec3::ZERO;
        if d_mf2 > 1e-10 * d_gf2 { f_e = r_mf * ( k * d_gf / d_mf2); }
        if d_ng2 > 1e-10 * d_gf2 { f_h = r_ng * (-k * d_gf / d_ng2); }
        let f_f = f_e * (f_rem - 1.0) - f_h * f_rgn;
        let f_g = -(f_e + f_f + f_h);

        result.energy += energy;
        result.forces[crossdih.a] += f_a;
        result.forces[crossdih.b] += f_b;
        result.forces[crossdih.c] += f_c;
        result.forces[crossdih.d] += f_d;
        result.forces[crossdih.e] += f_e;
        result.forces[crossdih.f] += f_f;
        result.forces[crossdih.g] += f_g;
        result.forces[crossdih.h] += f_h;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use gromos_core::configuration::Configuration;
    use gromos_core::math::Vec3;
    use gromos_core::topology::{Atom, CrossDihedral, DihedralParameters,
                                ImproperDihedralParameters, Topology};

    fn make_atom() -> Atom {
        Atom { name: "C".into(), residue_nr: 1, residue_name: "TEST".into(),
               iac: 0, mass: 12.0, charge: 0.0, is_perturbed: false,
               is_polarisable: false, is_coarse_grained: false }
    }

    #[test]
    fn test_crossdihedral_simple() {
        let mut topo = Topology::new();
        for _ in 0..8 { topo.moltypes[0].atoms.push(make_atom()); }
        topo.mass = vec![12.0; 8];
        topo.inverse_mass = vec![1.0 / 12.0; 8];
        topo.solute_cross_dihedrals().push(CrossDihedral {
            a: 0, b: 1, c: 2, d: 3, e: 4, f: 5, g: 6, h: 7,
            cross_dihedral_type: 0,
        });
        topo.dihedral_parameters.push(DihedralParameters { k: 5.0, cospd: 1.0, pd: 0.0, m: 1 });

        let mut conf = Configuration::new(8, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(-1.0,  0.5, 0.0);
        conf.current_mut().pos[1] = Vec3::new( 0.0,  0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new( 1.0,  0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new( 2.0, -0.5, 0.0);
        conf.current_mut().pos[4] = Vec3::new( 3.0,  0.5, 0.0);
        conf.current_mut().pos[5] = Vec3::new( 4.0,  0.0, 0.0);
        conf.current_mut().pos[6] = Vec3::new( 5.0,  0.0, 0.0);
        conf.current_mut().pos[7] = Vec3::new( 6.0, -0.5, 0.0);

        let result = calculate_crossdihedral_forces(&topo, &conf);
        assert!(!result.energy.is_nan() && !result.energy.is_infinite());
        assert!(result.energy >= 0.0 && result.energy <= 20.0);
    }

    #[test]
    fn test_crossdihedral_force_conservation() {
        let mut topo = Topology::new();
        for _ in 0..8 { topo.moltypes[0].atoms.push(make_atom()); }
        topo.mass = vec![12.0; 8];
        topo.inverse_mass = vec![1.0 / 12.0; 8];
        topo.solute_cross_dihedrals().push(CrossDihedral {
            a: 0, b: 1, c: 2, d: 3, e: 4, f: 5, g: 6, h: 7,
            cross_dihedral_type: 0,
        });
        topo.dihedral_parameters.push(DihedralParameters {
            k: 10.0, cospd: 0.0, pd: std::f64::consts::PI / 2.0, m: 2,
        });

        let mut conf = Configuration::new(8, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(-1.2,  0.3,  0.1);
        conf.current_mut().pos[1] = Vec3::new(-0.4, -0.1,  0.2);
        conf.current_mut().pos[2] = Vec3::new( 0.5,  0.2, -0.1);
        conf.current_mut().pos[3] = Vec3::new( 1.3, -0.2,  0.3);
        conf.current_mut().pos[4] = Vec3::new( 2.1,  0.4, -0.2);
        conf.current_mut().pos[5] = Vec3::new( 3.0, -0.1,  0.1);
        conf.current_mut().pos[6] = Vec3::new( 3.8,  0.3, -0.3);
        conf.current_mut().pos[7] = Vec3::new( 4.7, -0.3,  0.2);

        let result = calculate_crossdihedral_forces(&topo, &conf);
        let total: Vec3 = result.forces.iter().copied().sum();
        assert!(total.length() < 1e-4, "forces not conserved: {}", total.length());
    }

    #[test]
    fn test_crossdihedral_coupling() {
        let mut topo = Topology::new();
        for _ in 0..8 { topo.moltypes[0].atoms.push(make_atom()); }
        topo.mass = vec![12.0; 8];
        topo.inverse_mass = vec![1.0 / 12.0; 8];
        topo.solute_cross_dihedrals().push(CrossDihedral {
            a: 0, b: 1, c: 2, d: 3, e: 4, f: 5, g: 6, h: 7,
            cross_dihedral_type: 0,
        });
        topo.dihedral_parameters.push(DihedralParameters { k: 10.0, cospd: 1.0, pd: 0.0, m: 1 });

        let mut conf = Configuration::new(8, 1, 1);
        conf.current_mut().pos[0] = Vec3::new(-1.0,  0.5, 0.0);
        conf.current_mut().pos[1] = Vec3::new( 0.0,  0.0, 0.0);
        conf.current_mut().pos[2] = Vec3::new( 1.0,  0.0, 0.0);
        conf.current_mut().pos[3] = Vec3::new( 2.0, -0.5, 0.0);
        conf.current_mut().pos[4] = Vec3::new( 3.0,  0.5, 0.0);
        conf.current_mut().pos[5] = Vec3::new( 4.0,  0.0, 0.0);
        conf.current_mut().pos[6] = Vec3::new( 5.0,  0.0, 0.0);
        conf.current_mut().pos[7] = Vec3::new( 6.0, -0.5, 0.0);

        let result = calculate_crossdihedral_forces(&topo, &conf);
        // φ ≈ π, ψ ≈ π → φ+ψ ≈ 2π, cos(2π)=1, E = K*(1+1) = 20
        assert!(result.energy > 15.0 && result.energy < 21.0,
                "Expected ~20, got {}", result.energy);
    }
}
