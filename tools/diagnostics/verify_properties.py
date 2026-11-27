#!/usr/bin/env python3
"""Verify the installed extended properties C ABI; standard library only."""
from __future__ import annotations

import argparse
import ctypes as c
import json
import math
from pathlib import Path

CASES = [
    ('water', [8, 1, 1], [(0.0, 0.0, 0.0), (0.0, 1.43233673, 1.10866041), (0.0, -1.43233673, 1.10866041)], 0, 0),
    ('hydrogen_sulfide', [16, 1, 1], [(0.0, 0.0, 0.0), (0.0, 0.0, 2.5), (2.4, 0.0, -0.6)], 0, 0),
    ('ammonium', [7, 1, 1, 1, 1], [
        (0.0, 0.0, 0.0), (1.1, 1.1, 1.1), (-1.1, -1.1, 1.1), (-1.1, 1.1, -1.1), (1.1, -1.1, -1.1),
    ], 1, 0),
    ('hydroxyl', [8, 1], [(0.0, 0.0, 0.0), (0.0, 0.0, 1.8)], 0, 1),
]


class API:
    def __init__(self, prefix, library=None):
        self.prefix = Path(prefix)
        self.lib = c.CDLL(str(library or self.prefix / 'lib/libxtb.so'))
        self.handles = []
        h, dp, ip = c.c_void_p, c.POINTER(c.c_double), c.POINTER(c.c_int)
        signatures = {
            'newEnvironment': (h, []), 'setVerbosity': (None, [h, c.c_int]),
            'checkEnvironment': (c.c_int, [h]), 'getError': (None, [h, c.c_char_p, ip]),
            'newMolecule': (h, [h, ip, ip, dp, dp, ip, h, h]),
            'newCalculator': (h, []), 'newResults': (h, []), 'copyResults': (h, [h]),
            'singlepoint': (None, [h, h, h, h]),
            'setMaxIter': (None, [h, h, c.c_int]),
            'setElectronicTemp': (None, [h, h, c.c_double]),
            'setSolvent': (None, [h, h, c.c_char_p, ip, dp, ip, c.c_char_p]),
            'getEnergy': (None, [h, h, dp]),
            'getSolvationEnergy': (None, [h, h, dp]),
            'getDipole': (None, [h, h, dp]),
            'getCharges': (None, [h, h, dp]),
            'rdkit_setKeepESP': (None, [h, h, c.c_bool]),
            'rdkit_getESP': (c.c_int, [h, h, c.c_int, dp, c.c_int, dp]),
            'rdkit_getPropertiesAPIVersion': (c.c_int, []),
            'rdkit_getEnergyDecomposition': (c.c_int, [h, h, ip, c.c_int, dp]),
            'rdkit_getPolarizability': (c.c_int, [h, h, c.c_int, dp, dp, dp]),
            'rdkit_getMultipoles': (c.c_int, [h, h, c.c_int, dp, dp, dp]),
            'rdkit_getNuclearESP': (c.c_int, [h, h, c.c_int, dp]),
        }
        for method in (0, 1, 2):
            signatures[f'loadGFN{method}xTB'] = (None, [h, h, h, c.c_char_p])
        for kind in ('Environment', 'Molecule', 'Calculator', 'Results'):
            signatures['del' + kind] = (None, [c.POINTER(h)])
        for name, (restype, argtypes) in signatures.items():
            fn = getattr(self.lib, 'xtb_' + name)
            fn.restype, fn.argtypes = restype, argtypes
            setattr(self, name, fn)
        assert self.rdkit_getPropertiesAPIVersion() == 1
        self.env = self.own('Environment', self.newEnvironment())
        self.setVerbosity(self.env, 0)

    def own(self, kind, value):
        assert value
        self.handles.append((kind, value))
        return value

    def close(self):
        for kind, value in reversed(self.handles):
            getattr(self, 'del' + kind)(c.byref(c.c_void_p(value)))
        self.handles.clear()

    def error(self):
        if not self.checkEnvironment(self.env):
            return ''
        buf, size = c.create_string_buffer(4096), c.c_int(4096)
        self.getError(self.env, buf, c.byref(size))
        return buf.value.decode()

    def calculate(self, case, method, *, etemp=300.0, solvent=None, solv_model='alpb', keep_esp=False):
        _, numbers, coords, charge, spin = case
        n = c.c_int(len(numbers))
        q = c.c_double(charge)
        uhf = c.c_int(spin)
        at = (c.c_int * len(numbers))(*numbers)
        xyz = (c.c_double * (3 * len(numbers)))(*(v for row in coords for v in row))
        mol = self.own('Molecule', self.newMolecule(
            self.env, c.byref(n), at, xyz, c.byref(q), c.byref(uhf), None, None))
        calc = self.own('Calculator', self.newCalculator())
        res = self.own('Results', self.newResults())
        getattr(self, f'loadGFN{method}xTB')(
            self.env, mol, calc, str(self.prefix / f'share/xtb/param_gfn{method}-xtb.txt').encode())
        self.setElectronicTemp(self.env, calc, etemp)
        if solvent is not None:
            self.setSolvent(self.env, calc, solvent.encode(), None, None, None, solv_model.encode())
        if keep_esp:
            self.rdkit_setKeepESP(self.env, res, True)
        self.singlepoint(self.env, mol, calc, res)
        assert not (err := self.error()), err
        return mol, calc, res

    def read_energy_decomp(self, res, capacity=15):
        flags = c.c_int(-999)
        vals = (c.c_double * 15)(*([12345.0] * 15))
        ok = self.rdkit_getEnergyDecomposition(self.env, res, c.byref(flags), capacity, vals)
        return ok, flags.value, list(vals)

    def read_polarizability(self, res, natoms):
        mol_v = (c.c_double * 3)(*([12345.0] * 3))
        alpha = (c.c_double * natoms)(*([12345.0] * natoms))
        c6 = (c.c_double * natoms)(*([12345.0] * natoms))
        ok = self.rdkit_getPolarizability(self.env, res, natoms, mol_v, alpha, c6)
        return ok, list(mol_v), list(alpha), list(c6)

    def read_multipoles(self, res, natoms):
        dip = (c.c_double * (3 * natoms))(*([12345.0] * (3 * natoms)))
        qp = (c.c_double * (6 * natoms))(*([12345.0] * (6 * natoms)))
        mqp = (c.c_double * 6)(*([12345.0] * 6))
        ok = self.rdkit_getMultipoles(self.env, res, natoms, dip, qp, mqp)
        return ok, list(dip), list(qp), list(mqp)

    def read_nuclear_esp(self, res, natoms):
        out = (c.c_double * natoms)(*([12345.0] * natoms))
        ok = self.rdkit_getNuclearESP(self.env, res, natoms, out)
        return ok, list(out)


def verify(prefix, library=None):
    prefix = Path(prefix).resolve()
    rows = []
    for method in (0, 1, 2):
        for case in CASES:
            name, numbers, coords, _, spin = case
            natoms = len(numbers)
            api = API(prefix, library)
            try:
                mol, calc, res = api.calculate(case, method, keep_esp=(name == 'water'))
                e_ref = c.c_double()
                dip_ref = (c.c_double * 3)()
                q_ref = (c.c_double * natoms)()
                api.getEnergy(api.env, res, c.byref(e_ref))
                api.getDipole(api.env, res, dip_ref)
                api.getCharges(api.env, res, q_ref)
                assert not api.error()

                ok, flags, evals = api.read_energy_decomp(res)
                assert ok == 1 and not api.error()
                assert all(map(math.isfinite, evals))
                assert abs(evals[0] - e_ref.value) < 1e-12
                assert (flags & 256) == 0  # gas phase
                assert (flags & 512) != 0  # Fermi at 300 K
                if spin == 0:
                    assert abs(evals[9] - evals[10]) < 1e-12
                if method == 2:
                    assert (flags & (16 | 32)) == (16 | 32)
                    assert (flags & 128) == 0
                    assert abs((evals[1] + evals[2]) - evals[0]) < 1e-10
                elif method == 1:
                    assert (flags & (16 | 32)) == 0
                    assert (flags & 128) == 128
                    assert abs((evals[1] + evals[2] + evals[6] + evals[7]) - evals[0]) < 1e-10
                else:
                    assert (flags & (16 | 32)) == 0
                    assert (flags & 128) == 128
                    assert abs((evals[1] + evals[2] + evals[3] + evals[6] + evals[7]) - evals[0]) < 1e-10

                # Polarizability
                pol_ok, mol_v, alpha, c6 = api.read_polarizability(res, natoms)
                assert not api.error()
                if method in (0, 2):
                    assert pol_ok == 1
                    assert mol_v[0] > 0.0 and mol_v[1] > 0.0 and mol_v[2] > 0.0
                    assert all(a > 0.0 for a in alpha) and all(v > 0.0 for v in c6)
                    assert abs(sum(alpha) - mol_v[0]) < 1e-10
                else:
                    assert pol_ok == 0
                    assert mol_v == [12345.0] * 3 and alpha == [12345.0] * natoms and c6 == [12345.0] * natoms

                # Multipoles
                mp_ok, dip_a, qp_a, mqp = api.read_multipoles(res, natoms)
                assert not api.error()
                if method in (1, 2):
                    assert mp_ok == 1
                    assert abs(mqp[0] + mqp[2] + mqp[5]) < 1e-12
                    for i in range(natoms):
                        qpi = qp_a[6 * i: 6 * (i + 1)]
                        assert abs(qpi[0] + qpi[2] + qpi[5]) < 1e-12
                    # Dipole identity: sum_i (q_i * r_i + dipm_i) == total dipole
                    for axis in range(3):
                        recon = sum(q_ref[i] * coords[i][axis] + dip_a[3 * i + axis] for i in range(natoms))
                        assert abs(recon - dip_ref[axis]) < 1e-10, (method, name, axis, recon, dip_ref[axis])
                else:
                    assert mp_ok == 0
                    assert mqp == [12345.0] * 6

                # Nuclear ESP
                epn_ok, epn = api.read_nuclear_esp(res, natoms)
                assert epn_ok == 1 and not api.error()
                assert all(map(math.isfinite, epn))
                if name == 'water':
                    assert epn[0] < epn[1] < 0.0
                    assert abs(epn[1] - epn[2]) < 1e-10
                    # Limit check against esp_points around atom 0
                    delta = 1e-4
                    probes = [(delta, 0, 0), (-delta, 0, 0), (0, delta, 0), (0, -delta, 0), (0, 0, delta), (0, 0, -delta)]
                    pxyz = (c.c_double * (3 * len(probes)))(*(v for row in probes for v in row))
                    pout = (c.c_double * len(probes))()
                    assert api.rdkit_getESP(api.env, res, len(probes), pxyz, len(probes), pout) == 1
                    z_ox = 6.0  # valence core charge of oxygen in xTB
                    lim = sum(v - z_ox / delta for v in pout) / len(probes)
                    assert abs(lim - epn[0]) < 1e-3, (method, lim, epn[0])

                # Copy ownership & invalidation tests
                copied = api.own('Results', api.copyResults(res))
                assert api.read_energy_decomp(copied) == (ok, flags, evals)
                assert api.read_nuclear_esp(copied, natoms) == (1, epn)
                # Wrong capacity / atom count leaves buffers untouched
                bad_ok, bad_f, bad_v = api.read_energy_decomp(copied, capacity=10)
                assert bad_ok == 0 and bad_f == -999 and bad_v == [12345.0] * 15
                bad_epn_ok, bad_epn = api.read_nuclear_esp(copied, natoms + 1)
                assert bad_epn_ok == 0 and bad_epn == [12345.0] * (natoms + 1) and api.error()
                # Failed setup invalidates res but leaves copied intact
                api.singlepoint(api.env, None, calc, res)
                assert api.error()
                assert api.read_energy_decomp(res)[0] == 0
                assert api.read_nuclear_esp(res, natoms)[0] == 0 and api.error()
                assert api.read_energy_decomp(copied) == (ok, flags, evals)
                assert api.read_nuclear_esp(copied, natoms) == (1, epn)

                rows.append({
                    'method': method, 'case': name,
                    'energy_flags': flags, 'molpol_bohr3': mol_v[0] if pol_ok else None,
                    'epn_0_au': epn[0],
                })
            finally:
                api.close()

    # Solvation & zero-etemp checks
    for method in (1, 2):
        api = API(prefix, library)
        try:
            _, _, res_solv = api.calculate(CASES[0], method, solvent='water', solv_model='alpb')
            ok, flags, evals = api.read_energy_decomp(res_solv)
            assert ok == 1 and (flags & 256) == 256 and (flags & 1024) == 1024
            g_api = c.c_double()
            api.getSolvationEnergy(api.env, res_solv, c.byref(g_api))
            assert not api.error()
            assert abs(evals[8] - g_api.value) < 1e-12
            assert abs(evals[8] - sum(evals[11:15])) < 1e-10
            assert evals[8] < 0.0  # water in water is stabilized

            _, _, res_zero_t = api.calculate(CASES[0], method, etemp=0.0)
            ok0, flags0, _ = api.read_energy_decomp(res_zero_t)
            assert ok0 == 1 and (flags0 & 512) == 0  # Fermi smearing skipped at 0 K
            rows.append({'method': method, 'case': 'water_alpb_water', 'g_solv_ha': evals[8]})
        finally:
            api.close()

    # Rigid motion + atom permutation invariance
    for method in (0, 1, 2):
        api = API(prefix, library)
        try:
            case = CASES[0]
            _, _, res1 = api.calculate(case, method)
            _, _, evals1 = api.read_energy_decomp(res1)
            _, epn1 = api.read_nuclear_esp(res1, 3)
            p_ok1, mol_v1, alpha1, c61 = api.read_polarizability(res1, 3)
            m_ok1, dip1, qp1, _ = api.read_multipoles(res1, 3)

            def transform(row):
                x, y, z = row
                return (0.6 * x - 0.8 * y + 1.3, 0.8 * x + 0.6 * y - 0.7, z + 2.1)

            perm = [2, 0, 1]
            moved = ('moved', [case[1][i] for i in perm], [transform(case[2][i]) for i in perm], 0, 0)
            _, _, res2 = api.calculate(moved, method)
            _, _, evals2 = api.read_energy_decomp(res2)
            _, epn2 = api.read_nuclear_esp(res2, 3)
            p_ok2, mol_v2, alpha2, c62 = api.read_polarizability(res2, 3)
            m_ok2, dip2, qp2, _ = api.read_multipoles(res2, 3)

            assert max(abs(a - b) for a, b in zip(evals1, evals2)) < 1e-8
            assert max(abs(epn1[perm[i]] - epn2[i]) for i in range(3)) < 1e-7
            if p_ok1:
                assert p_ok2
                assert max(abs(a - b) for a, b in zip(mol_v1, mol_v2)) < 1e-8
                assert max(abs(alpha1[perm[i]] - alpha2[i]) for i in range(3)) < 1e-8
                assert max(abs(c61[perm[i]] - c62[i]) for i in range(3)) < 1e-8
            if m_ok1:
                assert m_ok2
                for i in range(3):
                    orig = perm[i]
                    d1 = math.sqrt(sum(v * v for v in dip1[3 * orig: 3 * (orig + 1)]))
                    d2 = math.sqrt(sum(v * v for v in dip2[3 * i: 3 * (i + 1)]))
                    assert abs(d1 - d2) < 1e-7
                    def qnorm(q):
                        return math.sqrt(q[0]**2 + q[2]**2 + q[5]**2 + 2 * (q[1]**2 + q[3]**2 + q[4]**2))
                    assert abs(qnorm(qp1[6 * orig: 6 * (orig + 1)]) - qnorm(qp2[6 * i: 6 * (i + 1)])) < 1e-7
            rows.append({'method': method, 'invariance_check': 'passed'})
        finally:
            api.close()
    return rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prefix', type=Path)
    parser.add_argument('--library', type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.prefix, args.library)))
