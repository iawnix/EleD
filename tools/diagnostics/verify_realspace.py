#!/usr/bin/env python3
"""Verify the installed real-space wavefunction C ABI; standard library only."""
from __future__ import annotations

import argparse
import ctypes as c
import json
import math
from pathlib import Path

CASES = [
    ('water', [8, 1, 1], [(0.0, 0.0, 0.0), (0.0, 1.43233673, 1.10866041), (0.0, -1.43233673, 1.10866041)], 0, 0, 8.0),
    ('hydrogen_sulfide', [16, 1, 1], [(0.0, 0.0, 0.0), (0.0, 0.0, 2.5), (2.4, 0.0, -0.6)], 0, 0, 8.0),
    ('ammonium', [7, 1, 1, 1, 1], [
        (0.0, 0.0, 0.0), (1.1, 1.1, 1.1), (-1.1, -1.1, 1.1), (-1.1, 1.1, -1.1), (1.1, -1.1, -1.1),
    ], 1, 0, 8.0),
    ('hydroxyl', [8, 1], [(0.0, 0.0, 0.0), (0.0, 0.0, 1.8)], 0, 1, 7.0),
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
            'getDipole': (None, [h, h, dp]),
            'rdkit_getRealSpaceAPIVersion': (c.c_int, []),
            'rdkit_evalRealSpace': (
                c.c_int,
                [h, h, c.c_int, dp, c.c_double, c.c_double, dp, dp, dp, dp],
            ),
        }
        for method in (0, 1, 2):
            signatures[f'loadGFN{method}xTB'] = (None, [h, h, h, c.c_char_p])
        for kind in ('Environment', 'Molecule', 'Calculator', 'Results'):
            signatures['del' + kind] = (None, [c.POINTER(h)])
        for name, (restype, argtypes) in signatures.items():
            fn = getattr(self.lib, 'xtb_' + name)
            fn.restype, fn.argtypes = restype, argtypes
            setattr(self, name, fn)
        assert self.rdkit_getRealSpaceAPIVersion() == 1
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

    def calculate(self, numbers, coords, charge, spin, method):
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
        self.singlepoint(self.env, mol, calc, res)
        assert not (err := self.error()), err
        return mol, calc, res

    def eval_points(self, res, points, *, occ_thr=1e-7, lea_win=0.0):
        npts = len(points)
        pts = (c.c_double * (3 * npts))(*(v for row in points for v in row))
        rho = (c.c_double * npts)(*([12345.0] * npts))
        alie = (c.c_double * npts)(*([12345.0] * npts))
        vrho = (c.c_double * npts)(*([12345.0] * npts))
        lea = (c.c_double * npts)(*([12345.0] * npts))
        ok = self.rdkit_evalRealSpace(
            self.env, res, npts, pts, occ_thr, lea_win, rho, alie, vrho, lea
        )
        return ok, list(rho), list(alie), list(vrho), list(lea)


def _rotate_translate(coords):
    theta = 0.71
    ct, st = math.cos(theta), math.sin(theta)
    out = []
    for x, y, z in coords:
        out.append((ct * x - st * y + 1.3, st * x + ct * y - 0.7, z + 2.1))
    return out


def verify(prefix, library=None):
    prefix = Path(prefix).resolve()
    probes = [(0.0, 0.0, 2.0), (1.5, 0.0, 0.5), (0.0, 2.2, 1.5), (0.0, 0.0, -2.2)]
    rows = []
    for method in (0, 1, 2):
        for name, numbers, coords, charge, spin, ne_expected in CASES:
            api = API(prefix, library)
            try:
                _, calc, res = api.calculate(numbers, coords, charge, spin, method)
                ok, rho, alie, vrho, lea = api.eval_points(res, probes)
                assert ok == 1 and not api.error()
                assert all(v >= 0.0 and math.isfinite(v) for v in rho)
                assert all(v >= 0.0 and math.isfinite(v) for v in vrho)
                assert all(math.isfinite(v) for v in alie)
                assert all(math.isfinite(v) for v in lea)
                assert max(rho) > 1e-4

                # Rigid motion invariance (covers s, p, and d shells in H2S)
                moved_coords = _rotate_translate(coords)
                moved_probes = _rotate_translate(probes)
                _, _, res_rot = api.calculate(numbers, moved_coords, charge, spin, method)
                ok_r, rho_r, alie_r, vrho_r, lea_r = api.eval_points(res_rot, moved_probes)
                assert ok_r == 1 and not api.error()
                for a, b in zip(rho, rho_r):
                    assert abs(a - b) < 1e-8, (method, name, a, b)
                for a, b in zip(alie, alie_r):
                    assert abs(a - b) < 1e-7, (method, name, a, b)
                for a, b in zip(vrho, vrho_r):
                    assert abs(a - b) < 1e-8, (method, name, a, b)
                for a, b in zip(lea, lea_r):
                    assert abs(a - b) < 1e-7, (method, name, a, b)

                # Quadrature check on water for GFN2 (verifies normalization and p-sign against dipole)
                if name == 'water' and method == 2:
                    dip = (c.c_double * 3)()
                    api.getDipole(api.env, res, dip)
                    step = 0.3
                    axis = [round(-5.4 + i * step, 10) for i in range(37)]
                    grid = [(x, y, z) for x in axis for y in axis for z in axis]
                    dv = step ** 3
                    ok_q, rho_q, _, _, _ = api.eval_points(res, grid)
                    assert ok_q == 1 and not api.error()
                    ne_num = sum(rho_q) * dv
                    assert abs(ne_num - ne_expected) < 0.01, (ne_num, ne_expected)
                    elec_z = sum(pt[2] * val for pt, val in zip(grid, rho_q)) * dv
                    nuc_z = 6.0 * coords[0][2] + 1.0 * coords[1][2] + 1.0 * coords[2][2]
                    assert abs((nuc_z - elec_z) - dip[2]) < 1e-2, (nuc_z - elec_z, dip[2])

                # Copy ownership & error guards
                copied = api.own('Results', api.copyResults(res))
                assert api.eval_points(copied, probes) == (1, rho, alie, vrho, lea)
                bad_ok, b_rho, b_alie, b_vrho, b_lea = api.eval_points(copied, probes, occ_thr=-1.0)
                assert bad_ok == 0 and api.error()
                assert b_rho == [12345.0] * len(probes) and b_alie == [12345.0] * len(probes)
                api.singlepoint(api.env, None, calc, res)
                assert api.error()
                assert api.eval_points(res, probes)[0] == 0 and api.error()
                assert api.eval_points(copied, probes) == (1, rho, alie, vrho, lea)

                rows.append({'method': method, 'case': name, 'rho_0_au': rho[0], 'alie_0_ev': alie[0]})
            finally:
                api.close()
    return {'api_version': 1, 'cases': rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prefix', type=Path)
    parser.add_argument('--library', type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.prefix, args.library), indent=2))


if __name__ == '__main__':
    main()
