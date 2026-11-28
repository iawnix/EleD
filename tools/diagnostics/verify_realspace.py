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
            'rdkit_evalRealSpaceDerivatives': (
                c.c_int,
                [h, h, c.c_int, dp, c.c_double, dp, dp, dp, dp, dp, dp, dp, dp, dp],
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
        assert self.rdkit_getRealSpaceAPIVersion() == 2
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

    def eval_derivs(self, res, points, *, occ_thr=1e-7):
        npts = len(points)
        pts = (c.c_double * (3 * npts))(*(v for row in points for v in row))
        rho = (c.c_double * npts)(*([12345.0] * npts))
        grad = (c.c_double * (3 * npts))(*([12345.0] * (3 * npts)))
        hess = (c.c_double * (6 * npts))(*([12345.0] * (6 * npts)))
        lap = (c.c_double * npts)(*([12345.0] * npts))
        tau = (c.c_double * npts)(*([12345.0] * npts))
        elf = (c.c_double * npts)(*([12345.0] * npts))
        lol = (c.c_double * npts)(*([12345.0] * npts))
        rdg = (c.c_double * npts)(*([12345.0] * npts))
        sl2 = (c.c_double * npts)(*([12345.0] * npts))
        ok = self.rdkit_evalRealSpaceDerivatives(
            self.env, res, npts, pts, occ_thr, rho, grad, hess, lap, tau, elf, lol, rdg, sl2
        )
        return (
            ok,
            list(rho),
            [tuple(grad[3 * i : 3 * i + 3]) for i in range(npts)],
            [tuple(hess[6 * i : 6 * i + 6]) for i in range(npts)],
            list(lap),
            list(tau),
            list(elf),
            list(lol),
            list(rdg),
            list(sl2),
        )


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

                # ABI v2 derivatives verification (grad, Hessian, Laplacian, tau, ELF, LOL, RDG, sign(lambda2)*rho)
                ok_d, d_rho, grad, hess, lap, tau, elf, lol, rdg, sl2 = api.eval_derivs(copied, probes)
                assert ok_d == 1 and not api.error()
                for a, b in zip(rho, d_rho):
                    assert abs(a - b) < 1e-12, (method, name, a, b)
                for i in range(len(probes)):
                    assert abs((hess[i][0] + hess[i][2] + hess[i][5]) - lap[i]) < 1e-12
                    assert tau[i] >= 0.0 and 0.0 <= elf[i] <= 1.0 + 1e-12 and 0.0 <= lol[i] <= 1.0 + 1e-12
                    assert rdg[i] >= 0.0 and abs(abs(sl2[i]) - rho[i]) < 1e-12

                # Analytical vs finite-difference gradient and Hessian check on H2S (includes d-shells)
                if name == 'hydrogen_sulfide' and method == 2:
                    h_step = 1e-4
                    p0 = probes[1]
                    fd_pts = []
                    for ax in range(3):
                        for sgn in (1.0, -1.0):
                            q = list(p0); q[ax] += sgn * h_step; fd_pts.append(tuple(q))
                    ok_fd, rho_fd, _, _, _ = api.eval_points(copied, fd_pts)
                    assert ok_fd == 1 and not api.error()
                    for ax in range(3):
                        g_num = (rho_fd[2 * ax] - rho_fd[2 * ax + 1]) / (2.0 * h_step)
                        assert abs(grad[1][ax] - g_num) < 1e-7, (ax, grad[1][ax], g_num)
                    # Finite-difference check on Hessian from analytical gradient
                    ok_hfd, _, grad_fd, _, _, _, _, _, _, _ = api.eval_derivs(copied, fd_pts)
                    assert ok_hfd == 1 and not api.error()
                    h_idx = {(0, 0): 0, (0, 1): 1, (1, 1): 2, (0, 2): 3, (1, 2): 4, (2, 2): 5}
                    for (a_i, b_i), slot in h_idx.items():
                        h_num = (grad_fd[2 * b_i][a_i] - grad_fd[2 * b_i + 1][a_i]) / (2.0 * h_step)
                        assert abs(hess[1][slot] - h_num) < 1e-6, ((a_i, b_i), hess[1][slot], h_num)
                    # Rotation invariance of v2 scalar fields on H2S
                    ok_dr, d_rho_r, grad_r, _, lap_r, tau_r, elf_r, lol_r, rdg_r, sl2_r = api.eval_derivs(
                        res_rot, moved_probes
                    )
                    assert ok_dr == 1 and not api.error()
                    for i in range(len(probes)):
                        gn0 = math.sqrt(sum(v * v for v in grad[i]))
                        gn1 = math.sqrt(sum(v * v for v in grad_r[i]))
                        assert abs(gn0 - gn1) < 1e-7
                        assert abs(lap[i] - lap_r[i]) < 1e-7
                        assert abs(tau[i] - tau_r[i]) < 1e-7
                        assert abs(elf[i] - elf_r[i]) < 1e-7
                        assert abs(lol[i] - lol_r[i]) < 1e-7
                        assert abs(rdg[i] - rdg_r[i]) < 1e-6
                        assert abs(sl2[i] - sl2_r[i]) < 1e-8

                # Bad-input guard leaves v2 buffers untouched
                bad_d = api.eval_derivs(copied, probes, occ_thr=-0.5)
                assert bad_d[0] == 0 and api.error()
                assert bad_d[1] == [12345.0] * len(probes)

                api.singlepoint(api.env, None, calc, res)
                assert api.error()
                assert api.eval_points(res, probes)[0] == 0 and api.error()
                assert api.eval_derivs(res, probes)[0] == 0 and api.error()
                assert api.eval_points(copied, probes) == (1, rho, alie, vrho, lea)

                rows.append({'method': method, 'case': name, 'rho_0_au': rho[0], 'alie_0_ev': alie[0], 'elf_0': elf[0]})
            finally:
                api.close()

    # Exact closed-shell 2-electron H2 check: Becke-Edgecombe ELF == 1.0 everywhere!
    api = API(prefix, library)
    try:
        _, _, res_h2 = api.calculate([1, 1], [(0.0, 0.0, 0.0), (0.0, 0.0, 1.4)], 0, 0, 2)
        ok_h2, _, _, _, _, _, elf_h2, _, _, _ = api.eval_derivs(res_h2, [(0.0, 0.0, 0.7), (0.5, 0.3, 0.7)])
        assert ok_h2 == 1 and not api.error()
        assert all(abs(v - 1.0) < 1e-10 for v in elf_h2), elf_h2
    finally:
        api.close()

    return {'api_version': 2, 'cases': rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prefix', type=Path)
    parser.add_argument('--library', type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.prefix, args.library), indent=2))


if __name__ == '__main__':
    main()
