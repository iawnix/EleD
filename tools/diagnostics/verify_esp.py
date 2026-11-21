#!/usr/bin/env python3
"""Verify the installed ESP C ABI against CLI --esp; standard library only."""
import argparse
import ctypes as c
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

BOHR_ANGSTROM = 0.529177210903
POINTS = [(4.1, 2.3, 3.7), (-3.2, 4.5, 2.1), (2.2, -4.3, 3.8), (0.4, 1.2, -5.1)]
CASES = [
    ('water', [8, 1, 1], [(0., 0., 0.), (0., 0., 1.8), (1.7, 0., -.6)], 0, 0),
    ('hydrogen_sulfide', [16, 1, 1], [(0., 0., 0.), (0., 0., 2.5), (2.4, 0., -.6)], 0, 0),
    ('ammonium', [7, 1, 1, 1, 1], [(0., 0., 0.), (1.1, 1.1, 1.1),
                                    (-1.1, -1.1, 1.1), (-1.1, 1.1, -1.1), (1.1, -1.1, -1.1)], 1, 0),
    ('hydroxyl', [8, 1], [(0., 0., 0.), (0., 0., 1.8)], 0, 1),
]


class API:
    def __init__(self, prefix, library=None):
        self.prefix = Path(prefix)
        self.lib = c.CDLL(str(library or self.prefix/'lib/libxtb.so'))
        self.handles = []
        h, dp, ip = c.c_void_p, c.POINTER(c.c_double), c.POINTER(c.c_int)
        signatures = {
            'newEnvironment': (h, []), 'setVerbosity': (None, [h, c.c_int]),
            'checkEnvironment': (c.c_int, [h]), 'getError': (None, [h, c.c_char_p, ip]),
            'newMolecule': (h, [h, ip, ip, dp, dp, ip, h, h]),
            'newCalculator': (h, []), 'newResults': (h, []), 'copyResults': (h, [h]),
            'singlepoint': (None, [h, h, h, h]), 'setMaxIter': (None, [h, h, c.c_int]),
            'getEnergy': (None, [h, h, dp]),
            'rdkit_getESPAPIVersion': (c.c_int, []),
            'rdkit_setKeepESP': (None, [h, h, c.c_bool]),
            'rdkit_getESP': (c.c_int, [h, h, c.c_int, dp, c.c_int, dp]),
        }
        for method in (0, 1, 2):
            signatures[f'loadGFN{method}xTB'] = (None, [h, h, h, c.c_char_p])
        for kind in ('Environment', 'Molecule', 'Calculator', 'Results'):
            signatures['del'+kind] = (None, [c.POINTER(h)])
        for name, (result, args) in signatures.items():
            f = getattr(self.lib, 'xtb_'+name)
            f.restype, f.argtypes = result, args
            setattr(self, name, f)
        assert self.rdkit_getESPAPIVersion() == 1
        self.env = self.own('Environment', self.newEnvironment())
        self.setVerbosity(self.env, 0)

    def own(self, kind, value):
        assert value
        self.handles.append((kind, value))
        return value

    def close(self):
        for kind, value in reversed(self.handles):
            getattr(self, 'del'+kind)(c.byref(c.c_void_p(value)))
        self.handles.clear()

    def error(self):
        if not self.checkEnvironment(self.env):
            return ''
        buf, size = c.create_string_buffer(4096), c.c_int(4096)
        self.getError(self.env, buf, c.byref(size))
        return buf.value.decode()

    def calculate(self, case, method, keep=True):
        _, numbers, coords, charge, spin = case
        n, q, uhf = c.c_int(len(numbers)), c.c_double(charge), c.c_int(spin)
        at = (c.c_int*len(numbers))(*numbers)
        xyz = (c.c_double*(3*len(numbers)))(*(v for row in coords for v in row))
        mol = self.own('Molecule', self.newMolecule(self.env, c.byref(n), at, xyz,
                                                 c.byref(q), c.byref(uhf), None, None))
        calc = self.own('Calculator', self.newCalculator())
        res = self.own('Results', self.newResults())
        getattr(self, f'loadGFN{method}xTB')(self.env, mol, calc,
            str(self.prefix/f'share/xtb/param_gfn{method}-xtb.txt').encode())
        self.rdkit_setKeepESP(self.env, res, keep)
        self.singlepoint(self.env, mol, calc, res)
        assert not (error := self.error()), error
        return mol, calc, res

    def read(self, res, points=POINTS, capacity=None, count=None, null_points=False, null_output=False):
        xyz = (c.c_double*(3*len(points)))(*(v for row in points for v in row))
        output = (c.c_double*len(points))(*([12345.]*len(points)))
        ok = self.rdkit_getESP(self.env, res, len(points) if count is None else count,
                              None if null_points else xyz,
                              len(points) if capacity is None else capacity,
                              None if null_output else output)
        return ok, list(output)

    def values(self, res, points=POINTS):
        ok, values = self.read(res, points)
        assert ok == 1, self.error()
        assert all(map(math.isfinite, values))
        return values


def difference(a, b):
    return max(abs(x-y) for x, y in zip(a, b))


def cli_values(prefix, executable, case, method, cwd):
    name, numbers, coords, charge, spin = case
    symbols = {1: 'H', 7: 'N', 8: 'O', 16: 'S'}
    xyz = '\n'.join([str(len(numbers)), name] + [
        f'{symbols[z]} ' + ' '.join(f'{v*BOHR_ANGSTROM:.15g}' for v in row)
        for z, row in zip(numbers, coords)])+'\n'
    (cwd/'molecule.xyz').write_text(xyz)
    (cwd/'esp_coord').write_text('\n'.join(' '.join(map(str, row)) for row in POINTS)+'\n')
    env = os.environ.copy()
    env.update(XTBPATH=str(prefix/'share/xtb'), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    run = subprocess.run([str(executable), 'molecule.xyz', '--gfn', str(method), '--chrg', str(charge),
                          '--uhf', str(spin), '--esp'], cwd=cwd, env=env,
                         capture_output=True, text=True, timeout=60)
    assert run.returncode == 0, run.stdout+run.stderr
    return [float(line.split()[3]) for line in (cwd/'xtb_esp.dat').read_text().splitlines()]


def verify(prefix, library=None, executable=None):
    prefix = Path(prefix).resolve()
    executable = Path(executable).resolve() if executable else prefix/'bin/xtb'
    rows = []
    with tempfile.TemporaryDirectory(prefix='eled-esp-') as temporary:
        work = Path(temporary)
        old_cwd = Path.cwd()
        os.chdir(work)
        try:
            for method in (0, 1, 2):
                for case in CASES:
                    api = API(prefix, library)
                    try:
                        mol, calc, res = api.calculate(case, method)
                        before = set(work.iterdir())
                        values = api.values(res)
                        assert set(work.iterdir()) == before, 'ESP API wrote files'
                        reference_dir = work/f'{method}-{case[0]}'
                        reference_dir.mkdir()
                        reference = cli_values(prefix, executable, case, method, reference_dir)
                        error = difference(values, reference)
                        assert error < 2e-7, (method, case[0], error, values, reference)
                        copied = api.own('Results', api.copyResults(res))
                        assert api.values(copied) == values
                        # A different molecule/calculator must not alter a retained snapshot.
                        api.calculate(CASES[(CASES.index(case)+1) % len(CASES)], (method+1) % 3)
                        assert api.values(copied) == values
                        api.rdkit_setKeepESP(api.env, res, False)
                        ok, untouched = api.read(res)
                        assert ok == 0 and untouched == [12345.]*len(POINTS)
                        assert api.error()
                        assert api.values(copied) == values
                        for kwargs in ({'capacity': 3}, {'count': -1}, {'count': 0},
                                       {'null_points': True}, {'null_output': True},
                                       {'points': [(math.nan, 0., 0.)]},
                                       {'points': [(math.inf, 0., 0.)]}, {'points': [case[2][0]]}):
                            ok, untouched = api.read(copied, **kwargs)
                            assert ok == 0 and all(v == 12345. for v in untouched)
                            assert api.error()
                        # Failed setup invalidates an otherwise valid snapshot.
                        api.singlepoint(api.env, None, calc, copied)
                        assert api.error()
                        assert api.read(copied)[0] == 0
                        assert api.error()
                        # Failed SCC must also invalidate, even when reusing wavefunction storage.
                        if method:
                            api.rdkit_setKeepESP(api.env, res, True)
                            api.setMaxIter(api.env, calc, 1)
                            api.singlepoint(api.env, mol, calc, res)
                            assert api.error()
                            assert api.read(res)[0] == 0
                            assert api.error()
                        rows.append({'method': method, 'case': case[0], 'max_cli_error_au': error})
                    finally:
                        api.close()
            # Rigid motion, atom permutation, and point ordering contracts.
            for method in (0, 1, 2):
                api = API(prefix, library)
                try:
                    case = CASES[0]
                    _, _, res = api.calculate(case, method)
                    values = api.values(res)
                    assert api.values(res, POINTS[::-1]) == values[::-1]
                    def transform(row):
                        x, y, z = row
                        return (.6*x-.8*y+1.3, .8*x+.6*y-.7, z+2.1)
                    moved = ('moved', case[1][::-1], list(map(transform, case[2][::-1])), 0, 0)
                    _, _, other = api.calculate(moved, method)
                    error = difference(values, api.values(other, list(map(transform, POINTS))))
                    # Cartesian integral screening is orientation dependent at CLI precision.
                    assert error < 2e-5, (method, error)
                    rows.append({'method': method, 'rigid_motion_and_permutation_error_au': error})
                finally:
                    api.close()
        finally:
            os.chdir(old_cwd)
    return rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prefix', type=Path)
    parser.add_argument('--library', type=Path)
    parser.add_argument('--executable', type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.prefix, args.library, args.executable)))
