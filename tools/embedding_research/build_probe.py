#!/usr/bin/env python3
"""Build a read-only embedding probe against the pinned RDKit release source.

The generated translation unit contains upstream Embedder.cpp plus logging hooks.
It links to an existing RDKit SDK; it never replaces the installed RDKit library.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

COMMIT = 'fece8caa860bdf6c9c82bdb9253f22a0cf87150c'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--sdk', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--compiler', default='c++')
    args = parser.parse_args()
    source, sdk, output = (p.resolve() for p in (args.source, args.sdk, args.output))
    revision = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    if revision != COMMIT:
        raise ValueError(f'expected RDKit {COMMIT}, found {revision}')
    upstream = source/'Code/GraphMol/DistGeomHelpers/Embedder.cpp'
    original = upstream.read_text()
    if subprocess.check_output(['git', '-C', str(source), 'diff', 'HEAD', '--',
                                'Code/GraphMol/DistGeomHelpers'], text=True):
        raise ValueError('upstream embedding source must be unmodified')
    decl = '''namespace probe {
void setup(const detail::EmbedArgs &);
void snapshot(const RDGeom::PointPtrVect &, const detail::EmbedArgs &, const char *);
inline void snapshot(const RDGeom::PointPtrVect *p, const detail::EmbedArgs &e, const char *s) {
  snapshot(*p, e, s);
}
}
'''
    anchor = 'constexpr std::array<unsigned int, 3> ALLOWED_ET_VERSIONS'
    assert original.count(anchor) == 1
    instrumented = original.replace(anchor, decl + anchor)
    pattern = r'embedParams\s*\.failures\[EmbedFailureCauses::(\w+)\]\+\+;'
    instrumented, count = re.subn(pattern, lambda m: m[0] +
        f'\nprobe::snapshot(positions, eargs, "{m[1]}");', instrumented)
    assert count == 20, count
    anchor = '    if (numThreads == 1) {\n      detail::embedHelper_'
    assert instrumented.count(anchor) == 1
    instrumented = instrumented.replace(anchor, '    probe::setup(eargs);\n' + anchor)
    anchor = '    // copy the coordinates into the correct conformer'
    assert instrumented.count(anchor) == 1
    instrumented = instrumented.replace(anchor,
        '    if (gotCoords) probe::snapshot(positions, *eargs, "ACCEPTED");\n' + anchor)
    driver = Path(__file__).with_name('probe_driver.cpp').resolve()
    instrumented += '\n#include ' + json.dumps(str(driver)) + '\n'
    output.mkdir(parents=True, exist_ok=True)
    generated = output/'Embedder.instrumented.cpp'
    generated.write_text(instrumented)
    libraries = ['DistGeomHelpers', 'SmilesParse', 'MolAlign', 'ForceFieldHelpers',
                 'SubstructMatch', 'GraphMol', 'DistGeometry', 'Alignment',
                 'ForceField', 'RDGeometryLib', 'RDGeneral']
    command = [args.compiler, '-std=c++20', '-O2', '-pthread',
        '-I'+str(upstream.parent), '-I'+str(sdk/'include/rdkit'),
        '-I'+str(sdk/'include'), '-I'+str(sdk/'include/eigen3'), str(generated),
        '-L'+str(sdk/'lib'), '-Wl,-rpath,'+str(sdk/'lib'),
        *['-lRDKit'+lib for lib in libraries], '-o', str(output/'embedding_probe')]
    subprocess.run(command, check=True)
    manifest = {'rdkit_source_commit': revision,
        'upstream_embedder_sha256': hashlib.sha256(upstream.read_bytes()).hexdigest(),
        'driver_sha256': hashlib.sha256(driver.read_bytes()).hexdigest(),
        'builder_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'generated_source_sha256': hashlib.sha256(generated.read_bytes()).hexdigest(),
        'binary_sha256': hashlib.sha256((output/'embedding_probe').read_bytes()).hexdigest(),
        'command': command, 'failure_hooks': count,
        'compiler': subprocess.check_output([args.compiler, '--version'], text=True)}
    (output/'build-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__':
    main()
