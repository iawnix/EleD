"""Build the optional RDKit-style adapter against the reviewed xTB runtime."""
import os
from pathlib import Path
import subprocess
import sys

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

ROOT = Path(__file__).resolve().parent


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        output = Path(self.get_ext_fullpath(ext.name)).resolve().parent
        build = Path(self.build_temp).resolve() / 'rdxtb'
        prefix = Path(os.environ.get('RDXTB_RDKIT_PREFIX', sys.prefix)).resolve()
        local_xtb = ROOT / '.xtb-native'
        if not local_xtb.is_dir():
            local_xtb = ROOT / '.xtb-rdkit'
        xtb = Path(os.environ.get('RDXTB_XTB_PREFIX',
                                 local_xtb if local_xtb.is_dir() else prefix)).resolve()
        args = ['cmake', '-S', str(ROOT), '-B', str(build), '-G', 'Ninja',
                '-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_PREFIX_PATH=' + str(prefix),
                '-DPython_EXECUTABLE=' + sys.executable,
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + str(output),
                '-DRDXTB_XTB_PREFIX=' + str(xtb),
                '-DCMAKE_BUILD_RPATH=' + str(prefix / 'lib'),
                '-DBUILD_TESTING=OFF']
        subprocess.run(args, check=True)
        subprocess.run(['cmake', '--build', str(build), '--parallel',
                        os.environ.get('CMAKE_BUILD_PARALLEL_LEVEL', '4')], check=True)


compiled = os.environ.get('RDKIT_ETH_BUILD_RDXTB') == '1'
setup(ext_modules=[Extension('rdkit_eth.rdXTBTools', sources=[])] if compiled else [],
      cmdclass={'build_ext': CMakeBuild} if compiled else {})
