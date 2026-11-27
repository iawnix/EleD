import unittest
from types import SimpleNamespace
from unittest.mock import patch

from eled.xtb import _frontier_energies


class AdapterTests(unittest.TestCase):
    def test_unknown_native_api_is_rejected_before_binding_calculations(self):
        from eled.backends.native import _NativeXTB
        from eled.errors import XTBNativeNotFoundError
        # Deliberately expose no calculator functions: an unknown ABI must
        # be rejected on its version before those functions are bound/used.
        library = SimpleNamespace(xtb_getAPIVersion=lambda: 30000)
        with patch('eled.backends.native.ctypes.CDLL', return_value=library):
            with self.assertRaisesRegex(XTBNativeNotFoundError, 'version 30000'):
                _NativeXTB('future-library')

    def test_frontier_energies_use_fractional_occupations(self):
        homo, lumo = _frontier_energies((-2.0, -1.0, 0.5), (2.0, 0.0001, 0.0))
        self.assertEqual(homo, -1.0)
        self.assertEqual(lumo, 0.5)


    def test_implicit_hydrogen_guard_is_kept_at_the_rdkit_boundary(self):
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError:
            self.skipTest("RDKit is optional for the parser-only test suite")
        from helpers import calculate_xtb

        mol = Chem.MolFromSmiles("CC")
        AllChem.EmbedMolecule(mol, randomSeed=1)
        with self.assertRaisesRegex(ValueError, "implicit hydrogens"):
            calculate_xtb(mol, xtbLibrary="/does/not/matter")


if __name__ == "__main__":
    unittest.main()
