import pytest

from rdkit_eth import XTBCalculator
from rdkit_eth.result import XTBResult


def test_calculator_owns_validated_options_and_repr():
    calculator = XTBCalculator(method='gfn1', electronic_temperature=300.)
    assert calculator.options.method == 'GFN1-xTB'
    assert calculator.to_options() == calculator.options
    assert 'GFN1-xTB' in repr(calculator)


def test_calculator_forwards_rdkit_result_and_per_call_overrides(monkeypatch):
    captured = {}
    expected = XTBResult('GFN2-xTB', None, -1., None, 1, 2, 0)

    def fake_run(mol, **options):
        captured['mol'] = mol
        captured['options'] = options
        return expected

    monkeypatch.setattr('rdkit_eth.calculator.run_xtb', fake_run)
    calculator = XTBCalculator(method='GFN2-xTB', electronic_temperature=300.)
    molecule = object()
    assert calculator.calculate(molecule, conf_id=2, keep_overlap=True) is expected
    assert captured['mol'] is molecule
    assert captured['options']['conf_id'] == 2
    assert captured['options']['keep_overlap'] is True
    assert captured['options']['electronic_temperature'] == 300.


def test_calculator_rejects_unknown_options():
    with pytest.raises(TypeError, match='unknown XTBCalculator option'):
        XTBCalculator(not_an_option=True)
    with pytest.raises(TypeError, match='unknown XTBCalculator override'):
        XTBCalculator().calculate(object(), not_an_option=True)
