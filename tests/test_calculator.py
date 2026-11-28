import pytest
import eled

from eled import XTBCalculator, XTBResults


def test_calculator_owns_validated_options_and_repr():
    calculator = XTBCalculator(method='gfn1', electronicTemperature=300.)
    assert calculator.GetOptions().method == 'GFN1-xTB'
    assert not hasattr(calculator, 'options') and not callable(calculator)
    assert 'GFN1-xTB' in repr(calculator)


def test_first_public_api_has_no_legacy_run_xtb_entrypoint():
    assert not hasattr(eled, 'run_xtb')


def test_calculator_forwards_rdkit_result_and_per_call_overrides(monkeypatch):
    captured = {}
    expected = XTBResults('GFN2-xTB', None, -1., None, 1, 2, 0)

    def fake_run(mol, **options):
        captured['mol'] = mol
        captured['options'] = options
        return expected

    monkeypatch.setattr('eled.calculator._run_xtb', fake_run)
    calculator = XTBCalculator(method='GFN2-xTB', electronicTemperature=300.)
    molecule = object()
    assert calculator.Calculate(molecule, confId=2, keepOverlap=True) is expected
    assert captured['mol'] is molecule
    assert captured['options']['confId'] == 2
    assert captured['options']['keepOverlap'] is True
    assert captured['options']['electronicTemperature'] == 300.
    assert calculator(molecule, confId=5) is expected
    assert captured['options']['confId'] == 5


def test_calculator_rejects_unknown_options():
    with pytest.raises(TypeError, match='unknown XTBCalculator option'):
        XTBCalculator(not_an_option=True)
    with pytest.raises(TypeError, match='unknown XTBCalculator override'):
        XTBCalculator().Calculate(object(), not_an_option=True)
