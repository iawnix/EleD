"""Test-only helpers for the first public calculator API."""

from eled import XTBCalculator


def calculate_xtb(mol, **options):
    return XTBCalculator(**options).calculate(mol)
