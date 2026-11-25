"""Test-only helpers for the public calculator API."""

from eled import XTBCalculator


def calculate_xtb(mol, **options):
    return XTBCalculator(**options).Calculate(mol)
