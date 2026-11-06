# Project positioning

`rdkit-eth` is an extension layer for RDKit's molecular model. RDKit remains
the source of truth for the molecular graph, atom and bond indices, conformers,
stereochemistry and coordinate units. `rdkit-eth` adds electronic-structure
calculation, analysis and provenance around that model.

The closest software reference is a descriptor library such as DScribe. ASE is
useful as an architectural reference for calculators, but `rdkit-eth` does not
aim to become a general atomistic simulation environment. Periodic systems,
molecular dynamics, trajectory management and broad solvent workflows are
outside the core scope.

## Extension boundary

```text
RDKit Mol / Conformer
        │
        ├── checked geometry preparation
        │
        ├── XTBCalculator → XTBResult
        │                         │
        │                         ├── molecular/atom/bond descriptors
        │                         ├── density and shell descriptors
        │                         ├── fixed-geometry reactivity descriptors
        │                         └── indexed RDKit annotations/reports
        │
        └── optional YAeHMOP comparison helper
```

The public layers should remain independent:

1. **RDKit boundary** validates `ROMol`, conformer IDs, atom order and
   provenance.
2. **Calculator layer** runs a selected electronic-structure backend and owns
   numerical settings and engine identity.
3. **Result layer** stores method-defined arrays, units and same-calculation
   basis data.
4. **Descriptor layer** consumes a result without invoking an engine again.
5. **Report layer** serializes units, missing reasons, quality flags and hashes.

This separation lets RDKit users use basic charge and bond features without
retaining a basis, while density and reactivity analyses request the additional
capabilities explicitly.

## Public API direction

The public calculation API is the reusable calculator object:

```python
from rdkit_eth import XTBCalculator
from rdkit_eth import calc_xtb_descriptors

calculator = XTBCalculator(method="GFN2-xTB", electronic_temperature=300.0)
result = calculator.calculate(mol, keep_overlap=True)
features = calc_xtb_descriptors(result)
```

Descriptor objects or registry entries should expose a common contract for
name, units, required result capabilities, invariance and missing-value
semantics. The current `descriptor_catalog()` and `descriptor_status()` are the
initial registry implementation.

## Packaging direction

The Python/RDKit descriptor API and the native xTB runtime should be separable:

- the core package contains result models, descriptor extraction, reports and
  RDKit mapping;
- an optional native extra or conda package supplies the reviewed extended
  `libxtb` and compiled `rdXTBTools` module;
- source installation remains available for developers and reproducible Linux
  research environments;
- the standard wheel should not contain a complete compiler and conda runtime.

This preserves a small importable RDKit extension while making the native
calculation profile explicit and testable.

## Scope of the next releases

The next API work should be spent on a calculator facade, descriptor registry,
backend capability discovery, platform packaging and examples. Geometry
optimization remains an optional coordinate-preparation tool. New descriptors
must answer a defined RDKit analysis need and pass unit, invariance, conservation
and ctypes/compiled parity checks before entering the recommended catalog.

The EAS benchmark is a validation application of this extension layer. It is
kept in research tooling and does not define the package's core API.
