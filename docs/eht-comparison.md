# Optional YAeHMOP/EHT comparison

The core `rdkit_eth` protocol uses the project-built extended `libxtb` for
GFN0, GFN1 and GFN2 calculations. RDKit's `rdEHTTools`/YAeHMOP implementation
is a separate pre-existing electronic-structure model.

The optional helpers in `rdkit_eth.descriptors` only normalize YAeHMOP scalar
outputs for side-by-side comparisons:

```python
from rdkit_eth.descriptors import calc_eht_descriptors, calc_orbital_descriptors

# `calc_eht_descriptors` consumes an existing rdEHTTools result.
# `calc_orbital_descriptors(..., backend="eht")` calls RDKit's YAeHMOP API.
```

YAeHMOP and xTB have different Hamiltonians, basis conventions, charge models,
energy definitions and reference zeros. Their total energies and populations
must not be mixed or treated as interchangeable. EHT fields are excluded from
`descriptor_catalog()`, the xTB report schema, the EAS feature protocol and the
extended-libxtb capability checks.

This adapter does not modify YAeHMOP and does not make it a dependency of the
xTB workflow. Users who only need the project xTB descriptors can ignore it.

