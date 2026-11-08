# Conda recipe status

This recipe is retained as packaging scaffolding and has not been validated
with `conda build` or published. The `xtb` package from conda-forge is an
unmodified upstream library and therefore cannot satisfy the runtime contract
of `eled` by itself: calculations require the project-built extended
`libxtb` with basis export, SCC diagnostics and the verified GFN0 state fix.

Use `python tools/install_native.py` from the repository root for the tested
full profile: extended libxtb, all supported descriptors and rdXTBTools.
A relocatable conda package for that full profile is future distribution work.
