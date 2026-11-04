# Official-library compatibility recipe

This historical recipe installs the pure Python adapter with official conda
xTB. It is not the standard 0.5 installation and does not provide the native
basis export or verified GFN0 capability. It has not been validated with
`conda build` or published.

Use `python tools/install_native.py` from the repository root for the tested
full profile: extended libxtb, all supported descriptors and rdXTBTools.
A relocatable conda package for that full profile is future distribution work;
no unpublished dependency is added to this recipe.
