# 0.5 native dependency migration and installation evidence

Validated locally on 2026-11-04 (Asia/Shanghai). Source revision `33fe381`
contains upstream's GFN0 repair; only our basis-export patch is applied.

- `migration-*.json`: new versus verified 0.4 engine, 33 molecules × three
  methods, separate processes. All recorded numerical differences are zero.
- `python-*.log`: full adapter tests with new, legacy and official libraries.
- `python-installed.log`: full tests of the installed wheel, isolated Python
  outside the checkout with source-path configuration disabled.
- `installed-wheel-gfn0.json`: fresh-venv wheel check of both interfaces.
- `installation.json`: six installed checks without xTB path overrides.
- `xtb-build.json` / `xtb-ctest.log`: exact engine, compiler and test identity.
- `cpp-core.log`: standalone C++ verification of the standard engine.

Paths identify local prefixes and are not portable installation instructions.
See `docs/install.md` to reproduce the complete profile. These records establish
API/numerical contracts, not excited-state prediction accuracy. They do not
publish binaries or certify arbitrary xTB versions.
