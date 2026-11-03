# Local 0.4.0 validation snapshot

These records preserve the pre-repository 0.4.0 build and installed-wheel
checks: 205 adapter tests passed with the extended engine, 122 with official
xTB, plus upstream and standalone C++ checks. They are historical local
results, not GitHub Actions results.

Paths identify the environments used at validation time. Temporary wheel-test
venvs have been removed after verification. The artifact hashes refer to the
original local `dist/0.4.0` source and wheels; those binaries and environments
are excluded from Git. Reproduce builds using `docs/install.md` and the locked
environments. Source metadata/documentation have since been updated to record
the upstream fix; these are not hashes of the initial Git commit.
