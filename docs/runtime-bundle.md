# Offline Linux runtime bundle (0.11)

The bundle contains the matched RDKit/Python runtime, compiled `rdXTBTools`,
extended libxtb, parameter files and its separate native dependency prefix.
Installation needs neither conda, a compiler, network access nor an xTB CLI.
It is a Linux x86-64 environment distribution, not a portable manylinux wheel.
Target CPU/glibc compatibility still needs testing on each supported platform.
Windows/macOS and package-index publication remain outside this release.

Given a trusted bundle directory containing `install.py` and `bundle.json`:

```bash
python3 /path/to/bundle/install.py --directory /absolute/new/eled
/absolute/new/eled/runtime/bin/python -m eled doctor --check compiled --keep-overlap
```

Use Python 3.10+ with tarfile's backported `data` extraction filter (or Python
3.12+) to run the installer. Choose a new destination. The installer verifies
payload checksums, extracts both prefixes, runs conda's relocation code and
repairs cross-prefix ELF RPATHs. It preserves the original binary/build-record
hashes in relocation provenance and fingerprints the relocated library anew.
It then checks GFN0/1/2, overlap/populations and GFN1/2 reactivity through both
interfaces from outside the checkout, with library-selection overrides removed.
The bundle builder also records hashes for the xTB executable and all GFN
parameter files, so the installed runtime contains the complete extended xTB
profile rather than only `libxtb.so`.
It also rejects unresolved or original-prefix dependencies in `ldd` output.

The resulting `installation.json` records the checks. Keep `runtime` and
`toolchain` together at their installed paths. To install somewhere else, rerun
the installer from the original bundle into a new directory. Failed installs
are left available for inspection and never overwritten automatically.

Python compiled calls now use the same active-runtime libxtb discovery as the
ctypes interface. An explicit `xtbLibrary` or `XTB_LIBRARY` remains authoritative.
Standalone C++ users should provide `XTBOptions.libraryPath` when relocating
their own linked application.

## Build the bundle

First create and verify a complete source installation. In a separate packaging
environment, install the build-only tools:

```bash
python -m pip install --require-hashes -r tools/distribution-requirements.txt
python tools/build_runtime_bundle.py \
  --runtime /path/to/installation/runtime \
  --toolchain /path/to/installation/toolchain \
  --output /new/bundle-directory
```

The builder retains conda metadata/licenses, supplies the exact patched xTB
source and its licenses, and includes the adapter's BSD-3-Clause license.
xTB-derived patches retain the upstream licensing terms. The bundle is large
because it preserves the locked SDK and dependency environment; it is not a
minimal runtime image. No artifact is uploaded by either tool.

Checksums detect accidental changes and associate the installation with the
selected artifact; they are not a publisher signature. Use artifacts from a
trusted source.
