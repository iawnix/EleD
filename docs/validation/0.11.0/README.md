# 0.11 mainline validation

Validated on 2026-11-05. This evidence covers the electronic-descriptor mainline; embedding algorithm patches remain paused.

| Check | Result |
|---|---|
| Source package and extended engine suite | 416 passed, 17 skipped |
| Relocated bundle runtime suite | 416 passed, 17 skipped |
| Focused worker/SCC/checked-QM8 tests | 14 passed, 4 skipped; worker/SCC installed subset 16 passed |
| C++ core test | 1 passed |
| QM8 stereo-refit | 974/1,024 accepted; 9 stored coordinates rejected by audit |
| Checked preparation QM8 | 999/1,024 accepted; 25 no accepted conformer |
| QM9 dipole comparison | 256/256 for GFN1 and GFN2 |
| Engine parity | 33 molecules × 3 methods; no observed numerical delta |

The detailed reports retain exact manifests, predictions, failures, split IDs, source and engine hashes.
