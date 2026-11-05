# Upstream GFN0 state fix — checked 2026-11-04 (Asia/Shanghai)

> 0.5 migration: the current source lock now pins the merged repair and applies
> the native basis-export and SCC-diagnostics patches. The historical status check below describes
> 0.4; see [current validation](../validation.md) for the migration results.

Our [issue #1486](https://github.com/grimme-lab/xtb/issues/1486) is closed as
completed. Maintainer Thomas Froitzheim's
[PR #1490](https://github.com/grimme-lab/xtb/pull/1490) was merged on
2026-11-02 at 22:26:51 Asia/Shanghai (14:26:51 UTC), at commit
`33fe381adde775c4992610c0f90d8914b7572bb5`. The maintainer explicitly confirmed
that obsolete parameter selection overwrote D4 reference charges.

The change removes the first `select case(mode)` block that mutated module
reference arrays, while retaining the model-owned `dispm%q` / local `tmp_hq`
selection. It adds a C API regression checking GFN2 energy and every orbital
energy before and after loading GFN0 with fresh handles, and a Fortran GFN2
test that initializes GFN0 first.

## Relationship to this project's patches

Both repairs remove global `refq`/`refh` mutation. They are not numerically
identical implementations: our patch retains mode-specific secondary charges
in local `secondary_q`, whereas upstream replaces `secq(is)+iz` by `iz` in
`zeta`. Therefore a migration needs fresh GFN0/1/2 numerical comparisons as
well as the method-switch regression; the absence of state contamination
alone does not prove all energies are unchanged.

The current 0.4 source lock still pins `a9040512`, which predates the fix, and
requires the reviewed local D4 patch. Its binary capability hashes remain
valid for that tested build. Updating the lock, builder and Python/C++
capability checks is a separate dependency migration. Do not apply the old
D4 patch blindly on top of the repaired revision or remove GFN0 guards on
version-string evidence alone.

The native overlap/AO/shell mapping extension is independent of this bug fix.
PR #1490 does not add our basis-export C API; population descriptors still
need that extension, even after migrating to an upstream-fixed engine.

At this check GitHub's latest published release remains **v6.7.1**, published
2024-07-23. A fix merged into main does not update already installed 6.7.1
binaries, and this check does not establish downstream package backports.
The merged revision was reviewed at source level, not independently built
or numerically revalidated here.

[Machine-readable status](xtb-gfn0-upstream-fix.json) and the
[reviewed PR diff](xtb-pr-1490.diff) are retained. The submitted issue body is
preserved unchanged in [the original report](xtb-gfn0-state-issue.md).
