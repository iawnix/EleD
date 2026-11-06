# Literature and positioning audit

This audit fixes the wording used for the new descriptors. It is a scope
check, not a claim that the descriptors are experimentally unique observables.

| Feature | Established baseline | xTB + ML status | Position in `rdkit-eth` |
| --- | --- | --- | --- |
| Shell-resolved s/p/d populations | Mulliken and Löwdin population analyses are standard; xTB supplies a minimal-basis wavefunction | Charge and orbital features are used in molecular ML, but a portable xTB C API for shell maps/populations is not a standard public feature | **new extraction layer**, validated by `trace(P S)` and shell-to-atom sums |
| π population / π bond order | π projections are standard in localized orbital and bond analysis | No general GFN1/GFN2 xTB API convention for local p-normal projection or RDKit bond mapping was found | **method-defined feature**, only reported when the plane and p shell are defined |
| Frontier orbital bond decomposition | Related signed overlap and bond analyses exist in MO/ETS/COHP workflows | No established xTB + ML benchmark or stable xTB C API field for per-RDKit-bond HOMO/LUMO decomposition | **exploratory feature**, group-projector invariant and explicitly not an excitation descriptor |
| Atomic conformer dispersion | Ensemble and conformer sensitivity statistics are common in cheminformatics | Usually reported as a workflow diagnostic rather than an xTB-specific learned feature | **reliability diagnostic**, with equal-weight and optional electronic-energy weighting |
| Fukui / local conceptual-DFT indices | Condensed Fukui indices and global softness/electrophilicity are established | GFN1/GFN2 ΔSCF and charge/shell baselines are used for reactivity ML; GFN0 EEQ charges are not equivalent to density populations | **analysis of existing ΔSCF results**, with units and state-failure semantics |

The primary validation target is atom-level EAS regioselectivity using the
public RegioML/RegioSQM records. This target can test whether shell and
frontier features add information beyond a GFN1 xTB partial-charge/shell baseline. QM8
E1 is retained as a molecular smoke test only: it cannot establish utility for
an atom-labelled target. CYP site-of-metabolism is the fallback target if the
EAS atom mapping cannot be made reproducible.

The implementation follows the definitions in Mulliken, *J. Chem. Phys.* 23,
1833 (1955); Löwdin, *J. Chem. Phys.* 18, 365 (1950); Mayer, *Chem. Phys.
Lett.* 97, 270 (1983); Yang and Mortier, *J. Am. Chem. Soc.* 108, 5708
(1986); Parr and Yang, *Density-Functional Theory of Atoms and Molecules*
(1989); and the GFN-xTB method papers by Grimme and co-workers (GFN1-xTB,
*J. Chem. Theory Comput.* 13, 1989 (2017); GFN2-xTB, *J. Chem. Theory
Comput.* 15, 1652 (2019)).

These references support the mathematical definitions and the xTB electronic
structure model. They do not establish that every new per-atom descriptor is
predictive. Predictivity is therefore an empirical question fixed by
[`protocol.json`](../protocol.json), with failures and negative results
reported alongside improvements.
