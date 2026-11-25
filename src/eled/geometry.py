"""Auditable SMILES-to-3D preparation with checks of specified stereochemistry."""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
import json
from typing import Any

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from .input import _geometry_hash, _ordered_structure_hash, _xyz_from_mol

_METHODS = {n.lower(): n for n in ('ETKDG', 'ETKDGv2', 'ETKDGv3', 'ETKDGv4', 'srETKDGv3', 'KDG')}
_DEFAULT_FALLBACKS = ('ETKDGv4', 'srETKDGv3')


@dataclass(frozen=True)
class ConformerGenerationResult:
    """Owned explicit-H molecule and JSON report; no conformers on failure."""

    molecule: Any
    report: dict

    @property
    def success(self) -> bool:
        return self.report['n_conformers'] > 0

    def GetMolecule(self) -> Any:
        return self.molecule

    def GetReport(self) -> dict:
        return self.report


def _integer(value, name, minimum=1, maximum=2147483647):
    if isinstance(value, bool) or not isinstance(value, Integral) or not minimum <= value <= maximum:
        raise ValueError(f'{name} must be an integer in [{minimum}, {maximum}]')
    return int(value)


def _method(name):
    if not isinstance(name, str) or name.strip().lower() not in _METHODS:
        raise ValueError(f'unknown embedding preset: {name!r}')
    return _METHODS[name.strip().lower()]


def _failure_counts(params):
    """Keep raw slots: new RDKit releases may add counters without enum names."""
    if not hasattr(params, 'GetFailureCounts'):
        return {'available': False, 'raw': [], 'named': {}}
    counts = list(map(int, params.GetFailureCounts()))
    names = {int(v): str(v) for v in AllChem.EmbedFailureCauses.values.values()}
    return {
        'available': True,
        'raw': counts,
        'named': {names.get(i, f'UNKNOWN_SLOT_{i}'): n for i, n in enumerate(counts) if n},
    }


def _stereo_constraints(mol):
    tetra = (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
    if mol.GetStereoGroups():
        raise ValueError('enhanced stereo groups are not supported by geometry preparation')
    if any(a.GetChiralTag() not in (*tetra, Chem.ChiralType.CHI_UNSPECIFIED) for a in mol.GetAtoms()):
        raise ValueError('only tetrahedral atom stereochemistry is supported')
    bonds = {}
    for b in mol.GetBonds():
        s = b.GetStereo()
        if s in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ):
            bonds[b.GetIdx()] = str(s)
        elif s != Chem.BondStereo.STEREONONE:
            raise ValueError('only specified E/Z double-bond stereochemistry is supported')
    return {
        'atoms': {a.GetIdx(): str(a.GetChiralTag()) for a in mol.GetAtoms() if a.GetChiralTag() in tetra},
        'bonds': bonds,
    }


def _check_geometry(mol, confId, constraints):
    coords = np.asarray(mol.GetConformer(confId).GetPositions())
    if not mol.GetConformer(confId).Is3D() or not np.isfinite(coords).all():
        return {'ok': False, 'reason': 'invalid_coordinates'}
    if len(coords) > 1:
        distances = np.linalg.norm(coords[:, None] - coords[None, :], axis=2)
        distances[np.diag_indices(len(coords))] = np.inf
        if distances.min() < 0.3:
            return {'ok': False, 'reason': 'atomic_collision', 'minimum_distance_angstrom': float(distances.min())}
    observed = Chem.Mol(mol)
    Chem.RemoveStereochemistry(observed)
    Chem.AssignStereochemistryFrom3D(observed, confId=confId, replaceExistingTags=True)
    atoms = [i for i, tag in constraints['atoms'].items() if str(observed.GetAtomWithIdx(i).GetChiralTag()) != tag]
    bonds = [i for i, stereo in constraints['bonds'].items() if str(observed.GetBondWithIdx(i).GetStereo()) != stereo]
    return {
        'ok': not atoms and not bonds,
        'reason': 'stereo_mismatch' if atoms or bonds else None,
        'mismatched_atom_indices': atoms,
        'mismatched_bond_indices': bonds,
    }


def PrepareConformers(
    source,
    *,
    numConformers=1,
    method='ETKDGv3',
    fallbackMethods=None,
    randomSeed=42,
    attemptsPerStrategy=2,
    retryRandomCoords=True,
    maxEmbeddingIterations=300,
    optimize='MMFF94',
    optimizationIterations=500,
    requireAssignedStereo=False,
    xtbMethod='GFN2-xTB',
    xtbLibrary=None,
    charge=None,
    numUnpairedElectrons=None,
):
    """Generate a checked conformer pool from SMILES or a copied RDKit molecule."""
    requested = _integer(numConformers, 'numConformers')
    attempts = _integer(attemptsPerStrategy, 'attemptsPerStrategy')
    seed = _integer(randomSeed, 'randomSeed', minimum=0, maximum=2147483646)
    embed_iterations = _integer(maxEmbeddingIterations, 'maxEmbeddingIterations')
    iterations = _integer(optimizationIterations, 'optimizationIterations')
    if not isinstance(retryRandomCoords, bool) or not isinstance(requireAssignedStereo, bool):
        raise ValueError('retryRandomCoords and requireAssignedStereo must be Boolean')
    if optimize not in (None, 'MMFF94', 'UFF', 'xTB'):
        raise ValueError("optimize must be None, 'MMFF94', 'UFF' or 'xTB'")
    primary = _method(method)
    if not hasattr(AllChem, primary):
        raise ValueError(f'installed RDKit lacks {primary}')
    defaults = fallbackMethods is None
    if isinstance(fallbackMethods, (str, bytes)):
        raise ValueError('fallbackMethods must be a sequence of preset names, not a string')
    methods, unavailable = [primary], []
    for value in _DEFAULT_FALLBACKS if defaults else fallbackMethods:
        name = _method(value)
        if not hasattr(AllChem, name):
            if not defaults:
                raise ValueError(f'installed RDKit lacks {name}')
            unavailable.append(name)
        elif name not in methods:
            methods.append(name)
    if isinstance(source, str):
        mol = Chem.MolFromSmiles(source)
        if mol is None:
            raise ValueError('SMILES could not be parsed by RDKit')
        source_kind, source_smiles = 'smiles', source
    elif isinstance(source, Chem.Mol):
        mol = Chem.Mol(source)
        Chem.SanitizeMol(mol)
        source_kind, source_smiles = 'molecule', Chem.MolToSmiles(mol, isomericSmiles=True)
    else:
        raise TypeError('source must be a SMILES string or an RDKit Mol')
    if not mol.GetNumAtoms() or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ValueError('input must contain real atoms and no dummy atoms')
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError(
            'disconnected fragments need an explicit relative geometry; supply a prepared molecule to XTBCalculator'
        )
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    unassigned = [
        {'type': str(s.type), 'centered_on': int(s.centeredOn)}
        for s in Chem.FindPotentialStereo(mol)
        if str(s.specified) != 'Specified'
    ]
    if unassigned and requireAssignedStereo:
        raise ValueError('input has unspecified stereochemistry')
    discarded = mol.GetNumConformers()
    mol.RemoveAllConformers()
    mol = Chem.AddHs(mol)
    constraints = _stereo_constraints(mol)
    report = {
        'schema_version': 1,
        'protocol_version': '1.0',
        'report_type': 'conformer_generation',
        'status': 'failed',
        'n_requested': requested,
        'n_conformers': 0,
        'n_missing': requested,
        'input_smiles': source_smiles,
        'canonical_smiles': Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True),
        'source_kind': source_kind,
        'discarded_input_conformers': discarded,
        'source_coordinates_used': False,
        'rdkit_version': rdBase.rdkitVersion,
        'specified_stereo': {
            k: {str(i): v for i, v in values.items()} for k, values in constraints.items()
        },
        'unspecified_stereo': unassigned,
        'warnings': ['unspecified_stereochemistry'] if unassigned else [],
        'retry_policy': {
            'methods': methods,
            'unavailable_default_fallbacks': unavailable,
            'random_seed': seed,
            'attempts_per_strategy': attempts,
            'retry_random_coords': retryRandomCoords,
            'max_embedding_iterations': embed_iterations,
            'prune_rms_angstrom': 0.2,
            'num_threads': 1,
            'enforce_chirality': True,
            'selection': 'first strategy with any accepted conformers',
            'require_assigned_stereo': requireAssignedStereo,
        },
        'optimization': optimize,
        'optimization_iterations': iterations if optimize else None,
        'charge': charge,
        'uhf': numUnpairedElectrons,
        'attempts': [],
        'conformers': [],
        'recommended_conformer_id': None,
        'error': None,
    }
    output = Chem.Mol(mol)
    if optimize in ('MMFF94', 'UFF'):
        has_parameters = (
            AllChem.MMFFHasAllMoleculeParams(mol)
            if optimize == 'MMFF94'
            else AllChem.UFFHasAllMoleculeParams(mol)
        )
        if not has_parameters:
            report['error'] = {'reason': 'force_field_parameters_unavailable', 'force_field': optimize}
            return ConformerGenerationResult(output, report)
    for preset in methods:
        for random_coords in ([False, True] if retryRandomCoords else [False]):
            for retry in range(attempts):
                remaining = requested - output.GetNumConformers()
                if remaining <= 0:
                    break
                current_seed = (seed + 104729 * len(report['attempts'])) % 2147483647
                params = getattr(AllChem, preset)()
                params.randomSeed = current_seed
                params.numThreads = 1
                params.pruneRmsThresh = 0.2
                params.enforceChirality = True
                params.useRandomCoords = random_coords
                params.maxIterations = embed_iterations
                params.trackFailures = True
                working = Chem.Mol(mol)
                ids = list(AllChem.EmbedMultipleConfs(working, numConfs=remaining, params=params))[:remaining]
                attempt = {
                    'method': preset,
                    'retry': retry,
                    'random_seed': current_seed,
                    'use_random_coords': random_coords,
                    'enforce_chirality': True,
                    'n_embedded': len(ids),
                    'n_accepted': 0,
                    'failure_counts': _failure_counts(params),
                    'candidates': [],
                    'status': 'embedding_failed' if not ids else 'no_valid_conformers',
                }
                report['attempts'].append(attempt)
                for conf_id in ids:
                    check = _check_geometry(working, conf_id, constraints)
                    record = {
                        'embedded_conformer_id': int(conf_id),
                        'before_optimization': check,
                        'optimization_status': None,
                        'energy_kcal_mol': None,
                        'accepted': False,
                    }
                    attempt['candidates'].append(record)
                    if not check['ok']:
                        continue
                    candidate = Chem.Mol(mol)
                    candidate.AddConformer(Chem.Conformer(working.GetConformer(conf_id)), assignId=True)
                    if optimize in ('MMFF94', 'UFF'):
                        outcomes = (
                            AllChem.MMFFOptimizeMoleculeConfs(
                                candidate, numThreads=1, maxIters=iterations, mmffVariant='MMFF94'
                            )
                            if optimize == 'MMFF94'
                            else AllChem.UFFOptimizeMoleculeConfs(candidate, numThreads=1, maxIters=iterations)
                        )
                        status, energy = outcomes[0]
                        if not np.isfinite(energy):
                            raise ValueError('force field returned nonfinite energy')
                        record.update(optimization_status=int(status), energy_kcal_mol=float(energy))
                        if status != 0:
                            record['rejection'] = 'force_field_not_converged'
                            continue
                    elif optimize == 'xTB':
                        from .optimization import OptimizeXTB
                        optimized = OptimizeXTB(
                            candidate,
                            confId=0,
                            method=xtbMethod,
                            charge=charge,
                            numUnpairedElectrons=numUnpairedElectrons,
                            xtbLibrary=xtbLibrary,
                            maxIterations=iterations,
                        )
                        status = 0 if optimized.success else 1
                        energy = optimized.report['final_energy_hartree']
                        record.update(
                            optimization_status=status,
                            energy_hartree=float(energy),
                            energy_kcal_mol=float(energy * 627.5094740631),
                            xtb_optimization=optimized.report,
                        )
                        if not optimized.success:
                            record['rejection'] = 'xtb_not_converged'
                            continue
                        candidate = optimized.molecule
                    check = _check_geometry(candidate, 0, constraints)
                    record['after_optimization'] = check
                    if not check['ok']:
                        continue
                    new_id = output.AddConformer(Chem.Conformer(candidate.GetConformer()), assignId=True)
                    record.update(accepted=True, conformer_id=new_id)
                    report['conformers'].append({
                        'conformer_id': new_id,
                        'attempt_index': len(report['attempts']) - 1,
                        'embedded_conformer_id': int(conf_id),
                        'energy_kcal_mol': record['energy_kcal_mol'],
                        'geometry_sha256': _geometry_hash(_xyz_from_mol(output, new_id)),
                    })
                    attempt['n_accepted'] += 1
                if attempt['n_accepted']:
                    attempt['status'] = 'accepted'
                    report.update(
                        n_conformers=output.GetNumConformers(),
                        n_missing=requested - output.GetNumConformers(),
                        status='complete' if output.GetNumConformers() == requested else 'partial',
                        ordered_structure_sha256=_ordered_structure_hash(output),
                    )
                    report['recommended_conformer_id'] = (
                        min(report['conformers'], key=lambda r: (r['energy_kcal_mol'], r['conformer_id']))[
                            'conformer_id'
                        ]
                        if optimize
                        else 0
                    )
                    json.dumps(report, allow_nan=False)
                    return ConformerGenerationResult(output, report)
    report['error'] = {'reason': 'no_accepted_conformer'}
    json.dumps(report, allow_nan=False)
    return ConformerGenerationResult(output, report)
