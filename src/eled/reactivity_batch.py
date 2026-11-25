"""Dataset-safe vertical reactivity: explicit jobs, successful cache, resumable rows."""
from __future__ import annotations

from contextlib import closing, nullcontext
from copy import deepcopy
from dataclasses import dataclass, replace
import hashlib
from numbers import Integral

from .charges import CHARGE_MODEL_CM5, CHARGE_MODEL_SCC
from .electronic_descriptors import _fragments
from .errors import DescriptorUnavailableError, XTBCacheError
from .input import _effective_charge_uhf, _input_metadata, _xyz_from_mol
from .options import XTBOptions
from .persistence import ReactivityCheckpoint, XTBDiskCache, _digest, _implementation_identity
from .reactivity import (
    REACTIVITY_PROTOCOL_VERSION,
    REACTIVITY_UNITS,
    RunReactivity,
    _empty_report,
    _population_model,
    _scheme,
)
from .reports import DESCRIPTOR_VERSION, SCHEMA_VERSION
from .runtime import _engine_identity, _parameter_context, _parameter_file
from .selections import NormalizeSelections
from .xtb import _resolve_options


@dataclass(frozen=True)
class ReactivityJob:
    """One indexed molecule and its per-molecule settings. IDs must be unique."""

    recordId: str
    mol: object
    options: XTBOptions | None = None
    ionizedUnpaired: int | None = None
    attachedUnpaired: int | None = None
    fragments: dict | None = None
    reactionSites: object = None


def _prepare_jobs(jobs, options, populationScheme, interface, initialization):
    from rdkit import Chem
    populationScheme = _scheme(populationScheme)
    if interface not in ('native', 'compiled'):
        raise ValueError('interface must be native or compiled')
    if initialization not in ('reference', 'independent'):
        raise ValueError('initialization must be reference or independent')
    if options is not None and not isinstance(options, XTBOptions):
        raise TypeError('options must be XTBOptions')
    jobs = list(jobs)
    ids = []
    for job in jobs:
        if not isinstance(job, ReactivityJob):
            raise TypeError('jobs must contain ReactivityJob objects')
        if not isinstance(job.recordId, str) or not job.recordId.strip():
            raise ValueError('recordId must be a nonempty string')
        ids.append(job.recordId)
    if len(set(ids)) != len(ids):
        raise ValueError('record IDs must be unique within a dataset')
    implementation = _implementation_identity(interface)
    prepared = []
    for job in jobs:
        config = job.options if job.options is not None else options or XTBOptions()
        if not isinstance(config, XTBOptions):
            raise TypeError('job options must be XTBOptions')
        if config.method not in ('GFN1-xTB', 'GFN2-xTB'):
            raise ValueError('reactivity requires GFN1/2')
        for spin in (job.ionizedUnpaired, job.attachedUnpaired):
            if spin is not None and (isinstance(spin, bool) or not isinstance(spin, Integral) or not 0 <= spin <= 2147483647):
                raise ValueError('adjacent spins must be nonnegative integers')
        config = _resolve_options(
            replace(config, keepOverlap=config.keepOverlap or populationScheme in ('mulliken', 'lowdin'))
        )
        engine = _engine_identity('native', config.xtbLibrary)
        if not engine['sha256']:
            raise ValueError('persistent protocols require a fingerprintable libxtb file')
        if config.keepOverlap and not engine['orbital_basis']:
            raise DescriptorUnavailableError('population analysis requires the verified native basis extension')
        parameters = _parameter_file(config.method, config.xtbLibrary)
        if job.mol is not None and not isinstance(job.mol, Chem.Mol):
            raise TypeError('job molecule must be an RDKit Mol or None')
        mol = None if job.mol is None else Chem.Mol(job.mol)
        if mol is not None and mol.GetNumAtoms():
            charge, spin = _effective_charge_uhf(mol, config.charge, config.numUnpairedElectrons)
            if spin and (job.ionizedUnpaired is None or job.attachedUnpaired is None):
                raise ValueError(f'{job.recordId}: open-shell references require explicit adjacent spins')
            _effective_charge_uhf(mol, charge + 1, 1 if job.ionizedUnpaired is None else job.ionizedUnpaired)
            _effective_charge_uhf(mol, charge - 1, 1 if job.attachedUnpaired is None else job.attachedUnpaired)
            groups = None if job.fragments is None else {
                k: v.tolist() for k, v in _fragments(job.fragments, mol.GetNumAtoms()).items()
            }
            if job.reactionSites is None:
                sites = None
            elif hasattr(job.reactionSites, 'reactionSites'):
                sites = {
                    k: list(v) for k, v in NormalizeSelections(mol, job.reactionSites).reactionSites.items()
                }
            else:
                sites = {
                    k: list(v)
                    for k, v in NormalizeSelections(mol, {'reactionSites': job.reactionSites}).reactionSites.items()
                }
        else:
            groups = deepcopy(job.fragments)
            sites = deepcopy(job.reactionSites)
        input_error = None
        try:
            if mol is None:
                raise ValueError('molecule is None (upstream parsing failed)')
            xyz = _xyz_from_mol(mol, config.confId)
            metadata = _input_metadata(mol, config.confId, charge, spin, xyz)
            config = replace(config, charge=charge, numUnpairedElectrons=spin, confId=metadata['conformer_id'])
        except ValueError as exc:
            input_error = {'type': type(exc).__name__, 'message': str(exc)}
            metadata = {
                'invalid_input_sha256': None if mol is None else hashlib.sha256(mol.ToBinary()).hexdigest()
            }
        protocol = {
            'cache_schema': 2,
            'protocol_version': REACTIVITY_PROTOCOL_VERSION,
            'schema_version': SCHEMA_VERSION,
            'descriptor_version': DESCRIPTOR_VERSION,
            'implementation': implementation,
            'input': metadata,
            'options': config.ToDict(),
            'engine': engine,
            'parameter_file': parameters,
            'parameter_context': _parameter_context(),
            'interface': interface,
            'initialization': initialization,
            'population_scheme': populationScheme,
            'charge_model': CHARGE_MODEL_CM5 if populationScheme == 'cm5' else CHARGE_MODEL_SCC,
            'population_model': _population_model(populationScheme),
            'ionized_unpaired': None if job.ionizedUnpaired is None else int(job.ionizedUnpaired),
            'attached_unpaired': None if job.attachedUnpaired is None else int(job.attachedUnpaired),
            'fragments': groups,
            'reaction_sites': sites,
        }
        prepared.append((job.recordId, mol, config, protocol, _digest(protocol), input_error))
    return prepared


def RunReactivityBatch(
    jobs,
    *,
    options=None,
    populationScheme='xtb_scc',
    interface='native',
    initialization='reference',
    cache=None,
    checkpoint=None,
    retryFailed=False,
    maxWorkers=1,
    timeout=None,
    isolate=False,
):
    """Run ordered ReactivityJob records with molecule-level checkpoint commits."""
    from .errors import XTBTimeoutError, XTBWorkerError
    from .workers import isolated_map, molecule_bytes, reactivity, validate_execution
    isolated = validate_execution(maxWorkers, timeout, isolate)
    if not isinstance(retryFailed, bool):
        raise ValueError('retryFailed must be Boolean')
    if cache is not None and not isinstance(cache, XTBDiskCache):
        raise TypeError('cache must be XTBDiskCache')
    populationScheme = _scheme(populationScheme)
    prepared = _prepare_jobs(jobs, options, populationScheme, interface, initialization)
    if checkpoint is not None and not isinstance(checkpoint, ReactivityCheckpoint):
        checkpoint = ReactivityCheckpoint(checkpoint)
    if cache is not None and checkpoint is not None and cache.path == checkpoint.path:
        raise ValueError('cache and checkpoint must use separate files')
    manifest = {
        'schema_version': 2,
        'rows': [{'record_id': p[0], 'key': p[4], 'protocol': p[3]} for p in prepared],
    }
    reports = [None] * len(prepared)
    tasks, pending = [], []

    def finish(index, report, source, attempt):
        identifier, mol, config, protocol, key, input_error = prepared[index]
        report['record_id'] = identifier
        report['provenance']['batch_protocol'] = protocol
        if source == 'calculation' and cache is not None and report['error'] is None:
            cache.Set(key, report)
        report['execution'] = {'source': source, 'key': key, 'attempt': attempt}
        if checkpoint is not None:
            checkpoint._save(identifier, key, report, attempt)
        reports[index] = report

    with checkpoint._locked() if checkpoint is not None else nullcontext():
        if checkpoint is not None:
            checkpoint._bind(manifest)
        for index, (identifier, mol, config, protocol, key, input_error) in enumerate(prepared):
            old = None if checkpoint is None else checkpoint._store.read('rows', identifier)
            if old is not None:
                if old['key'] != key:
                    raise XTBCacheError('checkpoint row does not match the manifest')
                if old['report']['error'] is None or not retryFailed:
                    report = deepcopy(old['report'])
                    report['execution']['source'] = 'checkpoint'
                    reports[index] = report
                    continue
            attempt = 1 if old is None else old['attempt'] + 1
            report = None if cache is None or input_error is not None else cache.Get(key)
            source = 'cache' if report is not None else 'calculation'
            if report is not None and (
                not isinstance(report, dict)
                or report.get('provenance', {}).get('batch_protocol') != protocol
            ):
                raise XTBCacheError('cached report protocol does not match the requested calculation')
            if input_error is not None:
                source = 'input_validation'
                report = _empty_report(identifier, populationScheme, {})
                report['error'] = dict(input_error, stage='input')
                report['quality'] = {
                    'status': 'error',
                    'flags': [{
                        'code': 'invalid_molecular_input',
                        'severity': 'error',
                        'detail': input_error['message'],
                    }],
                }
                report['missing'] = dict.fromkeys(REACTIVITY_UNITS, 'invalid_molecular_input')
                report['analysis_settings'].update(
                    protocol_version=REACTIVITY_PROTOCOL_VERSION,
                    initialization=initialization,
                    requested_options=config.ToDict(),
                    interface=interface,
                )
            elif report is None:
                kwargs = dict(
                    options=config,
                    ionizedUnpaired=protocol['ionized_unpaired'],
                    attachedUnpaired=protocol['attached_unpaired'],
                    fragments=protocol['fragments'],
                    reactionSites=protocol['reaction_sites'],
                    populationScheme=populationScheme,
                    interface=interface,
                    initialization=initialization,
                    returnErrors=True,
                )
                if isolated:
                    tasks.append((molecule_bytes(mol), kwargs))
                    pending.append((index, attempt))
                    continue
                report = RunReactivity(mol, **kwargs)
            finish(index, report, source, attempt)
        if isolated:
            with closing(isolated_map(reactivity, tasks, maxWorkers=maxWorkers, timeout=timeout)) as completed:
                for task_index, report in completed:
                    index, attempt = pending[task_index]
                    identifier, mol, config, protocol, key, input_error = prepared[index]
                    if isinstance(report, (XTBTimeoutError, XTBWorkerError)):
                        error = report
                        report = _empty_report(identifier, populationScheme, {})
                        report['error'] = {'type': type(error).__name__, 'message': str(error), 'stage': 'worker'}
                        report['quality'] = {
                            'status': 'error',
                            'flags': [{
                                'code': 'isolated_execution_failed',
                                'severity': 'error',
                                'detail': str(error),
                            }],
                        }
                        report['missing'] = dict.fromkeys(REACTIVITY_UNITS, 'isolated_execution_failed')
                        report['analysis_settings'].update(
                            protocol_version=REACTIVITY_PROTOCOL_VERSION,
                            initialization=initialization,
                            requested_options=config.ToDict(),
                            interface=interface,
                        )
                    elif isinstance(report, Exception):
                        raise report
                    finish(index, report, 'calculation', attempt)
    return reports
