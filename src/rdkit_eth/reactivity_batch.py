"""Dataset-safe vertical reactivity: explicit jobs, successful cache, resumable rows."""
from __future__ import annotations

from contextlib import nullcontext, closing
from copy import deepcopy
from dataclasses import dataclass, replace
from numbers import Integral

from .errors import DescriptorUnavailableError, XTBCacheError
from .input import _effective_charge_uhf, _input_metadata, _xyz_from_mol
from .options import XTBOptions
from .persistence import (XTBDiskCache, ReactivityCheckpoint, _digest,
                          _implementation_identity)
from .reactivity import (run_reactivity, _empty_report, _scheme, REACTIVITY_UNITS,
                         REACTIVITY_PROTOCOL_VERSION)
from .reports import DESCRIPTOR_VERSION, SCHEMA_VERSION
from .runtime import _engine_identity, _parameter_context, _parameter_file
from .xtb import _resolve_options
from .electronic_descriptors import _fragments


@dataclass(frozen=True)
class ReactivityJob:
    """One indexed molecule and its per-molecule settings. IDs must be unique.

    mol=None may represent an upstream parsing failure and becomes an input-error
    row. A real molecule requires explicit H atoms and a finite 3D conformer.
    Missing/invalid spin choices or fragment definitions are configuration errors.
    """
    record_id: str
    mol: object
    options: XTBOptions | None = None
    ionized_uhf: int | None = None
    attached_uhf: int | None = None
    fragments: dict | None = None


def _prepare_jobs(jobs, options, population_scheme, interface, initialization):
    from rdkit import Chem
    _scheme(population_scheme)
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
        if not isinstance(job.record_id, str) or not job.record_id.strip():
            raise ValueError('record_id must be a nonempty string')
        ids.append(job.record_id)
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
        for spin in (job.ionized_uhf, job.attached_uhf):
            if spin is not None and (isinstance(spin, bool) or not isinstance(spin, Integral) or not 0 <= spin <= 2147483647):
                raise ValueError('adjacent spins must be nonnegative integers')
        # Fail on installation/capability errors even if a stored row exists.
        config = _resolve_options(replace(config, keep_overlap=config.keep_overlap or population_scheme != 'charges'))
        engine = _engine_identity('native', config.xtb_library)
        if not engine['sha256']:
            raise ValueError('persistent protocols require a fingerprintable libxtb file')
        if config.keep_overlap and not engine['orbital_basis']:
            raise DescriptorUnavailableError('population analysis requires the verified native basis extension')
        parameters = _parameter_file(config.method, config.xtb_library)
        if job.mol is not None and not isinstance(job.mol, Chem.Mol):
            raise TypeError('job molecule must be an RDKit Mol or None')
        mol = None if job.mol is None else Chem.Mol(job.mol)
        if mol is not None and mol.GetNumAtoms():
            charge, spin = _effective_charge_uhf(mol, config.charge, config.uhf)
            if spin and (job.ionized_uhf is None or job.attached_uhf is None):
                raise ValueError(f'{job.record_id}: open-shell references require explicit adjacent spins')
            _effective_charge_uhf(mol, charge+1, 1 if job.ionized_uhf is None else job.ionized_uhf)
            _effective_charge_uhf(mol, charge-1, 1 if job.attached_uhf is None else job.attached_uhf)
            groups = None if job.fragments is None else {
                k:v.tolist() for k,v in _fragments(job.fragments, mol.GetNumAtoms()).items()}
        else:
            groups = deepcopy(job.fragments)
        input_error = None
        try:
            if mol is None:
                raise ValueError('molecule is None (upstream parsing failed)')
            xyz = _xyz_from_mol(mol, config.conf_id)
            metadata = _input_metadata(mol, config.conf_id, charge, spin, xyz)
            config = replace(config, charge=charge, uhf=spin, conf_id=metadata['conformer_id'])
        except ValueError as exc:
            input_error = {'type': type(exc).__name__, 'message': str(exc)}
            # Invalid geometries cannot use the normal finite XYZ fingerprint.
            import hashlib
            metadata = {'invalid_input_sha256': None if mol is None else hashlib.sha256(mol.ToBinary()).hexdigest()}
        protocol = {'cache_schema': 1, 'protocol_version': REACTIVITY_PROTOCOL_VERSION,
            'schema_version': SCHEMA_VERSION, 'descriptor_version': DESCRIPTOR_VERSION,
            'implementation': implementation, 'input': metadata,
            'options': config.to_dict(), 'engine': engine, 'parameter_file': parameters,
            'parameter_context': _parameter_context(), 'interface': interface,
            'initialization': initialization, 'population_scheme': population_scheme,
            'ionized_uhf': None if job.ionized_uhf is None else int(job.ionized_uhf),
            'attached_uhf': None if job.attached_uhf is None else int(job.attached_uhf),
            'fragments': groups}
        prepared.append((job.record_id, mol, config, protocol, _digest(protocol), input_error))
    return prepared


def run_reactivity_batch(jobs, *, options=None, population_scheme='charges', interface='native',
                         initialization='reference', cache=None, checkpoint=None, retry_failed=False,
                         max_workers=1, timeout=None, isolate=False):
    """Run ordered ReactivityJob records with molecule-level checkpoint commits.

    Configuration and runtime errors propagate. Invalid molecular geometry/None
    and genuine engine calculation failures become explicit null-feature rows.
    Completed checkpoint rows resume by default, including failures. Set
    retry_failed=True to retry failures under exactly the same protocol, retaining
    attempt history. Changing the protocol requires a new checkpoint path.
    An interrupted triplet restarts from its reference; no wavefunction is pickled.
    """
    from .workers import validate_execution, isolated_map, reactivity, molecule_bytes
    from .errors import XTBTimeoutError, XTBWorkerError
    isolated = validate_execution(max_workers, timeout, isolate)
    if not isinstance(retry_failed, bool):
        raise ValueError('retry_failed must be Boolean')
    if cache is not None and not isinstance(cache, XTBDiskCache):
        raise TypeError('cache must be XTBDiskCache')
    prepared = _prepare_jobs(jobs, options, population_scheme, interface, initialization)
    if checkpoint is not None and not isinstance(checkpoint, ReactivityCheckpoint):
        checkpoint = ReactivityCheckpoint(checkpoint)
    if cache is not None and checkpoint is not None and cache.path == checkpoint.path:
        raise ValueError('cache and checkpoint must use separate files')
    manifest = {'schema_version': 1, 'rows': [{'record_id': p[0], 'key': p[4], 'protocol': p[3]} for p in prepared]}
    reports = [None] * len(prepared)
    tasks, pending = [], []

    def finish(index, report, source, attempt):
        identifier, mol, config, protocol, key, input_error = prepared[index]
        report['record_id'] = identifier
        report['provenance']['batch_protocol'] = protocol
        if source == 'calculation' and cache is not None and report['error'] is None:
            cache.set(key, report)
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
                if old['report']['error'] is None or not retry_failed:
                    report = deepcopy(old['report'])
                    report['execution']['source'] = 'checkpoint'
                    reports[index] = report
                    continue
            attempt = 1 if old is None else old['attempt'] + 1
            report = None if cache is None or input_error is not None else cache.get(key)
            source = 'cache' if report is not None else 'calculation'
            if report is not None and (not isinstance(report, dict) or
                    report.get('provenance', {}).get('batch_protocol') != protocol):
                raise XTBCacheError('cached report protocol does not match the requested calculation')
            if input_error is not None:
                source = 'input_validation'
                report = _empty_report(identifier, population_scheme, {})
                report['error'] = dict(input_error, stage='input')
                report['quality'] = {'status': 'error', 'flags': [
                    {'code': 'invalid_molecular_input', 'severity': 'error', 'detail': input_error['message']}]}
                report['missing'] = dict.fromkeys(REACTIVITY_UNITS, 'invalid_molecular_input')
                report['analysis_settings'].update(protocol_version=REACTIVITY_PROTOCOL_VERSION,
                    initialization=initialization, requested_options=config.to_dict(), interface=interface)
            elif report is None:
                kwargs = dict(options=config, ionized_uhf=protocol['ionized_uhf'],
                    attached_uhf=protocol['attached_uhf'], fragments=protocol['fragments'],
                    population_scheme=population_scheme, interface=interface,
                    initialization=initialization, return_errors=True)
                if isolated:
                    tasks.append((molecule_bytes(mol), kwargs)); pending.append((index, attempt))
                    continue
                report = run_reactivity(mol, **kwargs)
            finish(index, report, source, attempt)
        if isolated:
            with closing(isolated_map(reactivity, tasks, max_workers=max_workers, timeout=timeout)) as completed:
                for task_index, report in completed:
                    index, attempt = pending[task_index]
                    identifier, mol, config, protocol, key, input_error = prepared[index]
                    if isinstance(report, (XTBTimeoutError, XTBWorkerError)):
                        error = report
                        report = _empty_report(identifier, population_scheme, {})
                        report['error'] = {'type': type(error).__name__, 'message': str(error), 'stage': 'worker'}
                        report['quality'] = {'status': 'error', 'flags': [
                            {'code': 'isolated_execution_failed', 'severity': 'error', 'detail': str(error)}]}
                        report['missing'] = dict.fromkeys(REACTIVITY_UNITS, 'isolated_execution_failed')
                        report['analysis_settings'].update(protocol_version=REACTIVITY_PROTOCOL_VERSION,
                            initialization=initialization, requested_options=config.to_dict(), interface=interface)
                    elif isinstance(report, Exception):
                        raise report
                    finish(index, report, 'calculation', attempt)
    return reports
