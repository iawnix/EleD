"""Transactional JSON storage for numerical results and batch checkpoints."""
from __future__ import annotations

from contextlib import closing
from dataclasses import asdict, fields
import hashlib
import json
from pathlib import Path
import sqlite3
import sys

from .errors import XTBCacheError
from .result import XTBResult
from .runtime import _fingerprint

STORAGE_VERSION = 1
APPLICATION_ID = 1380202571


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def _digest(value):
    return hashlib.sha256(_json(value).encode()).hexdigest()


def _implementation_identity(interface='native'):
    import numpy
    from rdkit import rdBase
    root = Path(__file__).parent
    identity = {'python': list(sys.version_info[:2]), 'rdkit': rdBase.rdkitVersion,
                'numpy': numpy.__version__,
                'python_sources': {str(p.relative_to(root)): _fingerprint(p)
                                   for p in sorted(root.rglob('*.py'))}}
    if interface == 'compiled':
        from . import rdXTBTools
        if not hasattr(rdXTBTools.XTBResults, 'failureKind'):
            raise ImportError('compiled adapter is outdated; rebuild rdXTBTools for failure classification')
        identity['compiled_sha256'] = _fingerprint(Path(rdXTBTools.__file__))
    return identity


class _JSONStore:
    """One short-lived connection per operation; each write is transactional."""
    def __init__(self, path):
        self.path = Path(path).resolve()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with closing(self._connect()) as db, db:
            db.execute('BEGIN IMMEDIATE')
            app = db.execute('PRAGMA application_id').fetchone()[0]
            version = db.execute('PRAGMA user_version').fetchone()[0]
            tables = db.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
            if not tables and app == 0 and version == 0:
                db.execute('CREATE TABLE records (namespace TEXT, key TEXT, payload TEXT NOT NULL, '
                           'sha256 TEXT NOT NULL, PRIMARY KEY(namespace,key))')
                db.execute(f'PRAGMA application_id={APPLICATION_ID}')
                db.execute(f'PRAGMA user_version={STORAGE_VERSION}')
            elif app != APPLICATION_ID or version != STORAGE_VERSION:
                raise XTBCacheError('incompatible or unrelated SQLite database; choose another path')

    def _connect(self):
        return sqlite3.connect(self.path, timeout=30)

    @staticmethod
    def _decode(row):
        if row is None:
            return None
        payload, checksum = row
        if hashlib.sha256(payload.encode()).hexdigest() != checksum:
            raise XTBCacheError('persistent record checksum mismatch')
        try:
            value = json.loads(payload)
            _json(value)
            return value
        except (ValueError, TypeError) as exc:
            raise XTBCacheError('invalid persistent JSON record') from exc

    def read(self, namespace, key):
        with closing(self._connect()) as db:
            return self._decode(db.execute('SELECT payload,sha256 FROM records WHERE namespace=? AND key=?',
                                           (namespace, key)).fetchone())

    def write_many(self, records):
        encoded = []
        for namespace, key, value in records:
            # Preserve ordered descriptor/unit columns; keys use canonical _json separately.
            payload = json.dumps(value, separators=(',', ':'), allow_nan=False)
            encoded.append((namespace, key, payload, hashlib.sha256(payload.encode()).hexdigest()))
        with closing(self._connect()) as db, db:
            db.executemany('INSERT OR REPLACE INTO records VALUES (?,?,?,?)', encoded)

    def all(self, namespace):
        with closing(self._connect()) as db:
            return {key: self._decode((payload, checksum)) for key, payload, checksum in db.execute(
                'SELECT key,payload,sha256 FROM records WHERE namespace=? ORDER BY key', (namespace,))}


def _encode_result(value):
    if isinstance(value, XTBResult):
        if value.metadata.get('calculation_status') != 'success':
            raise ValueError('only successful single points may be cached')
        return {'kind': 'single_point', 'data': asdict(value)}
    if isinstance(value, dict) and value.get('report_type') == 'vertical_reactivity':
        if (value.get('error') is not None or value.get('quality', {}).get('status') == 'error'
                or set(value.get('states', {})) != {'reference', 'ionized', 'attached'}
                or any(s.get('status') != 'success' for s in value['states'].values())):
            raise ValueError('only complete successful reactivity reports may be cached')
        return {'kind': 'vertical_reactivity', 'data': value}
    raise TypeError('expected XTBResult or a successful reactivity report')


def _decode_result(value):
    try:
        if value['kind'] == 'single_point':
            data = value['data']
            def tuples(v):
                return tuple(tuples(x) for x in v) if isinstance(v, list) else v
            for field in fields(XTBResult):
                if isinstance(field.default, tuple):
                    data[field.name] = tuples(data[field.name])
            result = XTBResult(**data)
        elif value['kind'] == 'vertical_reactivity':
            result = value['data']
        else:
            raise ValueError('unknown stored result kind')
        _encode_result(result)
        return result
    except (ValueError, KeyError, TypeError) as exc:
        raise XTBCacheError('invalid persistent result schema') from exc


class XTBDiskCache:
    """Successful single-point/report cache. JSON only; no pickle or native handles.

    Works with run_xtb_batch, run_descriptor_batch and run_reactivity_batch.
    SQLite commits survive process interruption. Concurrent writers are safe,
    but concurrent identical misses may both calculate; no distributed work queue.
    Corruption or storage errors propagate rather than being treated as misses.
    """
    def __init__(self, path):
        self._store = _JSONStore(path)
        self.path = self._store.path

    def get(self, key):
        value = self._store.read('cache', key)
        return None if value is None else _decode_result(value)

    def set(self, key, value):
        self._store.write_many([('cache', key, _encode_result(value))])

    def get_or_compute(self, key, calculate):
        value = self.get(key)
        if value is None:
            value = calculate()
            self.set(key, value)
        return value

    def __len__(self):
        with closing(self._store._connect()) as db:
            return db.execute("SELECT COUNT(*) FROM records WHERE namespace='cache'").fetchone()[0]


class ReactivityCheckpoint:
    """Ordered dataset manifest and committed attempts, including failures.

    A run holds a nonblocking exclusive process lock; competing runs fail clearly.
    Resume requires the exact IDs, input order, geometries, protocols and engine
    fingerprints. Use a new path for a changed dataset and reuse the shared cache.
    """
    def __init__(self, path):
        self._store = _JSONStore(path)
        self.path = self._store.path

    def history(self, record_id):
        entries = self._store.all('attempts').values()
        return sorted((v for v in entries if v['record_id'] == record_id), key=lambda v: v['attempt'])

    def reports(self):
        manifest = self._store.read('manifest', 'dataset')
        rows = self._store.all('rows')
        return [] if manifest is None else [rows[r['record_id']]['report'] for r in manifest['rows']
                                            if r['record_id'] in rows]

    def _bind(self, manifest):
        old = self._store.read('manifest', 'dataset')
        if old is None:
            self._store.write_many([('manifest', 'dataset', manifest)])
        elif old != manifest:
            raise XTBCacheError('checkpoint dataset/protocol differs; use a new checkpoint path')

    def _save(self, identifier, key, report, attempt):
        row = {'record_id': identifier, 'key': key, 'report': report, 'attempt': attempt}
        self._store.write_many([('rows', identifier, row),
                                ('attempts', _json([identifier, attempt]), row)])

    def _locked(self):
        from contextlib import contextmanager
        import fcntl
        @contextmanager
        def lock():
            with self.path.with_suffix(self.path.suffix + '.lock').open('a') as handle:
                try:
                    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError as exc:
                    raise XTBCacheError('checkpoint already has an active writer') from exc
                try:
                    yield self
                finally:
                    fcntl.flock(handle, fcntl.LOCK_UN)
        return lock()
