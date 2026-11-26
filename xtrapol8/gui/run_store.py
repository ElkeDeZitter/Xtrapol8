"""Persistent run registry. Scientific results are never stored or deleted here."""
import os
import sqlite3
from pathlib import Path

SCHEMA = '''CREATE TABLE IF NOT EXISTS runs (
 id TEXT PRIMARY KEY,
 name TEXT NOT NULL,
 directory TEXT NOT NULL,
 status TEXT NOT NULL,
 progress INTEGER NOT NULL DEFAULT 0,
 created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
 updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
)'''


def default_path():
    override = os.environ.get('XTRAPOL8_RUN_DB')
    if override:
        return Path(override).expanduser()
    return Path.home() / '.config' / 'xtrapol8' / 'run_manager.sqlite'


class RunStore:
    def __init__(self, path=None):
        self.path = Path(path) if path is not None else default_path()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self._connect() as con:
            con.execute(SCHEMA)

    def _connect(self):
        return sqlite3.connect(str(self.path), timeout=10)

    def all(self):
        with self._connect() as con:
            rows = con.execute('SELECT id,name,directory,status,progress FROM runs ORDER BY created_at, rowid').fetchall()
        return [dict(zip(('id','name','directory','status','progress_value'), row)) for row in rows]

    def save(self, entry):
        with self._connect() as con:
            con.execute('''INSERT INTO runs (id,name,directory,status,progress) VALUES (?,?,?,?,?)
                ON CONFLICT(id) DO UPDATE SET name=excluded.name, directory=excluded.directory,
                status=excluded.status, progress=excluded.progress, updated_at=CURRENT_TIMESTAMP''',
                (entry['id'], entry['name'], str(entry['directory']), entry['status'],
                 max(0, min(100, int(entry.get('progress_value', 0))))))

    def remove(self, run_id):
        with self._connect() as con:
            con.execute('DELETE FROM runs WHERE id=?', (run_id,))
