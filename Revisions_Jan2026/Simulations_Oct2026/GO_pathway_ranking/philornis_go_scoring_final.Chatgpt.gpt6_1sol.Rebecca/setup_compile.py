#!/usr/bin/env python3
"""Prepare, validate, back up, and compile independent ratings; never assign scores."""
import argparse, collections, csv, hashlib, json, random, shutil, zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'coord' / 'all_pathways.txt'
RUNS = ROOT / 'runs'
STATE = ROOT / 'coord' / 'run_state.json'

def terms():
    lines = SOURCE.read_text().splitlines()
    entries = [line.split(' ', 1) for line in lines[1:]]
    assert len(entries) == 1008 and len({x[0] for x in entries}) == 1008
    return lines[0], entries

def setup(source):
    shutil.copyfile(source, SOURCE)
    header, entries = terms()
    template = (ROOT / 'coord' / 'rater_template.txt').read_text()
    RUNS.mkdir(exist_ok=True)
    with (RUNS / 'run_manifest.tsv').open('w') as f:
        f.write('run\tshuffle_seed\tinput_sha256\n')
        for n in range(1, 31):
            d = RUNS / f'run_{n:03d}'
            d.mkdir(exist_ok=True)
            shuffled = [a + ' ' + b for a, b in entries]
            random.Random(1000 + n).shuffle(shuffled)
            p = d / 'input_terms.txt'
            p.write_text(header + '\n' + '\n'.join(shuffled) + '\n')
            prompt = template.format(INPUT=p, REVIEW=d/'review.md', SCORES=d/'scores.tsv', RUNDIR=d)
            (d/'prompt.txt').write_text(prompt)
            f.write(f'{n:03d}\t{1000+n}\t{hashlib.sha256(p.read_bytes()).hexdigest()}\n')
    STATE.write_text(json.dumps({f'{n:03d}': {'attempts': 0, 'token_use': 'unavailable', 'agents': []} for n in range(1,31)}, indent=2))

def validate(n):
    expected = {x[0] for x in terms()[1]}
    d = RUNS / f'run_{n:03d}'
    with (d/'scores.tsv').open() as f:
        reader = csv.reader(f, delimiter='\t')
        assert next(reader) == ['pathway','score'], 'wrong header'
        rows = list(reader)
    assert len(rows) == 1008, f'{len(rows)} rows'
    assert all(len(r) == 2 for r in rows), 'wrong column count'
    ids = [r[0] for r in rows]
    assert set(ids) == expected and len(set(ids)) == 1008, 'missing, extra, or repeated IDs'
    assert all(r[1] in {'1','2','3','4','5'} for r in rows), 'noninteger/out of range score'
    assert (d/'review.md').is_file() and (d/'review.md').stat().st_size > 0, 'missing review'
    return dict(rows)

def log():
    state = json.loads(STATE.read_text())
    with (ROOT/'run_log.tsv').open('w') as f:
        f.write('run\tshuffle_seed\ttoken_use\tre_run\tattempts\tstatus\tscore_1\tscore_2\tscore_3\tscore_4\tscore_5\n')
        for n in range(1,31):
            s = state[f'{n:03d}']
            try:
                scores = validate(n)
                counts = collections.Counter(scores.values())
                status = 'passed'
                count_fields = [str(counts[str(i)]) for i in range(1,6)]
            except (OSError, AssertionError, StopIteration) as e:
                status = 'pending' if s['attempts'] == 0 else 'incomplete'
                count_fields = ['']*5
            fields = [f'run_{n:03d}', str(1000+n), s['token_use'], str(s['attempts']>1).lower(), str(s['attempts']), status] + count_fields
            f.write('\t'.join(fields)+'\n')

def backup(batch, final=False):
    log()
    dest = ROOT/'backups'/('finch_ratings_final.zip' if final else f'finch_ratings_batch_{batch:02d}.zip')
    with zipfile.ZipFile(dest, 'w', zipfile.ZIP_DEFLATED) as z:
        for p in sorted(RUNS.rglob('*')):
            if p.is_file() and (final or p.name != 'input_terms.txt'):
                z.write(p, p.relative_to(ROOT))
        for p in [ROOT/'setup_compile.py', ROOT/'run_log.tsv', STATE, ROOT/'coord'/'rater_template.txt', SOURCE]:
            z.write(p,p.relative_to(ROOT))
        if (ROOT/'score_matrix.tsv').exists():
            z.write(ROOT/'score_matrix.tsv','score_matrix.tsv')
    print(dest)

def compile_all():
    data = [validate(n) for n in range(1,31)]
    with (ROOT/'score_matrix.tsv').open('w') as f:
        f.write('pathway\tterm_name\t'+'\t'.join(f'run_{n:03d}' for n in range(1,31))+'\n')
        for go, name in terms()[1]:
            f.write('\t'.join([go,name]+[r[go] for r in data])+'\n')
    backup(6,True)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action',choices=['setup','validate','attempt','backup','compile'])
    p.add_argument('value', nargs='?')
    p.add_argument('--agent',default='')
    a=p.parse_args()
    if a.action=='setup': setup(a.value)
    elif a.action=='validate':
        for n in map(int,a.value.split(',')):
            validate(n)
            print(f'run_{n:03d}: passed')
        log()
    elif a.action=='attempt':
        state=json.loads(STATE.read_text()); s=state[f'{int(a.value):03d}']; s['attempts']+=1; s['agents'].append(a.agent); STATE.write_text(json.dumps(state,indent=2))
    elif a.action=='backup': backup(int(a.value))
    elif a.action=='compile': compile_all()
