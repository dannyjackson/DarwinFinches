#!/usr/bin/env python3
"""Prepare independent rater inputs, validate outputs, and compile without scoring.

Examples (from the extracted archive root):
  python3 setup_compile.py --root . setup
  python3 setup_compile.py --root . validate 1 2 3 4 5
  python3 setup_compile.py --root . compile
This program never assigns or changes scores. The restore-inputs command
reconstructs inputs omitted from batch backups without changing any results.
"""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
import random
import re
import shutil
import zipfile

PROMPT = """You are acting as an independent expert rater for a study of Darwin's finches and the
invasive avian vampire fly, Philornis downsi.

TASK
The file {INPUT} lists 1,008 Gene Ontology (GO) biological-process terms (one per line:
GO ID, a space, then the term name; the first line is a header). Score every term from
1 to 5 for how physiologically related it is to a potential host response to P. downsi.
The goal is an unbiased a priori list: if a term were later detected as under selection,
would we reasonably infer it is potentially associated with a response to the fly?

STEPS
1. Review the published literature on P. downsi and its effects on Galápagos land birds
   (Darwin's finches, mockingbirds and other passerines), using web search. Cover what the
   larvae do to hosts, and which host defenses (resistance, tolerance, behavioral) have been
   documented or proposed.
2. Base your scoring on that review.
3. Score every term.

SCALE
5 = directly captures a host mechanism documented in the Philornis literature
4 = core component of such a mechanism, or a strong a priori prediction from it
3 = plausible downstream or supporting process
2 = weak, indirect or speculative link
1 = no specific link

RULES
- Use only your literature review and general biological knowledge. Do not use memory,
  past conversations, or any information about results from genomic studies of these birds.
  Do not use transcriptomic, gene-expression or genomic studies of Galápagos birds (e.g., 
  Knutie et al. bioRxiv 2023.03.06.531275); if you encounter one, do not use its results.
- Judge each term on its own merits; do not aim for any particular distribution of scores.
- Score all 1,008 terms. Do not skip, merge or rename any.

INDEPENDENCE (required)
- Read only the file {INPUT}. Do not list, open or search any other file or directory on
  this machine, and do not use any project, document, memory or artifact tools. Other raters'
  work exists on this machine and you must not see it.
- Write only the two output files named below. If you write a helper script, put it in
  {RUNDIR} (never in a shared or temporary directory) and give it a name unique to this run.

OUTPUT
1. Write a brief summary of the mechanisms you identified and the key sources (with URLs)
   to {REVIEW}.
2. Write the scores to {SCORES} as a tab-separated file with the header "pathway<TAB>score"
   and exactly one row per GO ID in the input, using each GO ID exactly as given. Writing a
   short script that maps GO IDs to scores is fine, but every score must be your own judgment.
3. Reply with only one line: "done: <number of rows written>".
"""

def source(root):
    p = root / 'coordinator' / 'all_pathways.txt'
    if not p.exists():
        p = root / 'all_pathways.txt'
    lines = p.read_text(encoding='utf-8').splitlines()
    terms = [line.split(' ', 1) for line in lines[1:]]
    assert len(terms) == 1008, f'Expected 1008 source terms, got {len(terms)}'
    assert all(len(t) == 2 and re.fullmatch(r'GO:\d{7}', t[0]) and t[1] for t in terms)
    assert len({t[0] for t in terms}) == 1008, 'Duplicate source IDs'
    return p, lines[0], lines[1:], terms

def state_path(root):
    return root / 'runs' / 'run_state.json'

def state(root):
    return json.loads(state_path(root).read_text())

def save_state(root, data):
    state_path(root).write_text(json.dumps(data, indent=2) + '\n')

def setup(root):
    src, header, lines, _ = source(root)
    root.joinpath('runs').mkdir(parents=True, exist_ok=True)
    if state_path(root).exists():
        raise RuntimeError('Already set up; refusing to overwrite existing run state')
    data = {}
    with root.joinpath('runs/run_manifest.tsv').open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['run_number', 'run', 'shuffle_seed', 'input_sha256', 'prompt_sha256'])
        for n in range(1, 31):
            name = f'run_{n:03d}'
            d = root / 'runs' / name
            d.mkdir(exist_ok=True)
            shuffled = lines.copy()
            random.Random(1000 + n).shuffle(shuffled)
            inp = header + '\n' + '\n'.join(shuffled) + '\n'
            prompt = PROMPT.format(INPUT=d/'input_terms.txt', REVIEW=d/'review.md', SCORES=d/'scores.tsv', RUNDIR=d)
            d.joinpath('input_terms.txt').write_text(inp, encoding='utf-8')
            d.joinpath('prompt.txt').write_text(prompt, encoding='utf-8')
            w.writerow([n, name, 1000+n, hashlib.sha256(inp.encode()).hexdigest(), hashlib.sha256(prompt.encode()).hexdigest()])
            data[name] = {'run_number': n, 'shuffle_seed': 1000+n, 'attempts': [], 'status': 'pending', 'token_use': None, 'token_use_status': 'not_exposed_by_agent_interface'}
    save_state(root, data)
    print(json.dumps({'setup_runs': 30, 'terms_per_run': 1008, 'root': str(root), 'source_sha256': hashlib.sha256(src.read_bytes()).hexdigest()}))

def register(root, n, agent_id):
    data = state(root)
    r = data[f'run_{n:03d}']
    r['attempts'].append({'attempt': len(r['attempts'])+1, 'agent_id': agent_id, 'status': 'running'})
    r['status'] = 'running'
    save_state(root, data)

def restore_inputs(root):
    _, header, lines, _ = source(root)
    for n in range(1, 31):
        d = root/'runs'/f'run_{n:03d}'
        shuffled = lines.copy()
        random.Random(1000+n).shuffle(shuffled)
        expected = header+'\n'+'\n'.join(shuffled)+'\n'
        p = d/'input_terms.txt'
        if p.exists():
            assert p.read_text(encoding='utf-8') == expected, f'Existing input differs: {p}'
        else:
            p.write_text(expected, encoding='utf-8')
    print('All 30 input files are present and match their original shuffle seeds.')

def validate_one(root, n):
    _, _, _, terms = source(root)
    expected = {t[0] for t in terms}
    d = root/'runs'/f'run_{n:03d}'
    with (d/'scores.tsv').open(newline='', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        assert header == ['pathway', 'score'], f'Incorrect header: {header!r}'
        rows = list(reader)
    assert len(rows) == 1008, f'Expected 1008 rows, got {len(rows)}'
    assert all(len(row) == 2 for row in rows), 'Wrong number of columns'
    ids = [r[0] for r in rows]
    assert len(set(ids)) == 1008, 'Duplicate GO IDs'
    assert set(ids) == expected, f'ID mismatch: missing={len(expected-set(ids))}, extra={len(set(ids)-expected)}'
    assert all(re.fullmatch('[1-5]', r[1]) for r in rows), 'Scores must be integers 1 through 5'
    assert (d/'review.md').is_file() and (d/'review.md').stat().st_size > 0, 'Missing literature review'
    return {i: int(s) for i, s in rows}

def validate(root, numbers):
    data = state(root)
    result = []
    for n in numbers:
        name = f'run_{n:03d}'
        r = data[name]
        try:
            scores = validate_one(root, n)
            counts = collections.Counter(scores.values())
            r.update(status='passed', score_counts={str(k): counts[k] for k in range(1,6)}, score_sha256=hashlib.sha256((root/'runs'/name/'scores.tsv').read_bytes()).hexdigest())
            report = {'run': name, 'status': 'passed', 'rows': len(scores)}
        except Exception as exc:
            r['status'] = 'failed'
            report = {'run': name, 'status': 'failed', 'error': str(exc)}
        if r['attempts']:
            r['attempts'][-1].update(report)
        result.append(report)
    save_state(root, data)
    print(json.dumps(result))

def retry(root, n):
    data = state(root)
    name = f'run_{n:03d}'
    r = data[name]
    assert r['status'] == 'failed', 'Only failed runs may be reset'
    d = root/'runs'/name
    archive = root/'coordinator'/'failed_attempts'/name/f"attempt_{len(r['attempts']):02d}"
    archive.mkdir(parents=True, exist_ok=False)
    for p in d.iterdir():
        if p.name not in ('input_terms.txt', 'prompt.txt'):
            shutil.move(str(p), str(archive/p.name))
    r['status'] = 'pending_retry'
    save_state(root, data)

def write_run_log(root):
    data = state(root)
    with root.joinpath('run_log.tsv').open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['run_number', 'run', 'shuffle_seed', 'token_use', 'token_use_status', 're_run', 'attempt_count', 'status', 'score_1_count', 'score_2_count', 'score_3_count', 'score_4_count', 'score_5_count'])
        for n in range(1,31):
            name=f'run_{n:03d}'
            r=data[name]
            w.writerow([n, name, r['shuffle_seed'], r['token_use'] if r['token_use'] is not None else 'NA', r['token_use_status'], 'yes' if len(r['attempts'])>1 else 'no', len(r['attempts']), r['status'], *[r.get('score_counts',{}).get(str(k),'') for k in range(1,6)]])

def compile_matrix(root):
    data=state(root)
    assert all(r['status']=='passed' for r in data.values()), 'All 30 runs must pass before compilation'
    results=[validate_one(root,n) for n in range(1,31)]
    _,_,_,terms=source(root)
    with root.joinpath('score_matrix.tsv').open('w', newline='') as f:
        w=csv.writer(f,delimiter='\t')
        w.writerow(['pathway','term_name',*[f'run_{n:03d}' for n in range(1,31)]])
        for go, name in terms:
            w.writerow([go,name,*[r[go] for r in results]])
    write_run_log(root)
    print(json.dumps({'matrix_rows':1008,'score_columns':30,'total_scores':30240}))

def backup(root, batch, final=False):
    write_run_log(root)
    dest=root/'deliverables'
    dest.mkdir(exist_ok=True)
    name='finch_scoring_final.zip' if final else f'finch_scoring_backup_batch_{batch:02d}.zip'
    target=dest/name
    with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted((root/'runs').rglob('*')):
            if p.is_file() and (final or p.name!='input_terms.txt'):
                z.write(p,p.relative_to(root))
        z.write(Path(__file__), 'setup_compile.py')
        z.write(source(root)[0], 'all_pathways.txt')
        z.write(root/'run_log.tsv','run_log.tsv')
        if final:
            z.write(root/'score_matrix.tsv','score_matrix.tsv')
        failures=root/'coordinator'/'failed_attempts'
        if failures.exists():
            for p in sorted(failures.rglob('*')):
                if p.is_file(): z.write(p,p.relative_to(root/'coordinator'))
    with zipfile.ZipFile(target) as z:
        assert z.testzip() is None
        assert final or not any(Path(n).name=='input_terms.txt' for n in z.namelist())
    print(json.dumps({'archive':str(target),'bytes':target.stat().st_size,'complete_runs':sum(r['status']=='passed' for r in state(root).values())}))

def main():
    p=argparse.ArgumentParser()
    here=Path(__file__).resolve().parent
    p.add_argument('--root',type=Path,default=here.parent if here.name=='coordinator' else here)
    s=p.add_subparsers(dest='command',required=True)
    s.add_parser('setup')
    s.add_parser('restore-inputs')
    a=s.add_parser('register'); a.add_argument('run',type=int); a.add_argument('agent_id')
    a=s.add_parser('validate'); a.add_argument('runs',type=int,nargs='+')
    a=s.add_parser('retry'); a.add_argument('run',type=int)
    s.add_parser('compile')
    a=s.add_parser('backup'); a.add_argument('batch',type=int); a.add_argument('--final',action='store_true')
    args=p.parse_args(); root=args.root.resolve()
    if args.command=='setup': setup(root)
    elif args.command=='restore-inputs': restore_inputs(root)
    elif args.command=='register': register(root,args.run,args.agent_id)
    elif args.command=='validate': validate(root,args.runs)
    elif args.command=='retry': retry(root,args.run)
    elif args.command=='compile': compile_matrix(root)
    elif args.command=='backup': backup(root,args.batch,args.final)

if __name__=='__main__': main()
