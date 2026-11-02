from pathlib import Path
import random, csv, json, zipfile, argparse

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'coordinator' / 'all_pathways.txt'
PROMPT = '''You are acting as an independent expert rater for a study of Darwin's finches and the
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
'''

def source():
    lines = SOURCE.read_text().splitlines()
    assert len(lines) == 1009
    ids = [line.split(' ', 1)[0] for line in lines[1:]]
    assert len(set(ids)) == 1008
    return lines, ids

def setup():
    lines, ids = source()
    runs = ROOT / 'runs'
    runs.mkdir(exist_ok=True)
    with (runs / 'run_manifest.tsv').open('w') as out:
        out.write('run_number\tshuffle_seed\n')
        for n in range(1, 31):
            d = runs / f'run_{n:03}'
            d.mkdir(exist_ok=True)
            terms = lines[1:].copy()
            random.Random(1000+n).shuffle(terms)
            (d/'input_terms.txt').write_text('\n'.join([lines[0]]+terms)+'\n')
            (d/'prompt.txt').write_text(PROMPT.format(INPUT=d/'input_terms.txt', REVIEW=d/'review.md', SCORES=d/'scores.tsv', RUNDIR=d))
            out.write(f'{n}\t{1000+n}\n')

def validate(n):
    _, ids = source()
    path = ROOT/'runs'/f'run_{n:03}'/'scores.tsv'
    with path.open() as f:
        rows=list(csv.reader(f, delimiter='\t'))
    assert rows[0] == ['pathway','score'], 'incorrect header'
    rows=rows[1:]
    assert len(rows)==1008, f'{len(rows)} rows'
    assert all(len(r)==2 for r in rows), 'column count'
    assert len({r[0] for r in rows})==1008, 'duplicate IDs'
    assert {r[0] for r in rows}==set(ids), 'ID mismatch'
    assert all(r[1] in ['1','2','3','4','5'] for r in rows), 'invalid scores'
    assert path.with_name('review.md').is_file(), 'missing review'
    return dict(rows)

def compile_results():
    lines, ids=source()
    scores=[validate(n) for n in range(1,31)]
    with (ROOT/'score_matrix.tsv').open('w') as f:
        writer=csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['pathway','term_name']+[f'run_{n:03}' for n in range(1,31)])
        for line, go in zip(lines[1:], ids):
            writer.writerow([go,line.split(' ',1)[1]]+[s[go] for s in scores])
    log_results()

def log_results():
    events_path=ROOT/'coordinator'/'events.json'
    events=json.loads(events_path.read_text()) if events_path.exists() else {}
    with (ROOT/'run_log.tsv').open('w') as f:
        writer=csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['run_number','shuffle_seed','token_use','re_run','score_1','score_2','score_3','score_4','score_5'])
        for n in range(1,31):
            try: s=validate(n)
            except (AssertionError,FileNotFoundError): continue
            event=events.get(str(n),{})
            writer.writerow([n,1000+n,event.get('token_use','unavailable'),event.get('re_run',False)]+[list(s.values()).count(str(k)) for k in range(1,6)])

def backup(name, final=False):
    log_results()
    dest=ROOT/name
    with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED) as z:
        for path in sorted((ROOT/'runs').rglob('*')):
            if path.is_file() and (final or path.name!='input_terms.txt'):
                z.write(path,path.relative_to(ROOT))
        for name in ['setup_compile.py','run_log.tsv','score_matrix.tsv','coordinator/events.json','coordinator/all_pathways.txt']:
            path=ROOT/name
            if path.exists(): z.write(path,path.relative_to(ROOT))
    print(dest)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('action',choices=['setup','validate','compile','backup','final'])
    p.add_argument('args',nargs='*')
    a=p.parse_args()
    if a.action=='setup': setup()
    elif a.action=='validate':
        for n in map(int,a.args):
            try: validate(n); print(f'run_{n:03}: PASS')
            except Exception as e: print(f'run_{n:03}: FAIL: {e}'); raise
    elif a.action=='compile': compile_results()
    else: backup(a.args[0],a.action=='final')
