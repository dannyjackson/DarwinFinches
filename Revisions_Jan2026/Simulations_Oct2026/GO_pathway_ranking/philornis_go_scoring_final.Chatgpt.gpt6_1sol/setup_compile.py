from pathlib import Path
import argparse, csv, json, random, zipfile

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent if SCRIPT_DIR.name == 'coordinator' else SCRIPT_DIR
SOURCE = ROOT / 'upload/all_pathways.txt'
if not SOURCE.exists():
    SOURCE = ROOT / 'all_pathways.txt'
RUNS = ROOT / 'runs'
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
3. Reply with only one line: "done: <number of rows written>".'''

def source():
    lines = SOURCE.read_text().splitlines()
    terms = [line.split(' ', 1) for line in lines[1:]]
    assert len(terms) == 1008 and len({x[0] for x in terms}) == 1008
    return lines, terms

def setup():
    lines, _ = source()
    (ROOT/'coordinator').mkdir(exist_ok=True)
    RUNS.mkdir(exist_ok=True)
    with (RUNS / 'run_manifest.tsv').open('w') as f:
        f.write('run_number\trun\tshuffle_seed\n')
        for n in range(1, 31):
            d = RUNS / f'run_{n:03d}'
            d.mkdir(exist_ok=True)
            terms = lines[1:].copy()
            random.Random(1000+n).shuffle(terms)
            (d / 'input_terms.txt').write_text('\n'.join([lines[0]]+terms)+'\n')
            (d / 'prompt.txt').write_text(PROMPT.format(INPUT=d/'input_terms.txt', REVIEW=d/'review.md', SCORES=d/'scores.tsv', RUNDIR=d)+'\n')
            f.write(f'{n}\trun_{n:03d}\t{1000+n}\n')

def validate(n):
    _, terms = source()
    expected = {x[0] for x in terms}
    d = RUNS / f'run_{n:03d}'
    with (d/'scores.tsv').open(newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        assert next(reader) == ['pathway','score'], 'incorrect header'
        rows = list(reader)
    assert len(rows) == 1008, f'{len(rows)} rows'
    assert all(len(r) == 2 for r in rows), 'incorrect field count'
    ids = [r[0] for r in rows]
    assert len(set(ids)) == 1008 and set(ids) == expected, 'ID mismatch or duplication'
    assert all(r[1] in {'1','2','3','4','5'} for r in rows), 'invalid score'
    assert (d/'review.md').is_file() and (d/'review.md').stat().st_size, 'missing review'
    return {r[0]:int(r[1]) for r in rows}

def log():
    metadata_path = ROOT/'coordinator/run_metadata.json'
    meta = json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
    with (ROOT/'run_log.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n')
        w.writerow(['run_number','shuffle_seed','token_use','rerun','score_1','score_2','score_3','score_4','score_5'])
        for n in range(1,31):
            try: scores=validate(n)
            except (AssertionError,FileNotFoundError,StopIteration): continue
            m=meta.get(str(n),{})
            w.writerow([n,1000+n,m.get('token_use','unavailable'),m.get('rerun',False),*[list(scores.values()).count(s) for s in range(1,6)]])

def compile_all():
    _, terms = source()
    scores = [validate(n) for n in range(1,31)]
    with (ROOT/'score_matrix.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n')
        w.writerow(['pathway','term_name',*[f'run_{n:03d}' for n in range(1,31)]])
        for go,name in terms: w.writerow([go,name,*[s[go] for s in scores]])
    log()

def backup(batch, final=False):
    log()
    out=ROOT/'coordinator'/('scoring_final.zip' if final else f'backup_batch_{batch:02d}.zip')
    with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(RUNS.rglob('*')):
            if p.is_file() and (final or p.name != 'input_terms.txt'):
                z.write(p,p.relative_to(ROOT))
        z.write(ROOT/'run_log.tsv','run_log.tsv')
        z.write(Path(__file__),'setup_compile.py')
        for name in ['run_metadata.json','validation_log.tsv']:
            p=ROOT/'coordinator'/name
            if p.exists(): z.write(p,'coordinator/'+name)
        z.write(SOURCE,'all_pathways.txt')
        if final: z.write(ROOT/'score_matrix.tsv','score_matrix.tsv')
    print(out)

if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('action',choices=['setup','validate','compile','backup','final'])
    p.add_argument('numbers',nargs='*',type=int)
    a=p.parse_args()
    if a.action=='setup': setup()
    elif a.action=='validate':
        for n in a.numbers:
            try:
                validate(n)
                print(f'run_{n:03d}: PASS')
                with (ROOT/'coordinator/validation_log.tsv').open('a') as f: f.write(f'{n}\tPASS\n')
            except Exception as e:
                print(f'run_{n:03d}: FAIL: {e}')
                with (ROOT/'coordinator/validation_log.tsv').open('a') as f: f.write(f'{n}\tFAIL\t{e}\n')
    elif a.action=='compile': compile_all()
    elif a.action=='backup': backup(a.numbers[0])
    elif a.action=='final': compile_all(); backup(6,True)
