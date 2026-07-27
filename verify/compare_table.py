#!/usr/bin/env python3
"""
compare_table.py <table-id> [--json] [--out DIR]

Convenience driver: given a table id like "14-3.01", render the ORIGINAL rtf
(from the pristine master worktree) and the NEW docx (from ./outputs) to PDF and
score visual fidelity. Thin wrapper over fidelity.compare().

    python verify/compare_table.py 14-3.01
"""
import sys, argparse, tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fidelity as F

REPO = Path(__file__).resolve().parent.parent
ORIG_WORKTREE = REPO.parent / "CDISC_pilot_replication-original"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table_id")
    ap.add_argument("--out", default=None)
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--ref", default=None, help="override reference rtf path")
    ap.add_argument("--cand", default=None, help="override candidate docx path")
    ap.add_argument("--reflow", action="store_true",
                    help="pagination-agnostic body pixel comparison (repaginated tables)")
    ap.add_argument("--content", action="store_true",
                    help="pagination- and pitch-agnostic body CONTENT comparison (values/rows/order)")
    ap.add_argument("--header", action="store_true",
                    help="header block line-structure comparison (catches a changed header wrap, "
                         "which the pixel and text-set gates cannot see)")
    args = ap.parse_args()

    ref = Path(args.ref) if args.ref else ORIG_WORKTREE / "outputs" / f"{args.table_id}.rtf"
    cand = Path(args.cand) if args.cand else REPO / "outputs" / f"{args.table_id}.docx"
    if not ref.exists():
        sys.exit(f"reference rtf not found: {ref}")
    if not cand.exists():
        sys.exit(f"candidate docx not found: {cand}")

    outdir = Path(args.out) if args.out else REPO / "verify" / "_artifacts" / args.table_id
    outdir.mkdir(parents=True, exist_ok=True)
    ref_pdf, cand_pdf = outdir / "reference.pdf", outdir / "candidate.pdf"
    F.to_pdf(ref, ref_pdf)
    F.to_pdf(cand, cand_pdf)

    import json
    if args.header:
        r = F.header_compare(ref_pdf, cand_pdf, outdir)
        print(f"table     : {args.table_id}  (header block line structure)")
        print(f"header lines: ref={r['ref_lines']}  cand={r['cand_lines']}")
        for i, ref_l, cand_l in r["diffs"]:
            print(f"  line {i + 1}:\n      REF : {ref_l[:90]}\n      CAND: {cand_l[:90]}")
        for ln in r["only_ref"]:  print(f"  REF only : {ln[:90]}")
        for ln in r["only_cand"]: print(f"  CAND only: {ln[:90]}")
        print(f"\n{'PASS - header wraps identical' if r['pass'] else 'FAIL / header wrap differs'}")
        sys.exit(0 if r["pass"] else 1)

    if args.content:
        r = F.content_compare(ref_pdf, cand_pdf, outdir)
        print(f"table     : {args.table_id}  (content lines / pagination-agnostic)")
        print(f"unique body/data lines: ref={r['ref_lines']}  cand={r['cand_lines']}")
        for ln in r["only_ref"]:  print(f"  REF only : {ln[:90]}")
        for ln in r["only_cand"]: print(f"  CAND only: {ln[:90]}")
        print(f"\n{'PASS - content identical' if r['pass'] else 'FAIL / review content'}")
        sys.exit(0 if r["pass"] else 1)

    if args.reflow:
        r = F.reflow_compare(ref_pdf, cand_pdf, outdir)
        print(f"table     : {args.table_id}  (reflow / pagination-agnostic)")
        print(f"artifacts : {outdir}")
        print(f"body height px: ref={r['ref_body_px']}  cand={r['cand_body_px']}")
        print(f"body pixel diff: {r['pixel_pct']}%   "
              f"{'PASS - body visually identical' if r['pass'] else 'FAIL / review'}")
        sys.exit(0 if r["pass"] else 1)

    res = F.compare(ref_pdf, cand_pdf, outdir)
    if args.json:
        print(json.dumps(res, indent=2)); sys.exit(0 if res["overall_pass"] else 1)
    print(f"table     : {args.table_id}")
    print(f"reference : {ref}  ({res['ref_pages']} pages)")
    print(f"candidate : {cand}  ({res['cand_pages']} pages)")
    print(f"artifacts : {outdir}")
    print(f"{'page':>5} {'pixel%':>8} {'grid_px':>8} {'geom':>6} {'text':>6}  verdict")
    for p in res["pages"]:
        print(f"{p['page']:>5} {p['pixel_pct']:>8.3f} {p['gridline_px']:>8} "
              f"{'ok' if p['geom_ok'] else 'BAD':>6} {'ok' if p['text_set_ok'] else 'DIFF':>6}  "
              f"{'pass' if p['pass'] else 'FAIL'}")
    v = res["overall_pass"]
    print(f"\nVERDICT: {'PASS - visually identical' if v else 'FAIL / review'} "
          f"(page-count {'match' if res['page_count_match'] else 'MISMATCH'}, "
          f"text-set {'all-match' if res['text_set_all_ok'] else 'differs'})")
    sys.exit(0 if res["overall_pass"] else 1)

if __name__ == "__main__":
    main()
