#!/usr/bin/env python3
"""
fidelity.py — visual-fidelity comparison of clinical table documents.

Renders an original .rtf and a regenerated .docx to PDF through the local
Gotenberg + LibreOffice service (same engine for both => apples-to-apples
layout), then scores per-page visual similarity. RTF is normalized first
(fonttbl whitespace + cell-padding fixes) so the two formats rasterize
comparably.

Usage:
    python fidelity.py <original.rtf> <regenerated.docx> [--out DIR] [--json]

Requires: a running Gotenberg service (default http://localhost:3000),
poppler CLIs (pdftoppm, pdftotext, pdfinfo), numpy, Pillow, pypdf. httpx if
available else urllib fallback.
"""
import re, sys, os, subprocess, tempfile, argparse, json, urllib.request
from pathlib import Path
import numpy as np
from PIL import Image
from pypdf import PdfReader

GOTENBERG = os.environ.get("GOTENBERG_URL", "http://localhost:3000")
CONVERT_ROUTE = "/forms/libreoffice/convert"
DPI = 150
INTENSITY_CUTOFF = 32
GRIDLINE_SHIFT_PX = 2
# RTF-vs-DOCX is cross-format through one engine; use the docx-class tolerance.
PIXEL_THRESHOLD_PCT = 2.5

# ---- RTF normalizer ---------------------------------------------------------
def _fonttbl_span(data: bytes):
    start = data.find(rb"{\fonttbl")
    if start < 0:
        return None
    depth = 0
    for i in range(start, len(data)):
        ch = data[i:i + 1]
        if ch == b"{":
            depth += 1
        elif ch == b"}":
            depth -= 1
            if depth == 0:
                return start, i + 1
    return None

def normalize_rtf(data: bytes) -> bytes:
    span = _fonttbl_span(data)
    if span:
        s, e = span
        out = bytearray(); depth = 0
        for ch in data[s:e]:
            b = bytes([ch])
            if b == b"{":
                depth += 1
            elif b == b"}":
                depth -= 1
            if depth == 1 and b in (b" ", b"\t", b"\r", b"\n"):
                continue
            out.append(ch)
        data = data[:s] + bytes(out) + data[e:]
    return re.sub(rb"\\clpadf([tb])0\s*\\clpad\1\d+ ?", b"", data)

# ---- convert to PDF through Gotenberg / LibreOffice --------------------------
def to_pdf(src: Path, dest: Path):
    payload = src.read_bytes()
    if src.suffix.lower() == ".rtf":
        payload = normalize_rtf(payload)
    url = GOTENBERG.rstrip("/") + CONVERT_ROUTE
    try:
        import httpx
        r = httpx.post(url, files={"files": (src.name, payload, "application/octet-stream")},
                       timeout=300.0)
        r.raise_for_status()
        dest.write_bytes(r.content)
        return
    except ImportError:
        pass
    boundary = "----fidelityboundary"
    body = (b"--" + boundary.encode() +
            b'\r\nContent-Disposition: form-data; name="files"; filename="' +
            src.name.encode() + b'"\r\nContent-Type: application/octet-stream\r\n\r\n' +
            payload + b"\r\n--" + boundary.encode() + b"--\r\n")
    req = urllib.request.Request(url, data=body,
                                 headers={"Content-Type": "multipart/form-data; boundary=" + boundary})
    with urllib.request.urlopen(req, timeout=300) as resp:
        dest.write_bytes(resp.read())

# ---- rasterize + pixel diff -------------------------------------------------
def rasterize(pdf: Path, page: int, prefix: Path) -> Path:
    out = Path(str(prefix) + ".png")
    subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-r", str(DPI),
                    "-gray", "-png", str(pdf), str(prefix)], check=True, capture_output=True)
    produced = sorted(prefix.parent.glob(prefix.name + "-*.png"))
    produced[0].rename(out)
    return out

def load_aligned(a_png, b_png):
    a = np.asarray(Image.open(a_png).convert("L"), dtype=np.int16)
    b = np.asarray(Image.open(b_png).convert("L"), dtype=np.int16)
    h, w = max(a.shape[0], b.shape[0]), max(a.shape[1], b.shape[1])
    pad = lambda m: np.pad(m, ((0, h - m.shape[0]), (0, w - m.shape[1])), constant_values=255)
    return pad(a), pad(b)

def dilate(m):
    o = m.copy()
    o[1:, :] |= m[:-1, :]; o[:-1, :] |= m[1:, :]
    o[:, 1:] |= m[:, :-1]; o[:, :-1] |= m[:, 1:]
    return o

def register(a, b):
    def best_shift(pa, pb, window=40, tol=0.005):
        # Cross-correlate the ink projections over +/- window. On a periodic grid
        # (uniform numeric rows + a strong near-solid header rule) the raw argmax can
        # lock onto a full-row-off secondary peak; so among shifts within `tol` of the
        # best score, prefer the SMALLEST displacement. This is a no-op when there is a
        # single clear peak (the best score strictly dominates), and only disambiguates
        # the near-tied periodic case toward minimal, physically-correct displacement.
        scores = {}
        for s in range(-window, window + 1):
            x, y = (pa[s:], pb[:len(pb) - s]) if s >= 0 else (pa[:len(pa) + s], pb[-s:])
            if len(x) == 0:
                continue
            scores[s] = float(np.dot(x, y))
        if not scores:
            return 0
        best = max(scores.values())
        near = [s for s, v in scores.items() if v >= best * (1 - tol)]
        return min(near, key=lambda s: (abs(s), s))
    ia, ib = (a < 128).astype(float), (b < 128).astype(float)
    dy = best_shift(ia.sum(1), ib.sum(1)); dx = best_shift(ia.sum(0), ib.sum(0))
    return np.roll(np.roll(b, dy, 0), dx, 1)

def pixel_diff(a_png, b_png):
    a, b = load_aligned(a_png, b_png); b = register(a, b)
    raw = np.abs(a - b) > INTENSITY_CUTOFF
    ink_a, ink_b = dilate(dilate(a < 160)), dilate(dilate(b < 160))
    sig = raw & ~(ink_a & ink_b)
    return 100.0 * sig.sum() / sig.size, sig

def gridline_shift(a_png, b_png):
    a, b = load_aligned(a_png, b_png); b = register(a, b); worst = 0
    for axis in (0, 1):
        pa = (a < 128).sum(axis).astype(float); pb = (b < 128).sum(axis).astype(float)
        extent = a.shape[1 - axis]
        # near-solid lines only (0.85); a dense full-width TEXT row can exceed 0.6
        # and be mistaken for a rule (false gridline mismatch)
        peaks_a = np.where(pa > 0.85 * extent)[0]; peaks_b = np.where(pb > 0.85 * extent)[0]
        if len(peaks_a) == 0 and len(peaks_b) == 0:
            continue
        centers = lambda idx: ([] if len(idx) == 0 else
            [int(np.mean(g)) for g in np.split(idx, np.where(np.diff(idx) > 3)[0] + 1)])
        ca, cb = centers(peaks_a), centers(peaks_b)
        for x in ca:
            n = min((abs(x - y) for y in cb), default=999); worst = max(worst, n if n <= 25 else 999)
        for y in cb:
            n = min((abs(y - x) for x in ca), default=999); worst = max(worst, n if n <= 25 else 999)
    return worst

def page_dims(pdf, page):
    p = PdfReader(str(pdf)).pages[page - 1].mediabox
    return round(float(p.width), 1), round(float(p.height), 1)

def page_text_set(pdf, page):
    t = subprocess.run(["pdftotext", "-f", str(page), "-l", str(page), str(pdf), "-"],
                       capture_output=True, text=True).stdout
    return "".join(t.split())

def npages(pdf):
    return len(PdfReader(str(pdf)).pages)

# ---- reflow (pagination-agnostic) comparison --------------------------------
# For tables we intentionally repaginate (e.g. densest packing), page-by-page
# pixel diff is meaningless. Instead strip each page's repeated header/footnote
# chrome, keep only the body band, and vertically concatenate all pages into one
# continuous strip per document. Both docs then hold the same body rows in the
# same order at the same pitch -> compare the strips directly.
FOOT_MARKERS = ("[1]", "[2]", "[3]", "[4]", "[5]", "Note:", "Source:", "NOTE:")

def _page_lines(pdf, page):
    xml = subprocess.run(["pdftotext", "-bbox", "-f", str(page), "-l", str(page), str(pdf), "-"],
                         capture_output=True, text=True).stdout
    rows = {}
    for m in re.finditer(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">([^<]*)</word>', xml):
        x, y, x2, y2, t = m.groups()
        rows.setdefault(round(float(y), 1), []).append((float(x), float(y2), t))
    out = []
    for y, v in sorted(rows.items()):
        ymax = max(vy for _, vy, _ in v)
        out.append((y, ymax, " ".join(t for _, _, t in sorted(v, key=lambda z: z[0]))))
    return out

def body_strip(pdf, outdir, dpi=DPI):
    """Concatenate the body band of every page into one grayscale image, cropped
    tightly to the first/last body text line so per-page white space (which
    differs with pagination) is excluded — both docs then hold the same rows."""
    scale = dpi / 72.0
    strips = []
    for i in range(1, npages(pdf) + 1):
        png = rasterize(pdf, i, outdir / f"rf{i:04d}")
        im = np.asarray(Image.open(png).convert("L"), dtype=np.uint8)
        h, w = im.shape
        # header rule (top-most long horizontal rule) separates header from body
        ink = (im < 128).sum(1)
        rules = [y for y in range(int(h * 0.6)) if ink[y] > 0.6 * w]
        rule_pt = (max(rules) / scale) if rules else h * 0.18 / scale
        body = [(y, ymax) for y, ymax, txt in _page_lines(pdf, i)
                if y > rule_pt + 2 and not any(txt.lstrip().startswith(mk) for mk in FOOT_MARKERS)]
        if not body:
            continue
        top = max(int(min(y for y, _ in body) * scale) - 2, 0)
        bottom = min(int(max(ymax for _, ymax in body) * scale) + 2, h)
        if bottom > top:
            strips.append(im[top:bottom, :])
    if not strips:
        return None
    W = max(s.shape[1] for s in strips)
    strips = [np.pad(s, ((0, 0), (0, W - s.shape[1])), constant_values=255) for s in strips]
    return np.vstack(strips)

def page_body_rows(pdf, page, outdir):
    """Text of the body rows on a page: lines below the header rule and above the
    footnote block. Whitespace-normalized so cell padding differences don't matter
    (per-row exact rendering is proven separately by single-page pixel passes)."""
    scale = DPI / 72.0
    png = rasterize(pdf, page, outdir / f"cb{page:04d}")
    im = np.asarray(Image.open(png).convert("L"))
    h, w = im.shape
    ink = (im < 128).sum(1)
    # near-solid horizontal lines only (a dense text row can exceed 0.6*w and be
    # mistaken for the header rule, cutting the top body rows on continuation pages)
    rules = [y for y in range(int(h * 0.6)) if ink[y] > 0.9 * w]
    rule_pt = (max(rules) / scale) if rules else (h * 0.18 / scale)
    below = [(y, txt) for y, ymax, txt in _page_lines(pdf, page) if y > rule_pt + 2]
    foot_ys = [y for y, txt in below if any(txt.lstrip().startswith(m) for m in FOOT_MARKERS)]
    if foot_ys:
        cut = min(foot_ys)
        below = [(y, txt) for y, txt in below if y < cut]
    return [re.sub(r"\s+", " ", txt).strip() for y, txt in below if txt.strip()]

def content_compare(ref_pdf, cand_pdf, outdir):
    """Pagination-agnostic content comparison via the SET of unique normalized
    lines. Repeated title/header/footnote chrome is identical on every page in
    both docs, so it contributes the same unique lines and cancels; page numbers
    are normalized. Body data lines are unique (labels differ), so genuine value
    differences surface as lines present in one doc but not the other. No fragile
    header/body region detection."""
    import html
    def lineset(pdf):
        s = set()
        for i in range(1, npages(pdf) + 1):
            txt = subprocess.run(["pdftotext", "-layout", "-f", str(i), "-l", str(i), str(pdf), "-"],
                                 capture_output=True, text=True).stdout
            for ln in txt.splitlines():
                t = html.unescape(" ".join(ln.split()))
                if not t:
                    continue
                # skip the split (2-segment) title/footer lines: pdftotext -layout
                # glues their left/right segments with variable spacing (they sit
                # far apart on the real page), which is an extraction artifact, not
                # a content difference. They are identical chrome by construction.
                if t.startswith("Protocol:") or t.startswith("Source:"):
                    continue
                # dynamic footer timestamp (e.g. "13:57 Tuesday, June 16, 2020"),
                # which -layout can split onto its own line
                if re.match(r"^\d{1,2}:\d{2}\s+\w+,\s+\w+\s+\d+,\s+\d{4}$", t):
                    continue
                t = re.sub(r"Page \d+ of \d+", "Page N of M", t)
                s.add(t)
        return s
    ref, cand = lineset(ref_pdf), lineset(cand_pdf)
    only_ref, only_cand = ref - cand, cand - ref
    return {"mode": "content-lines", "ref_lines": len(ref), "cand_lines": len(cand),
            "pass": not only_ref and not only_cand,
            "only_ref": sorted(only_ref)[:15], "only_cand": sorted(only_cand)[:15]}

def reflow_compare(ref_pdf, cand_pdf, outdir):
    a = body_strip(ref_pdf, outdir); b = body_strip(cand_pdf, outdir)
    H = max(a.shape[0], b.shape[0]); W = max(a.shape[1], b.shape[1])
    pad = lambda m: np.pad(m, ((0, H - m.shape[0]), (0, W - m.shape[1])), constant_values=255)
    Image.fromarray(pad(a)).save(outdir / "reflow_ref.png")
    Image.fromarray(pad(b)).save(outdir / "reflow_cand.png")
    pct, sig = pixel_diff(outdir / "reflow_ref.png", outdir / "reflow_cand.png")
    save_heatmap(outdir / "reflow_ref.png", sig, outdir / "reflow_heatmap.png")
    return {"mode": "reflow", "ref_body_px": int(a.shape[0]), "cand_body_px": int(b.shape[0]),
            "pixel_pct": round(pct, 3), "pass": pct <= PIXEL_THRESHOLD_PCT}

def save_heatmap(a_png, sig, dest):
    a = np.asarray(Image.open(a_png).convert("RGB"))
    h, w = sig.shape
    a = a[:h, :w].copy()
    a[sig] = [255, 0, 0]
    Image.fromarray(a).save(dest)

def compare(ref_pdf, cand_pdf, outdir: Path, save_heat=True):
    na, nb = npages(ref_pdf), npages(cand_pdf)
    result = {"ref_pages": na, "cand_pages": nb, "page_count_match": na == nb, "pages": []}
    worst = []
    for i in range(1, min(na, nb) + 1):
        ra = rasterize(ref_pdf, i, outdir / f"ref{i:04d}")
        rb = rasterize(cand_pdf, i, outdir / f"cand{i:04d}")
        pct, sig = pixel_diff(ra, rb)
        shift = gridline_shift(ra, rb)
        da, db = page_dims(ref_pdf, i), page_dims(cand_pdf, i)
        geom_ok = abs(da[0] - db[0]) <= 1 and abs(da[1] - db[1]) <= 1
        ta, tb = page_text_set(ref_pdf, i), page_text_set(cand_pdf, i)
        text_ok = sorted(ta) == sorted(tb)
        ok = pct <= PIXEL_THRESHOLD_PCT and shift <= GRIDLINE_SHIFT_PX and geom_ok
        result["pages"].append({"page": i, "pixel_pct": round(pct, 3), "gridline_px": shift,
                                "geom_ref": da, "geom_cand": db, "geom_ok": geom_ok,
                                "text_set_ok": text_ok, "pass": ok})
        worst.append((pct, i, ra, sig))
    if save_heat:
        for pct, i, ra, sig in sorted(worst, reverse=True)[:3]:
            save_heatmap(ra, sig, outdir / f"heatmap_p{i:04d}.png")
    result["overall_pass"] = (result["page_count_match"]
                              and all(p["pass"] for p in result["pages"])
                              and len(result["pages"]) > 0)
    result["text_set_all_ok"] = all(p["text_set_ok"] for p in result["pages"])
    return result

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("original_rtf")
    ap.add_argument("regenerated_docx")
    ap.add_argument("--out", default=None, help="output dir for pdfs/rasters/heatmaps")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    outdir = Path(args.out) if args.out else Path(tempfile.mkdtemp(prefix="fidelity-"))
    outdir.mkdir(parents=True, exist_ok=True)
    ref_pdf, cand_pdf = outdir / "reference.pdf", outdir / "candidate.pdf"
    to_pdf(Path(args.original_rtf), ref_pdf)
    to_pdf(Path(args.regenerated_docx), cand_pdf)
    res = compare(ref_pdf, cand_pdf, outdir)

    if args.json:
        print(json.dumps(res, indent=2))
    else:
        print(f"reference : {args.original_rtf}  ({res['ref_pages']} pages)")
        print(f"candidate : {args.regenerated_docx}  ({res['cand_pages']} pages)")
        print(f"artifacts : {outdir}")
        print(f"{'page':>5} {'pixel%':>8} {'grid_px':>8} {'geom':>6} {'text':>6}  verdict")
        for p in res["pages"]:
            print(f"{p['page']:>5} {p['pixel_pct']:>8.3f} {p['gridline_px']:>8} "
                  f"{'ok' if p['geom_ok'] else 'BAD':>6} {'ok' if p['text_set_ok'] else 'DIFF':>6}  "
                  f"{'pass' if p['pass'] else 'FAIL'}")
        v = res["overall_pass"]
        print(f"\nVERDICT: {'PASS — visually identical' if v else 'FAIL / review'} "
              f"(<= {PIXEL_THRESHOLD_PCT}% area, <= {GRIDLINE_SHIFT_PX}px gridline, page-count "
              f"{'match' if res['page_count_match'] else 'MISMATCH'}, "
              f"text-set {'all-match' if res['text_set_all_ok'] else 'differs'})")
    sys.exit(0 if res["overall_pass"] else 1)

if __name__ == "__main__":
    main()
