I have everything needed. Here are the complete notes.

---

# Visual-Fidelity Comparison Recipe — `doc-preview` (RTF/DOCX → PDF → similarity)

Scope note up front: `doc-preview` is a **preview viewer** product. Its shipped fidelity tooling (`tools/fidelity/`) is a **page-count / engine-identity gate**, not a pixel comparator. The actual **visual similarity engine** (rasterize + pixel-diff + gridline detector) lives in the Phase-0 spike at `scratch/spike/harness.py`. Both are documented below; your "render both to PDF and score" goal maps onto the spike harness plus the shared conversion path. There is **no SSIM anywhere in this repo** (`grep -i ssim|skimage|structural_similarity|opencv` returns nothing) — the metric is a custom masked pixel-area diff.

---

## 1. Conversion engine(s) — RTF→PDF and DOCX→PDF

**One engine for both formats: LibreOffice (`soffice`), driven over HTTP by Gotenberg.** There is no pandoc, no Node/JS conversion library, and no per-format engine — LibreOffice imports both `.rtf` and `.docx` and both go through the identical Gotenberg route. (`pdf.js` / `pdfjs-dist` appears in the repo, but only for *browser rasterization of the already-converted PDF* — it is not part of RTF/DOCX→PDF conversion.)

- **Route:** `POST /forms/libreoffice/convert`, multipart form, field name **`files`**, returns `application/pdf`.
  - Confirmed in `tools/fidelity/gate.py` (`_LO_CONVERT = "/forms/libreoffice/convert"`), `scratch/spike/convert_corpus.py`, and `server/packages/doc_preview_api/src/doc_preview_api/converter/gotenberg.py`.
  - **Gotenberg routes by file extension**, so the uploaded part's filename must end in `.rtf` or `.docx` (the server sends `document.<fmt>`; the spike sends the source's real name).
- **Pinned engine identity** (from `conversion-service/VERIFICATION.md` + `Dockerfile`): Gotenberg **8.34.0**, LibreOffice **26.2.4.2**, Debian 13 trixie, glibc 2.41, CPython 3.13.5, OpenJDK 21 (pdftk). Upstream digest `gotenberg/gotenberg@sha256:3c23aeb3...db2a` (= `gotenberg/gotenberg:8-libreoffice`).
- Runtime flags this repo deploys: `--api-timeout=120s --libreoffice-restart-after=10 --libreoffice-auto-start=true --libreoffice-max-queue-size=10`. One conversion in flight per replica (scale horizontally).
- Fonts (current build): **substitute-only** — LibreOffice's bundled **Liberation** family (metric-compatible for Arial/Times New Roman/Courier New), **no Microsoft core fonts**, no mscorefonts EULA needed for internal use. Real fonts are deferred (OQ-2); `preseed/mscorefonts.seed` is retained to re-enable them. Aptos has no substitute (surfaces a `font_substitution` notice).

### Runtime dependencies, precisely
- **Conversion:** a running Gotenberg-with-LibreOffice HTTP service. Nothing else.
- **Comparison (Python):** `numpy`, `Pillow` (PIL), `pypdf`, and (for `harness.py` as-is) `pyyaml`. Plus **poppler-utils CLIs invoked as subprocesses**: `pdftoppm` (150-DPI raster), `pdftotext` (text-set), `pdfinfo`/`pdffonts` used by other probes. PyMuPDF is deliberately rejected (AGPL); poppler is called as a CLI to keep licensing clean.

### State of THIS machine (verified live)
- `soffice`/`libreoffice`: **NOT on PATH**, no `/Applications/LibreOffice.app`. `qpdf`, ImageMagick: not installed.
- `pdftoppm`, `pdftotext`, `pdfinfo`, `pdffonts`: **present** (`/opt/homebrew/bin`). `python3` with `numpy 2.3.3`, `PIL`, `pypdf 6.14.2`: present. `docker`: present.
- The task said "no running docker service yet" — **that is not currently true**. A conversion container is **already up and healthy**:
  - `docker ps` → `doc-preview-conversion-1` … `Up 10 days (healthy)` … `0.0.0.0:3000->3000/tcp`.
  - `curl http://localhost:3000/health` → `{"status":"up","details":{"libreoffice":{"status":"up",...}}}`
  - `curl http://localhost:3000/version` → `8.34.0`; `soffice.bin --version` inside → `LibreOffice 26.2.4.2 620(Build:2)`; `uname -m` → `aarch64`.
  - Local images present: `gotenberg/gotenberg:8-libreoffice`, `gotenberg/gotenberg:8`, `doc-preview-conversion:latest`, `doc-preview-conversion:spike`.

**If you need to start one yourself** (either the daemon is down or you want a clean instance):
```bash
# Simplest, matches the repo engine (image already cached locally):
docker run --rm -p 3000:3000 gotenberg/gotenberg:8-libreoffice \
  gotenberg --api-timeout=120s --libreoffice-restart-after=10 \
  --libreoffice-auto-start=true --libreoffice-max-queue-size=10
# or the org build:
#   cd /Users/mstackhouse/Documents/repos/doc-preview/conversion-service
#   docker build -t doc-preview-conversion:dev .
#   docker run --rm -p 3000:3000 doc-preview-conversion:dev
curl -sf http://localhost:3000/health   # wait for "up"
```
**No-Docker fallback** (native LibreOffice CLI; only if you can't use Gotenberg):
```bash
brew install --cask libreoffice          # installs /Applications/LibreOffice.app
SOFFICE=/Applications/LibreOffice.app/Contents/MacOS/soffice
"$SOFFICE" --headless --convert-to pdf --outdir OUTDIR INPUT.docx
# poppler (already installed here): brew install poppler
```
Note: the CLI path bypasses Gotenberg's queueing/normalizer; you must apply the RTF normalizer (section 3) yourself before converting RTF.

---

## 2. Fidelity comparison method

### The metric (from `scratch/spike/harness.py`)
Not SSIM. Per document, per page, five checks combine into one pass/fail:

1. **Page-count** — `pypdf` page count of preview vs golden must be equal.
2. **Per-page text-content-set equality** — `pdftotext -f n -l n <pdf> -` per page; whitespace-stripped into a **sorted character multiset** (`tuple(sorted("".join(text.split())))`). Compares *which characters are on page n*, tolerating cell-extraction order/run-boundary noise but catching a row/footnote/header shifted across a page boundary or missing content. Symbol-font PUA bullets are normalized (`PUA_EQUIV = {"\uf0b7":"\u2022", "\uf0a7":"\u25aa", "\uf0d8":"\u27a2"}`).
3. **Geometry** — per-page `MediaBox` (w,h) within **±1 pt**.
4. **Rasterized pixel diff @ 150 DPI grayscale** (the visual score) — see below.
5. **Gridline-shift structural detector** — see below.

**Rasterization:** `pdftoppm -f N -l N -r 150 -gray -png <pdf> <prefix>` → `<prefix>-<page>.png` (renamed to `<prefix>.png`). DPI is hard-coded **150**, grayscale.

**Page alignment (two stages):**
- *Size:* `load_aligned()` opens both as grayscale `int16`, pads both to the common `max(H)×max(W)` with **255 (white)**, top-left anchored.
- *Registration:* `register(a,b)` computes ink-projection profiles (`pixels < 128`), then finds the best **integer global translation** (dx,dy) by cross-correlation over a **±40 px** window on each axis, and `np.roll`s `b` onto `a`. This absorbs the few-px block offset between Word and LibreOffice; relative offsets *within* the page survive. No scaling, no rotation.

**Pixel-diff score (`pixel_diff`):**
```
raw = |a - b| > 32                       # INTENSITY_CUTOFF = 32
ink_a = dilate(dilate(a < 160));  ink_b = dilate(dilate(b < 160))   # 3x3 binary dilation, x2
significant = raw & ~(ink_a & ink_b)     # ignore anti-aliasing/hinting near shared ink
pct = 100.0 * significant.sum() / significant.size     # % of page area that differs
```
Returns `(pct, mask)`; the mask feeds a red-overlay heatmap for the 3 worst pages.

**Gridline-shift detector (`gridline_shift`):** after registration, on each axis take ink projection, mark "gridlines" as runs where ink `> 0.6 × extent`, cluster into centers, match each golden gridline to the nearest preview gridline. Returns the max residual displacement in px; a gridline with **no counterpart within 25 px** returns **999** (structural failure).

**Thresholds (in `harness.py`, keyed by document class):**
```python
PIXEL_THRESHOLD = {"rtf-prepaginated": 1.0, "rtf-dynamic": 1.0,
                   "rtf-dynamic-real": 1.0, "docx": 2.5}   # percent of page area
INTENSITY_CUTOFF = 32
GRIDLINE_SHIFT_PX = 2        # a shift > 2 px fails regardless of pixel %
```
So: **≤1.0 %** differing area for all RTF classes, **≤2.5 %** for DOCX; **default 2.5 %** for unknown classes. Gridline shift must be **≤ 2 px**.

**Overall verdict (per document):**
```python
row["pass"] = (page_count_match and not text_bad and not geom_bad
               and not pixel_fail_pages and not shift_fail_pages)
```
i.e. ALL of: equal page count AND every page text-set-equal AND every page geometry within 1 pt AND no page over the pixel threshold AND no page over the gridline shift. Raster diff is capped at the first 60 pages in the spike.

### The shipped gate vs the golden manifest (`tools/fidelity/gate.py`, `golden-manifest.json`)
The **production** gate is deliberately narrower — it validates **pagination only**, not pixels:
- `gate.py` converts each manifest doc through the *running* service and asserts `pypdf` page count equals `expected_pages`. Pure functions `count_pages(pdf_bytes)` and `compare_page_counts(expected, actual)` are unit-tested without a service. A missing/failed conversion counts as a mismatch (not a skip); an **empty comparison set is a FAILURE, not a pass**.
- It **refuses to run unless `golden-manifest.json` has `"seeded": true`** (a real boolean) and a non-empty `documents` list; otherwise exit code 2.
- **`golden-manifest.json` is currently a SCAFFOLD**: `"seeded": false`, `engine.libreoffice_version: "REPLACE_WITH_VALIDATED_LIBREOFFICE_VERSION"`, all `expected_pages: 0`, `fonts_required: ["Liberation Sans","Liberation Serif","Liberation Mono"]`. It carries **no pixel thresholds** — the only numeric gate here is per-document `expected_pages`. To seed: build the image, copy `libreoffice_version` from `/usr/share/doc-preview/build-manifest.json`, set each `expected_pages` to Microsoft Word's page count, set `seeded=true`.
- Run: `python tools/fidelity/gate.py --gotenberg http://localhost:3000 --golden tools/fidelity/golden-manifest.json --corpus corpus`.

### `engine_skew.py` explained
This is a **drift alarm on the engine's identity**, not a visual comparison. Pagination reproducibility is keyed on *exactly* which LibreOffice + fonts render the doc, so:
- `diff_manifests(golden, runtime)` compares the deployed `build-manifest.json` (`/usr/share/doc-preview/build-manifest.json`, mounted as a ConfigMap at `/etc/doc-preview/...`) against the golden baseline. It flags drift if **`libreoffice_version` differs exactly**, or if any **`fonts_required`** entry is missing from `runtime.fonts` (extra runtime fonts are fine). Non-empty reasons → exit 1.
- Run: `python tools/fidelity/engine_skew.py --golden tools/fidelity/golden-manifest.json --runtime /etc/doc-preview/build-manifest.json`.
- **Scope limit (stated in the file):** the `build-manifest.json` is generated in and COPY'd *frozen from* the Docker `engine` stage, so `engine_skew.py` validates only the **engine-stage identity**. It **cannot** see font/lib drift introduced by the runtime/base image — that is caught separately in CI by the `fc-list` diff and the `gate.py` page-count gate, both run against the *built runtime image* (see `VERIFICATION.md` §B-4/§B-5). So `engine_skew.py` is the "did the layout engine change under us?" tripwire; `gate.py` is the "did pagination actually move?" tripwire.

---

## 3. RTF preprocessing — `rtf_normalize.py` and the fonttbl fix

RTF must be normalized **before** it is sent to LibreOffice, or it renders unfaithfully. The identical normalizer exists in both the spike (`scratch/spike/rtf_normalize.py`) and production (`server/packages/doc_preview_api/src/doc_preview_api/rtf_normalize.py`); the server converter calls it automatically for `fmt == "rtf"`. Public signature:

```python
def normalize_rtf(data: bytes) -> bytes   # bytes in, bytes out; deterministic, content-preserving
```

It fixes two characterized **LibreOffice 26.2 RTF-import defects** (minimal repros documented in `scratch/plans/research-findings.md` "Phase 0 memo"):

1. **fonttbl whitespace (the headline bug).** If any whitespace sits at **`{\fonttbl ...}` depth-1** (after an entry's closing `}` or before the table's closing `}`), LibreOffice **silently fails to resolve ANY font name for the whole document**; Word tolerates it. Repro: `{\fonttbl {\f1 Courier New;}}` resolves ✓ but `{\fonttbl {\f1 Courier New;} }` fails ✗. This matters acutely for clinical TLFs — the pharmaRTF/Tplyr/r2rtf generators emit exactly this indented form, and losing Courier reflows every monospace table (measured ~24→23 px pitch @150 DPI, spilling 1–2 rows/page). Fix: `_fonttbl_span()` finds the balanced `{\fonttbl…}` group, then all `space/tab/CR/LF` at **depth==1** inside it are dropped (entry-internal whitespace preserved). `scratch/spike/fonttbl_fix_test.py` is the verification harness that proved page-count parity restored after this transform.
2. **`\clpadfX0` ignored-value semantics.** Per RTF spec, cell-padding unit-type `0` means "null — ignore the companion `\clpadX`"; Word ignores it, LibreOffice applies it, inflating row heights. Real TLF cells carry `\clpadft0\clpadt80 \clpadfb0\clpadb80`. Fix: `re.sub(rb"\\clpadf([tb])0\s*\\clpad\1\d+ ?", b"", data)` strips only the **vertical** (top/bottom) pairs; horizontal (l/r) left untouched (no divergence observed).

Not fixable in markup (documented, out of scope): LibreOffice ignores `\clNoWrap` and wraps long cell text.

**Practical rule for your project:** run every `.rtf` through `normalize_rtf()` before POSTing to Gotenberg. DOCX needs no preprocessing.

---

## 4. Provenance of `scratch/spike/preview/real-docx-t-14-*.pdf`

These are **LibreOffice/Gotenberg renderings of user-supplied real CDISC-pilot table DOCX**, used as Phase-0 fidelity fixtures. Concretely:

- **What they are:** the 9 `real-docx-t-14-*.pdf` files are the *preview* (engine-rendered) outputs for the 9 private source DOCX in `corpus/private/documents/docx/`:
  `T_14_1_1_1`, `T_14_1_1_2`, `T_14_1_1_3`, `T_14_1_2_1`, `T_14_2_1_1`, `T_14_2_1_2`, `T_14_3_1_1`, `T_14_3_1_7`, `T_14_3_4_1_B` (`.docx`).
- **Who produced the previews:** `pdfinfo`/XMP on `real-docx-t-14-1-1-1.pdf` → `Producer: LibreOffice 26.2.4.2 (AARCH64)`, `CreatorTool: Writer`, `CreationDate D:20260712...` — i.e. generated by a locally-run (arm64) Gotenberg+LibreOffice during the spike. Page size **A4 landscape (841.89 × 595.30 pt)**, 1 page.
- **Their goldens:** `corpus/private/golden/real-docx-t-14-*.wordmac-user.pdf`, produced by **Microsoft Word** — `Author(Mike Stackhouse) /Creator(Microsoft Word)`, `CreationDate D:20260712...`, A4 (842 × 595.25 pt). Manifest `golden_method: word-mac-16.110-export (user-confirmed, same machine)`.
- **Source DOCX metadata (real study files, Word-native):** `docProps/core.xml` → `lastModifiedBy: michael.collins`, `revision 9`, `created 2017-02-28T11:18:00Z`, `modified 2026-05-07`; `docProps/app.xml` → `Application: Microsoft Office Word`, `Template: Normal.dotm`, heading style **"Titre"** (French). These are authored/edited in desktop Word, not generated by any code in this repo.
- **Registration:** `corpus/private/manifest.yaml`, id `real-docx-t-14-1-1-1` … `real-docx-t-14-3-4-1-b`, `class: docx`, note *"user-supplied real table DOCX (A4 landscape by section definition; companion golden matches A4)."*
- **Confidentiality:** everything under `corpus/private/` is **gitignored and never committed** (`git check-ignore` confirms; `git ls-files` shows no `t_14`/`14-` entries). README: *"Real documents are never committed (user directive, 2026-07-11)."*

**Answer to "was a DOCX replication of these exact CDISC pilot tables already attempted?"** — The `T_14_x_x_x` naming is the CDISC pilot TLF numbering (Table 14.x…), and these DOCX are real renditions of those tables used here purely as fidelity fixtures. **Important nuance:** `doc-preview` did **not generate** these DOCX from RTF — it only *consumes* them. There is **no code in this repo that converts/replicates RTF into DOCX**; the DOCX are Word-native, supplied by the user (michael.collins/Mike Stackhouse). So what exists is: (a) the CDISC-pilot table DOCX themselves, (b) their Word-export golden PDFs, (c) their LibreOffice-rendered preview PDFs — i.e. a complete before/after fidelity dataset, but not a generator.

Alongside them sit **32 real RTF TLFs** `corpus/private/documents/rtf/14-*.rtf` (ids `real-rtf-14-1-01` … `real-rtf-14-7-04`, `class: rtf-dynamic-real`), described as *"user-supplied real TLF (pharmaRTF/Tplyr idiom: no explicit \page, landscape section, dynamic pagination)"*, each with its own Word golden and preview. These are the same CDISC-pilot output family in RTF form — directly relevant to your original-RTF-vs-new-DOCX goal, since here the RTFs and DOCX are the *same tables* rendered two ways.

(The public corpus additionally contains a CDISC-derived Define-XML fixture, `documents/datasets/tdf_adam_define.xml`, noted as "PHUSE TestDataFactory ADaM define.xml … derived from the CDISC pilot" — unrelated to the T_14 tables but confirming the CDISC-pilot lineage of this corpus.)

---

## 5. Minimal, copy-pasteable harness (run from any repo)

Prereqs on this machine are already satisfied: Gotenberg is up at `http://localhost:3000`, poppler + numpy/PIL/pypdf are installed. If Gotenberg is down, start it with the `docker run` in section 1.

### Option A — self-contained script (no dependency on doc-preview being importable)
This inlines the RTF normalizer and the exact pixel-diff/registration/gridline logic from `harness.py`, so it runs anywhere `numpy`, `Pillow`, `pypdf`, `httpx`, and `pdftoppm`/`pdftotext` are available.

```python
#!/usr/bin/env python3
"""fidelity_check.py — render an RTF and a DOCX to PDF via Gotenberg, then score
per-page visual similarity. Verdict logic + thresholds mirror doc-preview's
scratch/spike/harness.py. Usage:
    python fidelity_check.py original.rtf regenerated.docx
"""
import re, sys, subprocess, tempfile
from pathlib import Path
import numpy as np, httpx
from PIL import Image
from pypdf import PdfReader

GOTENBERG = "http://localhost:3000/forms/libreoffice/convert"
DPI = 150
INTENSITY_CUTOFF = 32
GRIDLINE_SHIFT_PX = 2
PIXEL_THRESHOLD_PCT = 2.5   # RTF-vs-DOCX cross-engine; use 1.0 for RTF-vs-RTF pre-paginated

# ---- RTF normalizer (verbatim from doc-preview/rtf_normalize.py) -------------
def _fonttbl_span(data: bytes):
    start = data.find(rb"{\fonttbl")
    if start < 0: return None
    depth = 0
    for i in range(start, len(data)):
        ch = data[i:i+1]
        if ch == b"{": depth += 1
        elif ch == b"}":
            depth -= 1
            if depth == 0: return start, i+1
    return None

def normalize_rtf(data: bytes) -> bytes:
    span = _fonttbl_span(data)
    if span:
        s, e = span
        out = bytearray(); depth = 0
        for ch in data[s:e]:
            b = bytes([ch])
            if b == b"{": depth += 1
            elif b == b"}": depth -= 1
            if depth == 1 and b in (b" ", b"\t", b"\r", b"\n"): continue
            out.append(ch)
        data = data[:s] + bytes(out) + data[e:]
    return re.sub(rb"\\clpadf([tb])0\s*\\clpad\1\d+ ?", b"", data)

# ---- (a)(b) convert to PDF through Gotenberg/LibreOffice ---------------------
def to_pdf(src: Path, dest: Path):
    payload = src.read_bytes()
    if src.suffix.lower() == ".rtf":
        payload = normalize_rtf(payload)
    # filename extension drives Gotenberg's engine routing
    files = {"files": (src.name, payload, "application/octet-stream")}
    r = httpx.post(GOTENBERG, files=files, timeout=300.0)
    r.raise_for_status()
    dest.write_bytes(r.content)

# ---- rasterize + pixel diff (verbatim logic from harness.py) -----------------
def rasterize(pdf: Path, page: int, prefix: Path) -> Path:
    out = Path(str(prefix) + ".png")
    subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-r", str(DPI),
                    "-gray", "-png", str(pdf), str(prefix)], check=True,
                   capture_output=True)
    produced = sorted(prefix.parent.glob(prefix.name + "-*.png"))
    produced[0].rename(out)
    return out

def load_aligned(a_png, b_png):
    a = np.asarray(Image.open(a_png).convert("L"), dtype=np.int16)
    b = np.asarray(Image.open(b_png).convert("L"), dtype=np.int16)
    h, w = max(a.shape[0], b.shape[0]), max(a.shape[1], b.shape[1])
    pad = lambda m: np.pad(m, ((0, h-m.shape[0]), (0, w-m.shape[1])), constant_values=255)
    return pad(a), pad(b)

def dilate(m):
    o = m.copy()
    o[1:, :] |= m[:-1, :]; o[:-1, :] |= m[1:, :]
    o[:, 1:] |= m[:, :-1]; o[:, :-1] |= m[:, 1:]
    return o

def register(a, b):
    def best_shift(pa, pb, window=40):
        best, arg = -1.0, 0
        for s in range(-window, window+1):
            x, y = (pa[s:], pb[:len(pb)-s]) if s >= 0 else (pa[:len(pa)+s], pb[-s:])
            if len(x) == 0: continue
            sc = float(np.dot(x, y))
            if sc > best: best, arg = sc, s
        return arg
    ia, ib = (a < 128).astype(float), (b < 128).astype(float)
    dy = best_shift(ia.sum(1), ib.sum(1)); dx = best_shift(ia.sum(0), ib.sum(0))
    return np.roll(np.roll(b, dy, 0), dx, 1)

def pixel_diff(a_png, b_png):
    a, b = load_aligned(a_png, b_png); b = register(a, b)
    raw = np.abs(a - b) > INTENSITY_CUTOFF
    ink_a, ink_b = dilate(dilate(a < 160)), dilate(dilate(b < 160))
    sig = raw & ~(ink_a & ink_b)
    return 100.0 * sig.sum() / sig.size

def gridline_shift(a_png, b_png):
    a, b = load_aligned(a_png, b_png); b = register(a, b); worst = 0
    for axis in (0, 1):
        pa = (a < 128).sum(axis).astype(float); pb = (b < 128).sum(axis).astype(float)
        extent = a.shape[1-axis]
        peaks_a = np.where(pa > 0.6*extent)[0]; peaks_b = np.where(pb > 0.6*extent)[0]
        if len(peaks_a) == 0 and len(peaks_b) == 0: continue
        centers = lambda idx: ([] if len(idx) == 0 else
            [int(np.mean(g)) for g in np.split(idx, np.where(np.diff(idx) > 3)[0]+1)])
        ca, cb = centers(peaks_a), centers(peaks_b)
        for x in ca:
            n = min((abs(x-y) for y in cb), default=999); worst = max(worst, n if n <= 25 else 999)
        for y in cb:
            n = min((abs(y-x) for x in ca), default=999); worst = max(worst, n if n <= 25 else 999)
    return worst

def npages(pdf): return len(PdfReader(str(pdf)).pages)

# ---- (c) score + pass/fail ---------------------------------------------------
def main(rtf, docx):
    tmp = Path(tempfile.mkdtemp(prefix="fidelity-"))
    a_pdf, b_pdf = tmp/"original.pdf", tmp/"regenerated.pdf"
    to_pdf(Path(rtf), a_pdf); to_pdf(Path(docx), b_pdf)
    na, nb = npages(a_pdf), npages(b_pdf)
    print(f"pages: original={na}  regenerated={nb}")
    page_count_match = (na == nb)
    overall = page_count_match
    for i in range(1, min(na, nb)+1):
        ra = rasterize(a_pdf, i, tmp/f"a{i:04d}")
        rb = rasterize(b_pdf, i, tmp/f"b{i:04d}")
        pct = pixel_diff(ra, rb); shift = gridline_shift(ra, rb)
        ok = pct <= PIXEL_THRESHOLD_PCT and shift <= GRIDLINE_SHIFT_PX
        overall &= ok
        print(f"  page {i:>3}: diff={pct:6.3f}%  gridline_shift={shift:>3}px  "
              f"{'ok' if ok else 'FAIL'}")
    print(f"\nVERDICT: {'PASS — visually identical' if overall else 'FAIL'} "
          f"(threshold {PIXEL_THRESHOLD_PCT}% area, {GRIDLINE_SHIFT_PX}px gridline; "
          f"page-count {'match' if page_count_match else 'MISMATCH'})")
    sys.exit(0 if overall else 1)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
```
Run:
```bash
pip install numpy Pillow pypdf httpx        # poppler already present via `brew install poppler`
python fidelity_check.py original.rtf regenerated.docx
```
Add the text-content-set gate (recommended for clinical tables, catches row/footnote shifts pixels miss) by also comparing, per page, `tuple(sorted("".join(subprocess.run(["pdftotext","-f",str(i),"-l",str(i),PDF,"-"],capture_output=True,text=True).stdout.split())))` between the two PDFs — verdict then also requires text-set equality (exactly as `harness.py` does).

### Option B — reuse doc-preview's modules directly
The pure functions are importable with no service and no repo build:

```python
import sys
DP = "/Users/mstackhouse/Documents/repos/doc-preview"
sys.path.insert(0, f"{DP}/tools/fidelity")     # gate.py, engine_skew.py
sys.path.insert(0, f"{DP}/scratch/spike")      # harness.py, rtf_normalize.py

# pagination gate helpers (lazy-import pypdf/httpx internally; no service needed for these two)
from gate import count_pages, compare_page_counts

# the visual comparator (needs numpy, Pillow, pypdf, pyyaml installed; importing does NOT run main())
import harness   # exposes rasterize, load_aligned, register, pixel_diff, gridline_shift,
                 # page_dims, page_text, PIXEL_THRESHOLD, INTENSITY_CUTOFF, GRIDLINE_SHIFT_PX

from rtf_normalize import normalize_rtf
```
- Convert via Gotenberg exactly as the repo does — reuse `gate.convert(gotenberg_url, doc_path)` for DOCX; for RTF, normalize first (`gate.convert` does **not** normalize):
  ```python
  import httpx
  from rtf_normalize import normalize_rtf
  def convert_rtf(url, rtf_path):
      files = {"files": (rtf_path.name, normalize_rtf(rtf_path.read_bytes()), "application/octet-stream")}
      r = httpx.post(url.rstrip("/") + "/forms/libreoffice/convert", files=files, timeout=300.0)
      r.raise_for_status(); return r.content
  ```
- Then score with `harness.rasterize(...)` → `harness.pixel_diff(a_png, b_png)` → compare against `harness.PIXEL_THRESHOLD[cls]` and `harness.GRIDLINE_SHIFT_PX`, reproducing the verdict formula in section 2.

**Threshold guidance for your specific goal (original RTF vs newly generated DOCX):** you are comparing two *different source formats* rendered by the *same* engine, so expect slightly more registration noise than the repo's RTF-vs-RTF case. Use the DOCX threshold **2.5 %** area + **2 px** gridline as the pass bar, keep the **±1 pt geometry** and **page-count** gates hard, and add the **text-content-set** gate — that combination is exactly what `harness.py` treats as "visually identical" for the `docx` class. Tighten to **1.0 %** only if both are pre-paginated and you want the RTF-class strictness.