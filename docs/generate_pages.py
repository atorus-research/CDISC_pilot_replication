#!/usr/bin/env python3
"""
generate_pages.py — generate the per-output and per-code Quarto pages for the docs site.

Reads the table programs (tables/t_*.R) for each output's id / title / population /
method, and the R/ helper + builder files, then writes:
  docs/outputs/<id>.qmd   — one page per output, with an embedded PDF viewer
  docs/code/<slug>.qmd    — one page per R file, description + rendered source

Idempotent: safe to re-run (it overwrites the generated pages). Wired as the Quarto
pre-render step in _quarto.yml so the site stays in sync with the code, and the pages
are also committed so a plain `quarto render` works without running this first.

Run:  python3 docs/generate_pages.py     (from anywhere — paths are resolved from __file__)
"""
import os
import re
import glob

DOCS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(DOCS)

# Branch the GitHub source/download links point at.
BRANCH = "master"

# Curated one-line descriptions for the shared R files (setup/config + builders).
R_DESCRIPTIONS = {
    "setup.R": "Session setup: library loads, data paths, `read_adam()` for reading ADaM / SDTM "
               "`.xpt` datasets, and the fidelity run-date.",
    "helpers.R": "Shared formatting helpers: fixed-width number formatting (`num_fmt`), wrapped "
                 "`(N=n)` arm labels, blank spacer rows, and p-value formatters (ANOVA / Fisher's "
                 "exact / chi-square).",
    "clinify_defaults.R": "clinify house style — US-Letter landscape, 1in margins, Courier New 10pt, "
                          "header bottom-rule, and the central body / title / footnote row pitch.",
    "tplyr2_defaults.R": "Project-wide tplyr2 rounding and format defaults applied to every spec.",
    "titles.R": "Titles / footnotes engine: reads `data/titles.xlsx`, substitutes page / file / date "
                "tokens, and applies per-line alignment through clinify.",
    "efficacy.R": "Shared ANCOVA efficacy builder (`build_efficacy_table()`): descriptive blocks plus "
                  "an `lm` / `emmeans` ANCOVA with pairwise contrasts, feeding the 14-3.x efficacy tables.",
    "ae.R": "Shared adverse-event builder (`build_ae_table()`): nested SOC / Preferred-Term "
            "distinct-subject counts with Fisher's exact p-values, feeding tables 14-5.01 and 14-5.02.",
}

SETUP_CONFIG = ["setup.R", "helpers.R", "clinify_defaults.R", "tplyr2_defaults.R", "titles.R"]
BUILDERS = ["efficacy.R", "ae.R"]

HEADER_RE = re.compile(
    r"Table\s+(?P<id>[\d.\-]+):\s*(?P<title>.*?)\s*\(Population:\s*(?P<pop>[^)]*)\)\s*"
    r"Produces:\s*\S+\s*(?P<method>.*)",
    re.DOTALL,
)


def parse_table_header(path):
    """Parse a table program's top comment block into id / title / population / method."""
    lines = []
    with open(path) as fh:
        for line in fh:
            s = line.rstrip("\n")
            if s.startswith("#"):
                lines.append(s.lstrip("#").strip())
            else:
                break
    # lines[0] is the filename comment; join the rest so a wrapped title is captured.
    joined = " ".join(lines[1:])
    m = HEADER_RE.search(joined)
    fname = os.path.basename(path)
    fid = fname[2:-2]  # t_14_1_01.R -> 14_1_01
    a, b, c = fid.split("_")
    dash_id = f"{a}-{b}.{c}"  # -> 14-1.01, matches the docx name
    if not m:
        return {"id": dash_id, "title": dash_id, "pop": "", "method": "", "file": fname}
    return {
        "id": dash_id,
        "title": m.group("title").strip(),
        "pop": m.group("pop").strip(),
        "method": m.group("method").strip(),
        "file": fname,
    }


def write_output_page(meta):
    """Write docs/outputs/<id>.qmd with an embedded PDF viewer + downloads."""
    tid = meta["id"]
    slug_r = "t_" + tid.replace("-", "_").replace(".", "_") + ".R"
    body = f"""---
title: "{tid} — {yaml_escape(meta['title'])}"
description: "{yaml_escape(meta['title'])} (Population: {meta['pop']})"
categories: ["{meta['pop']}"]
---

**Population:** {meta['pop']}

{meta['method']}

::: {{.pdf-viewer}}
<iframe src="../pdf/{tid}.pdf" title="{tid} — {html_escape(meta['title'])}"></iframe>
:::

[Download PDF](../pdf/{tid}.pdf){{.btn .btn-outline-secondary .btn-sm}}
[Download DOCX](https://github.com/atorus-research/CDISC_pilot_replication/raw/{BRANCH}/outputs/{tid}.docx){{.btn .btn-outline-secondary .btn-sm}}
[View program on GitHub](https://github.com/atorus-research/CDISC_pilot_replication/blob/{BRANCH}/tables/{slug_r}){{.btn .btn-outline-secondary .btn-sm}}

See the program that builds this output: [`tables/{slug_r}`](../code/{slug_r[:-2]}.qmd).
"""
    with open(os.path.join(DOCS, "outputs", f"{tid}.qmd"), "w") as fh:
        fh.write(body)


def write_code_page(path, description, section, github_path):
    """Write docs/code/<slug>.qmd: description + the file's source in a highlighted block."""
    fname = os.path.basename(path)
    slug = fname[:-2] if fname.endswith(".R") else fname
    with open(path) as fh:
        source = fh.read().rstrip("\n")
    body = f"""---
title: "{fname}"
description: "{yaml_escape(description)}"
categories: ["{section}"]
---

*{description}*

[View on GitHub](https://github.com/atorus-research/CDISC_pilot_replication/blob/{BRANCH}/{github_path}){{.btn .btn-outline-secondary .btn-sm}}

```r
{source}
```
"""
    with open(os.path.join(DOCS, "code", f"{slug}.qmd"), "w") as fh:
        fh.write(body)


def yaml_escape(s):
    return s.replace('"', "'")


def html_escape(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def main():
    os.makedirs(os.path.join(DOCS, "outputs"), exist_ok=True)
    os.makedirs(os.path.join(DOCS, "code"), exist_ok=True)

    table_files = sorted(glob.glob(os.path.join(ROOT, "tables", "t_*.R")))
    metas = [parse_table_header(p) for p in table_files]
    for meta in metas:
        write_output_page(meta)

    # Shared R files: setup/config + builders.
    for fname in SETUP_CONFIG:
        write_code_page(os.path.join(ROOT, "R", fname), R_DESCRIPTIONS[fname],
                        "Setup & Config", f"R/{fname}")
    for fname in BUILDERS:
        write_code_page(os.path.join(ROOT, "R", fname), R_DESCRIPTIONS[fname],
                        "Shared Builders", f"R/{fname}")
    # Table programs.
    for p, meta in zip(table_files, metas):
        desc = f"Table {meta['id']}: {meta['title']} (Population: {meta['pop']})"
        write_code_page(p, desc, "Table Programs", f"tables/{os.path.basename(p)}")

    print(f"generated {len(metas)} output pages and "
          f"{len(SETUP_CONFIG) + len(BUILDERS) + len(metas)} code pages")


if __name__ == "__main__":
    main()
