> **RESOLVED** — filed as tplyr2 issues #13/#14 and fixed in merged PR #15
> (atorus-research/tplyr2). This project now uses the native behaviour/settings.

# tplyr2 FR: "less-than" percent + zero-count display for count layers

**Package:** tplyr2 · **Priority:** low/medium (nice-to-have; today handled by post-processing)

The legacy CDISC pilot displays use two common regulatory conventions that tplyr2
does not express natively:

1. **`<1%`** — a nonzero count whose percent rounds to 0 is shown as `<1` rather
   than `0`. (e.g. 1/254 → `1 ( <1%)`, not `1 (  0%)`.)
2. **Bare zero** — a category with a 0 count in an arm is shown as just ` 0`
   (no `( 0%)`).

We currently post-process the formatted strings (`fix_count_quirks()`), which is
brittle (regex on the display string). Proposals:

- An `f_str` modifier or `layer_settings` option for a **less-than threshold**,
  e.g. `pct_lt = 1` → render `<1` (and symmetrically `>99`), analogous to how AE
  Fisher p-values cap at `>.99`.
- A `zero_count_display = c("full","count_only","blank")` option so zero cells can
  render as ` 0` (count only) or blank, matching common SAP conventions.

Open question for the pilot rebuild: whether to reproduce these *legacy* display
quirks at all, or adopt tplyr2's cleaner defaults for the "modern" version. See
the decision noted in the project status.
