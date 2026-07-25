# R/titles.R
# -----------------------------------------------------------------------------
# Read title/footnote lines for a table out of data/titles.xlsx and return them
# as clinify-ready lists (for clin_add_titles() / clin_add_footnotes()).
#
# The spreadsheet follows the pilot's pharmaRTF convention:
#   columns: table_number, index, type(title|footnote), text1, text2, align, bold, italic, font
#   align: "split" (text1 left / text2 right), "left", "center"
#   dynamic tokens embedded in the text:
#     "PAGE_FORMAT: Page %s of %s"        -> "Page {PAGE} of {NUMPAGES}"  (clinify Word fields)
#     "FILE_PATH: Source: %s"             -> sprintf("Source: %s", source_path)
#     "DATE_FORMAT: <strftime pattern>"   -> formatted date (or a literal string)
#
# clinify alignment rules honoured:
#   * split line -> 2-element vector (left/right)
#   * title, single element  -> centered  ; title left -> duplicate trick c(t, t)
#   * footnote, single element -> left-aligned
#
# NOTE (feature request drafted for clinify): a spreadsheet->titles/footnotes
# reader like this belongs in clinify itself with pluggable token substitution;
# see notes/feature-requests/.
# -----------------------------------------------------------------------------

library(readxl)

.sub_tokens <- function(txt, source_path, date) {
  if (is.na(txt) || !nzchar(txt)) return("")
  if (startsWith(txt, "PAGE_FORMAT:")) {
    body <- trimws(sub("^PAGE_FORMAT:", "", txt))
    body <- sub("%s", "{PAGE}", body, fixed = TRUE)
    body <- sub("%s", "{NUMPAGES}", body, fixed = TRUE)
    return(body)
  }
  if (startsWith(txt, "FILE_PATH:")) {
    body <- trimws(sub("^FILE_PATH:", "", txt))   # e.g. "Source: %s"
    return(sprintf(body, if (is.null(source_path)) "" else source_path))
  }
  if (startsWith(txt, "DATE_FORMAT:")) {
    patt <- trimws(sub("^DATE_FORMAT:", "", txt))
    if (is.null(date))            return(format(Sys.time(), patt))
    if (is.character(date))       return(date)          # literal (fidelity mode)
    return(format(date, patt))
  }
  txt
}

.line_vec <- function(row, type, source_path, date) {
  t1 <- .sub_tokens(row$text1, source_path, date)
  if (identical(row$align, "split")) {
    t2 <- .sub_tokens(row$text2, source_path, date)
    return(c(t1, t2))
  }
  c(t1)                            # single element; alignment via .line_align()
}

# Per-line alignment for clin_add_titles(align=)/clin_add_footnotes(align=). NA = clinify
# default (title single -> center, footnote single -> left, 2-element -> split). A left title
# is now an explicit "left" (clinify #98) instead of the old duplicate-text merge trick.
.line_align <- function(row, type) {
  if (identical(type, "title") && identical(row$align, "left")) "left" else NA_character_
}

read_titles <- function(table_number,
                        path = "data/titles.xlsx",
                        source_path = NULL,
                        date = NULL) {
  col_types <- c("text", "numeric", "text", "text", "text", "text",
                 "logical", "logical", "text")
  df <- readxl::read_excel(path, col_types = col_types)
  df <- df[df$table_number == table_number, , drop = FALSE]
  if (nrow(df) == 0) stop("No titles/footnotes found for table ", table_number)

  mk <- function(kind) {
    sub <- df[df$type == kind, , drop = FALSE]
    sub <- sub[order(sub$index), , drop = FALSE]
    list(ls    = lapply(seq_len(nrow(sub)), function(i) .line_vec(sub[i, ], kind, source_path, date)),
         align = vapply(seq_len(nrow(sub)), function(i) .line_align(sub[i, ], kind), character(1)))
  }
  list(titles = mk("title"), footnotes = mk("footnote"))
}

# Extract the footer timestamp ("HH:MM Weekday, Month DD, YYYY") baked into a
# table's reference RTF, so a fidelity build reproduces it exactly.
ref_timestamp <- function(table_number, ref_dir = REF_DIR) {
  f <- file.path(ref_dir, paste0(table_number, ".rtf"))
  if (!file.exists(f)) return(NULL)
  txt <- readChar(f, file.info(f)$size, useBytes = TRUE)
  m <- regmatches(txt, regexpr("[0-9]{1,2}:[0-9]{2} [A-Za-z]+, [A-Za-z]+ [0-9]{1,2}, [0-9]{4}", txt))
  if (length(m)) m[1] else NULL
}

# Convenience wrapper: read titles for `table_number` and attach to a clintable
# or clindoc `x`. When `date` is NULL, the table's reference-RTF timestamp is used
# (fidelity mode).
add_titles_footnotes <- function(x, table_number,
                                 path = "data/titles.xlsx",
                                 source_path = NULL, date = NULL) {
  if (is.null(date)) date <- ref_timestamp(table_number)
  tf <- read_titles(table_number, path, source_path, date)
  # House row pitch (clinify #97 clin_row_height): body 15.35pt, titles/footnotes 11.4pt,
  # rule="atleast" (floor -> single lines compact, wrapped cells grow). One central call
  # replaces the former per-styler hrule/height_all and reaches the title/footnote blocks
  # (which the option stylers alone could not pitch).
  x |>
    clinify::clin_add_titles(tf$titles$ls, align = tf$titles$align) |>
    clinify::clin_add_footnotes(tf$footnotes$ls, align = tf$footnotes$align) |>
    clinify::clin_row_height(body = 15.35, title = 11.4, footnote = 11.4,
                             rule = "atleast", unit = "pt")
}
