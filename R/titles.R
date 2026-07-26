# R/titles.R — read a table's titles/footnotes from data/titles.xlsx as clinify-ready lists

library(readxl)

#' Substitute a dynamic token in a title/footnote cell
#'
#' Expands the spreadsheet's `PAGE_FORMAT:`, `FILE_PATH:` and `DATE_FORMAT:`
#' tokens; plain text is returned unchanged.
#'
#' @param txt Cell text, possibly beginning with a token prefix.
#' @param source_path Source path substituted into a `FILE_PATH:` token.
#' @param date Date for a `DATE_FORMAT:` token; a character value is used
#'   verbatim (fidelity mode) and `NULL` uses the current time.
#' @return The substituted string.
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

#' Build the aligned text vector for one title/footnote line
#'
#' @param row One spreadsheet row.
#' @param type Line type, `"title"` or `"footnote"`.
#' @param source_path Source path passed to token substitution.
#' @param date Date passed to token substitution.
#' @return A length-2 vector for a split line, otherwise length 1.
.line_vec <- function(row, type, source_path, date) {
  t1 <- .sub_tokens(row$text1, source_path, date)
  if (identical(row$align, "split")) {
    t2 <- .sub_tokens(row$text2, source_path, date)
    return(c(t1, t2))
  }
  c(t1)                            # single element; alignment via .line_align()
}

#' Resolve the explicit alignment for one title/footnote line
#'
#' @param row One spreadsheet row.
#' @param type Line type, `"title"` or `"footnote"`.
#' @return `"left"` for a left-aligned title, otherwise `NA` (clinify default).
.line_align <- function(row, type) {
  if (identical(type, "title") && identical(row$align, "left")) "left" else NA_character_
}

#' Read a table's titles and footnotes from the spreadsheet
#'
#' @param table_number Table identifier to filter on.
#' @param path Path to the titles spreadsheet.
#' @param source_path Source path passed to token substitution.
#' @param date Date passed to token substitution.
#' @return A list with `titles` and `footnotes`, each a list of line vectors
#'   (`ls`) and per-line alignments (`align`).
read_titles <- function(table_number,
                        path = "data/titles.xlsx",
                        source_path = NULL,
                        date = NULL) {
  col_types <- c("text", "numeric", "text", "text", "text", "text",
                 "logical", "logical", "text")
  df <- readxl::read_excel(path, col_types = col_types)
  df <- df[df$table_number == table_number, , drop = FALSE]
  if (nrow(df) == 0) stop("No titles/footnotes found for table ", table_number)

  #' Build the line and alignment lists for one line type
  #'
  #' @param kind Line type to select, `"title"` or `"footnote"`.
  #' @return A list with `ls` (line vectors) and `align` (alignments).
  mk <- function(kind) {
    sub <- df[df$type == kind, , drop = FALSE]
    sub <- sub[order(sub$index), , drop = FALSE]
    list(ls    = lapply(seq_len(nrow(sub)), function(i) .line_vec(sub[i, ], kind, source_path, date)),
         align = vapply(seq_len(nrow(sub)), function(i) .line_align(sub[i, ], kind), character(1)))
  }
  list(titles = mk("title"), footnotes = mk("footnote"))
}

#' Extract a table's reference-RTF footer timestamp
#'
#' @param table_number Table identifier (matches `<id>.rtf` in `ref_dir`).
#' @param ref_dir Directory holding the reference RTFs.
#' @return The `"HH:MM Weekday, Month DD, YYYY"` string, or `NULL` if absent.
ref_timestamp <- function(table_number, ref_dir = REF_DIR) {
  f <- file.path(ref_dir, paste0(table_number, ".rtf"))
  if (!file.exists(f)) return(NULL)
  txt <- readChar(f, file.info(f)$size, useBytes = TRUE)
  m <- regmatches(txt, regexpr("[0-9]{1,2}:[0-9]{2} [A-Za-z]+, [A-Za-z]+ [0-9]{1,2}, [0-9]{4}", txt))
  if (length(m)) m[1] else NULL
}

#' Attach a table's titles, footnotes and row pitch to a clinify object
#'
#' @param x A clintable or clindoc.
#' @param table_number Table identifier.
#' @param path Path to the titles spreadsheet.
#' @param source_path Source path passed to token substitution.
#' @param date Footer date; `NULL` uses the table's reference-RTF timestamp.
#' @return `x` with titles, footnotes and row height applied.
add_titles_footnotes <- function(x, table_number,
                                 path = "data/titles.xlsx",
                                 source_path = NULL, date = NULL) {
  if (is.null(date)) date <- ref_timestamp(table_number)
  tf <- read_titles(table_number, path, source_path, date)
  # House row pitch: body 15.35pt, titles/footnotes 11.4pt; "atleast" floors single
  # lines compact and lets wrapped cells grow.
  x |>
    clinify::clin_add_titles(tf$titles$ls, align = tf$titles$align) |>
    clinify::clin_add_footnotes(tf$footnotes$ls, align = tf$footnotes$align) |>
    clinify::clin_row_height(body = 15.35, title = 11.4, footnote = 11.4,
                             rule = "atleast", unit = "pt")
}
