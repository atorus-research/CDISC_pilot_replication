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

#' Read a table's titles and footnotes from the spreadsheet
#'
#' Returns clinify's long-form titles/footnotes contract (a data frame with
#' `type`, `text1`, `text2`, `align`), token-substituted and ready to hand to
#' `clin_add_titles()` / `clin_add_footnotes()` (each verb selects its own
#' `type` rows and ignores the surface it doesn't own). `text2` is populated
#' only for split lines; `align` is `"left"` for a left-aligned title and `NA`
#' (clinify default) otherwise, matching the reference presentation.
#'
#' @param table_number Table identifier to filter on.
#' @param path Path to the titles spreadsheet.
#' @param source_path Source path passed to token substitution.
#' @param date Date passed to token substitution.
#' @return A data frame with columns `type`, `text1`, `text2`, `align`.
read_titles <- function(table_number,
                        path = "data/titles.xlsx",
                        source_path = NULL,
                        date = NULL) {
  col_types <- c("text", "numeric", "text", "text", "text", "text",
                 "logical", "logical", "text")
  df <- readxl::read_excel(path, col_types = col_types)
  df <- df[df$table_number == table_number, , drop = FALSE]
  if (nrow(df) == 0) stop("No titles/footnotes found for table ", table_number)
  df <- df[order(match(df$type, c("title", "footnote")), df$index), , drop = FALSE]

  split <- !is.na(df$align) & df$align == "split"
  sub1 <- vapply(df$text1, .sub_tokens, character(1), source_path = source_path, date = date)
  sub2 <- vapply(seq_len(nrow(df)),
                 function(i) if (split[i]) .sub_tokens(df$text2[i], source_path, date) else NA_character_,
                 character(1))
  data.frame(
    type  = df$type,
    text1 = unname(sub1),
    text2 = sub2,
    align = ifelse(df$type == "title" & !is.na(df$align) & df$align == "left", "left", NA_character_),
    stringsAsFactors = FALSE
  )
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
#' @param header_leading Optional header line-spacing multiple for tables with
#'   multi-line header cells (e.g. wrapped arm labels). Passed through to the single
#'   `clin_row_height()` call below, because a second call would replace it rather
#'   than add to it.
#' @return `x` with titles, footnotes and row height applied.
add_titles_footnotes <- function(x, table_number,
                                 path = "data/titles.xlsx",
                                 source_path = NULL, date = NULL,
                                 header_leading = NULL) {
  if (is.null(date)) date <- ref_timestamp(table_number)
  tf <- read_titles(table_number, path, source_path, date)
  # One long-form data frame drives both verbs (clinify selects each surface's
  # rows and leaves the other alone). House row pitch: body 15.35pt,
  # titles/footnotes 11.4pt; "atleast" floors single lines compact and lets
  # wrapped cells grow.
  x |>
    clinify::clin_add_titles(tf) |>
    clinify::clin_add_footnotes(tf) |>
    clinify::clin_row_height(body = 15.35, title = 11.4, footnote = 11.4,
                             header_leading = header_leading,
                             rule = "atleast", unit = "pt")
}
