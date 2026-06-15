# Shared small helpers for the elderly breast cancer signature workflow.

`%!in%` <- function(x, y) !(x %in% y)

progress_index <- function(i, label = "Processing item") {
  message(label, ": ", i)
}
