project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)

  repeat {
    if (file.exists(file.path(current, "README.md")) &&
        file.exists(file.path(current, "environment.yml"))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find the elderly project root.", call. = FALSE)
    }
    current <- parent
  }
}

set_project_root <- function() {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  start <- if (length(script_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE))
  } else {
    getwd()
  }

  root <- project_root(start)
  setwd(root)
  invisible(root)
}

project_file <- function(...) {
  file.path(getwd(), ...)
}

data_file <- function(...) {
  project_file("data_bitbucket", ...)
}

source_script <- function(...) {
  source(project_file("R", ...), chdir = FALSE)
}

require_data_dir <- function() {
  data_dir <- data_file()
  if (!dir.exists(data_dir)) {
    stop(
      paste0(
        "Missing data directory: ", data_dir, "\n",
        "Download the study data and place it at data_bitbucket/ in the repository root."
      ),
      call. = FALSE
    )
  }
  invisible(data_dir)
}
