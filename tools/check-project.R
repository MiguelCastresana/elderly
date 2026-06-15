script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

setwd(repo_root)

required_files <- c(
  "README.md",
  "environment.yml",
  "R/main_signatures.R",
  "R/survival_analysis.R",
  "R/paths.R"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "), call. = FALSE)
}

r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
parse_results <- lapply(r_files, function(path) {
  tryCatch(
    {
      parse(path)
      NULL
    },
    error = function(err) {
      paste(path, conditionMessage(err), sep = ": ")
    }
  )
})

parse_errors <- Filter(Negate(is.null), parse_results)
if (length(parse_errors) > 0) {
  stop(paste(parse_errors, collapse = "\n"), call. = FALSE)
}

if (!dir.exists("data_bitbucket")) {
  message("data_bitbucket/ is not present. Full analysis run is skipped.")
}

message("Project structure and R syntax checks passed.")
