# Runs the full pipeline of the paper "Data Leakage and the Apparent Predictability
# of Exchange Rates with Statistical Learning: Evidence from Developed and BRICS
# Economies": data (02), models (03), evaluation (04, 06), tables (05) and the
# aggregate numbers quoted in the paper (07), for the leak-free design (main) and
# the four designs that reopen one leakage channel (window, nolags, filters,
# scaling), twelve economies and horizons 1 to 12.
#
# Usage (from this folder):
#   Rscript run_pipeline.R [options]
# Options:
#   --models-dir=DIR   folder for the estimated models (large .rds files),
#                      one subfolder per design. Default: models
#   --cores=N          cores used by 03 (default: all but one)
#   --designs=a,b      subset of main,window,nolags,filters,scaling (default: all)
#   --countries=a,b    subset of ca,dk,fr,gb,jp,se,br,cn,id,in,ru,za (default: all)
#   --tuning           also export the yearly tuning choices (export_tuning.R)
#   --dry-run          print the jobs without running them
#
# Each job runs in its own R process; its log is logs/<job>.log. A finished job
# writes logs/done/<job>, and is skipped when the pipeline is run again, so an
# interrupted run resumes where it stopped. Delete logs/done/<job> to run a job
# again. The pipeline stops at the first job that fails.

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
flag <- function(name) any(args == paste0("--", name))

all_designs <- c("main", "window", "nolags", "filters", "scaling")
sets <- list(cp = c("ca", "dk", "fr", "gb", "jp", "se"), brics = c("br", "cn", "id", "in", "ru", "za"))
models_dir <- opt("models-dir", "models")
cores      <- as.integer(opt("cores", max(1L, parallel::detectCores() - 1L)))
designs    <- strsplit(opt("designs", paste(all_designs, collapse = ",")), ",")[[1]]
countries  <- strsplit(opt("countries", paste(unlist(sets), collapse = ",")), ",")[[1]]
dry_run    <- flag("dry-run")
stopifnot(all(designs %in% all_designs), all(countries %in% unlist(sets)))

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
dir.create(file.path("logs", "done"), recursive = TRUE, showWarnings = FALSE)

jobs <- list()
add <- function(name, script, script_args = character(0), design = "main")
  jobs[[length(jobs) + 1]] <<- list(name = name, script = script, args = script_args, design = design)

# 1. Data. The main job of each country set downloads the series (IMF, FRED,
#    DBnomics; reads the OECD M1 file and fred_key.txt) and keeps a copy of the raw
#    data in outputs/raw_<set>.rds. The nolags and filters panels are built from
#    that copy. Window and scaling use the main panels.
used_sets <- names(sets)[vapply(sets, function(s) any(s %in% countries), logical(1))]
for (s in used_sets) add(paste0("02_data_", s, "_main"), "02_data_BRICS_plus_US_ref.r", s)
for (d in intersect(c("nolags", "filters"), designs))
  for (s in used_sets) add(paste0("02_data_", s, "_", d), "02_data_BRICS_plus_US_ref.r", s, d)

# 2. Models and evaluation, design by design.
for (d in designs) {
  md <- file.path(models_dir, d)
  for (c in countries) add(paste0("03_models_", d, "_", c), "03_models.R", c(c, md, cores), d)
  for (c in countries) add(paste0("04_analysis_", d, "_", c), "04_analysis_single.R", c(c, md), d)
  for (c in countries) add(paste0("06_mcs_", d, "_", c), "06_mcs_single.R", c(c, md), d)
  add(paste0("05_export_", d), "05_export_results.r", character(0), d)
  add(paste0("07_aggregate_", d), "aggregate_stats.R", character(0), d)
  if (flag("tuning"))
    for (c in countries) add(paste0("tuning_", d, "_", c), "export_tuning.R", c(d, c, md), d)
}
if (flag("tuning")) add("tuning_summary", "export_tuning.R", "summary")

cat(length(jobs), "jobs. Models in", normalizePath(models_dir, mustWork = FALSE), "\n")
for (j in jobs) {
  done <- file.path("logs", "done", j$name)
  if (file.exists(done)) { cat("[skip] ", j$name, "\n", sep = ""); next }
  cmd <- c("--no-save", "--no-restore", shQuote(j$script), vapply(j$args, shQuote, ""))
  if (dry_run) {
    cat(sprintf("[dry-run] RUN_VARIANT=%s  Rscript %s\n", j$design, paste(cmd[-(1:2)], collapse = " ")))
    next
  }
  dir.create(file.path(models_dir, j$design), recursive = TRUE, showWarnings = FALSE)
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " start ", j$name, "\n", sep = "")
  # The child R process inherits these variables (system2's env argument is not
  # used because it does not work on Windows). ONLY_HORIZONS empty: horizons 1 to
  # 12. MAX_FORECASTS empty: every forecast origin.
  Sys.setenv(RUN_VARIANT = j$design, ONLY_HORIZONS = "", MAX_FORECASTS = "")
  log <- file.path("logs", paste0(j$name, ".log"))
  status <- system2(rscript, cmd, stdout = log, stderr = log)
  if (!identical(as.integer(status), 0L))
    stop("Job ", j$name, " failed. See logs/", j$name, ".log", call. = FALSE)
  writeLines(format(Sys.time()), done)
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ok    ", j$name, "\n", sep = "")
}
cat("Finished.\n")
