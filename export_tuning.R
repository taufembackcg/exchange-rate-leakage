# Exporta as escolhas de tuning do re-tuning anual (M2 de reviews_v7.md) para
# CSV, para checar se o LASSO (e o Ridge) param no limite da grade de lambda
# (1e-4 a 10): com lambda = 10 todos os coeficientes sao zero e a previsao e a
# variacao media da janela (um drift).
#
# Uso:
#   Rscript export_tuning.R <variante> <pais> <pasta_dos_modelos>
#       le output_models_<pais>.rds e grava outputs/tuning/tuning_<variante>_<pais>.csv
#       (formato longo: uma linha por bloco anual, metodo e parametro)
#   Rscript export_tuning.R summary
#       junta os CSV em outputs/tuning/tuning_all.csv e grava o resumo do
#       LASSO e do Ridge em outputs/tuning/lambda_boundary_summary.csv e .txt
#
# Um processo R por arquivo de modelos (o .bat chama um por pais e variante),
# para nao acumular memoria.

suppressMessages(library(dplyr))
get_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", ca, value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  getwd()
}
setwd(get_script_dir())
out_dir <- file.path("outputs", "tuning")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Uso: export_tuning.R <variante> <pais> <pasta> | summary", call. = FALSE)

lambda_grid <- 10^seq(-4, 1, length.out = 50)   # grade do LASSO e do Ridge no 03

if (args[1] != "summary") {
  if (length(args) < 3) stop("Uso: export_tuning.R <variante> <pais> <pasta>", call. = FALSE)
  variant <- args[1]; cty <- args[2]
  f <- file.path(args[3], paste0("output_models_", cty, ".rds"))
  if (!file.exists(f)) { cat("sem arquivo, ignorado: ", f, "\n", sep = ""); quit(save = "no", status = 0) }
  models <- readRDS(f)
  rows <- list()
  for (key in names(models)) {
    parts <- regmatches(key, regexec("^[a-z]{2}_rolling_(dif|ecm)_(h[0-9]+)_form_(.+)$", key))[[1]]
    for (m in names(models[[key]])) {
      tb <- models[[key]][[m]]$tuning_by_year
      if (is.null(tb) || !nrow(tb)) next
      pars <- setdiff(names(tb), c("block", "first_rowIndex"))
      for (p in pars) {
        rows[[length(rows) + 1]] <- data.frame(
          variant = variant, country = cty, type = parts[2], horizon = parts[3], model = parts[4],
          method = m, block = tb$block, first_rowIndex = tb$first_rowIndex,
          param = p, value = as.character(tb[[p]]), stringsAsFactors = FALSE)
      }
    }
  }
  rm(models); invisible(gc())
  res <- if (length(rows)) do.call(rbind, rows) else data.frame()
  out <- file.path(out_dir, paste0("tuning_", variant, "_", cty, ".csv"))
  write.csv(res, out, row.names = FALSE)
  cat(variant, " ", cty, ": ", nrow(res), " linhas em ", out, "\n", sep = "")
  quit(save = "no", status = 0)
}

# ---- resumo ----
files <- list.files(out_dir, "^tuning_.*_[a-z]{2}[.]csv$", full.names = TRUE)
if (!length(files)) stop("Nenhum tuning_<variante>_<pais>.csv em ", out_dir, call. = FALSE)
all <- bind_rows(lapply(files, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  if (nrow(d)) d$value <- as.character(d$value)
  d
}))
write.csv(all, file.path(out_dir, "tuning_all.csv"), row.names = FALSE)

lam <- all %>% filter(method %in% c("lasso", "ridge"), param == "lambda") %>%
  mutate(lambda = as.numeric(value),
         at_max = abs(lambda - max(lambda_grid)) < 1e-9,
         at_min = abs(lambda - min(lambda_grid)) < 1e-12)
summ <- function(d, ...) d %>% group_by(...) %>%
  summarise(blocks = n(), pct_at_max = round(100 * mean(at_max), 1),
            pct_at_min = round(100 * mean(at_min), 1),
            median_lambda = signif(median(lambda), 3), .groups = "drop")
by_cty  <- summ(lam, variant, method, country)
by_type <- summ(lam, variant, method, country, type)
by_h    <- summ(lam, variant, method, horizon) %>%
  mutate(hn = as.integer(sub("h", "", horizon))) %>% arrange(variant, method, hn) %>% select(-hn)
write.csv(by_type, file.path(out_dir, "lambda_boundary_summary.csv"), row.names = FALSE)

txt <- file.path(out_dir, "lambda_boundary_summary.txt")
sink(txt)
cat("LASSO and Ridge: yearly lambda choices at the limits of the grid (1e-4 to 10)\n")
cat("pct_at_max: share of yearly blocks with lambda = 10 (all coefficients 0 in LASSO:\n")
cat("the forecast is the window's mean change, a drift)\n\n")
cat("=== By variant, method and country ===\n"); print(as.data.frame(by_cty), row.names = FALSE)
cat("\n=== By variant, method and horizon ===\n"); print(as.data.frame(by_h), row.names = FALSE)
cat("\n=== By variant, method, country and specification (DIF/ECM) ===\n"); print(as.data.frame(by_type), row.names = FALSE)
sink()
cat("Resumo gravado em ", txt, "\n", sep = "")
