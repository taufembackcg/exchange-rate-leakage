# READ-ONLY aggregation of already-computed results for the paper's text.
# Does not retrain or re-estimate anything; only reads the outputs of 04 and 06:
#   outputs/Resultados_Todos_<cty>.xlsx       (04: RMSE, relative RMSE, DM on the
#                                              common evaluation window)
#   outputs/tables/mcs_rw_pvalue_<cty>.csv    (06: random walk MCS p-value and
#                                              survival, per horizon)
#
# CORRECTION (point 21 of review_problems.md):
#   - MCS results now come from mcs_rw_pvalue_<cty>.csv. The previous version
#     parsed the old table_mcs_top10_<cty>.tex and treated "Random Walk appears
#     in the table" as survival, but that table lists every tested model,
#     eliminated or not. Survival now comes from 06 (MCS Info$included).
#   - Runs from the script's own folder (no hard-coded setwd).
#   - All previous aggregates are kept unchanged. Added: first horizon at which
#     the random walk is eliminated from the MCS, common evaluation windows
#     (sheet Janela_Avaliacao of 04), and model names (Colombo and Pelagatti
#     labels, point 13).

get_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", ca, value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  }
  getwd()
}
setwd(get_script_dir())
# Variante da rodada (M1 de reviews_v7.md; RUN_VARIANT no .bat): resultados em
# outputs/ (main) ou outputs_<variante>/.
OUT <- local({ v <- tolower(trimws(Sys.getenv("RUN_VARIANT", "main")))
               if (!nzchar(v) || v == "main") "outputs" else paste0("outputs_", v) })

library(readxl)
library(dplyr)
library(jsonlite)
library(stringr)

# Paises: os que tem outputs/Resultados_Todos_<pais>.xlsx (ponto 36: serve ao
# conjunto BRICS e ao de Colombo e Pelagatti, cada um na sua pasta de execucao).
countries <- sub("^Resultados_Todos_(.*)[.]xlsx$", "\\1",
                 list.files(OUT, pattern = "^Resultados_Todos_.*[.]xlsx$"))
if (length(countries) == 0) countries <- c("br", "cn", "id", "in", "ru", "za")
horizons_order <- paste0("h", 1:12)

model_labels <- c(
  monetary_standard      = "Monetary model",
  monetary_sticky        = "Monetary model, sticky prices",
  monetary_sticky_uip    = "Monetary model, sticky prices, UIP deviations",
  monetary_flexible_uip  = "Monetary model, UIP deviations",
  taylor                 = "Taylor rule model"
)

results_all <- purrr::map_dfr(countries, function(cty) {
  f <- paste0(OUT, "/Resultados_Todos_", cty, ".xlsx")
  if (!file.exists(f)) { warning("Missing: ", f); return(NULL) }
  read_excel(f, sheet = "Resultados_Todos")
})

cat("Total rows:", nrow(results_all), "\n")
cat("Countries:", paste(unique(results_all$country), collapse=", "), "\n")
cat("Horizons:", paste(sort(unique(results_all$horizon)), collapse=", "), "\n")
cat("Model types:", paste(unique(results_all$model_type), collapse=", "), "\n")
cat("Methods:", paste(sort(unique(results_all$method)), collapse=", "), "\n")
cat("Models:", paste(unique(results_all$model), collapse=", "), "\n\n")

# Ponto 30: Clark-West ao lado do DM (NA se o 04 for anterior a essa versao)
if (!"cw_test_p_value" %in% names(results_all)) results_all$cw_test_p_value <- NA_real_

results_all <- results_all %>%
  mutate(
    beats_rw = relative_rmse < 1,
    sig05 = dm_test_p_value < 0.05,
    sig10 = dm_test_p_value < 0.10,
    cw_sig05 = cw_test_p_value < 0.05,
    cw_sig10 = cw_test_p_value < 0.10
  )

n_total <- nrow(results_all)
n_beat  <- sum(results_all$beats_rw, na.rm = TRUE)
n_sig05_among_beat <- sum(results_all$beats_rw & results_all$sig05, na.rm = TRUE)

out <- list()
out$n_total <- n_total
out$n_beat <- n_beat
out$pct_beat <- round(100*n_beat/n_total,1)
out$n_sig05_among_beat <- n_sig05_among_beat
out$pct_sig05_among_beat <- round(100*n_sig05_among_beat/n_beat,1)
out$n_cw_sig05 <- sum(results_all$cw_sig05, na.rm = TRUE)
out$pct_cw_sig05 <- round(100*out$n_cw_sig05/n_total,1)
out$n_cw_sig05_among_beat <- sum(results_all$beats_rw & results_all$cw_sig05, na.rm = TRUE)
out$pct_cw_sig05_among_beat <- round(100*out$n_cw_sig05_among_beat/n_beat,1)

# by horizon
by_h <- results_all %>%
  group_by(horizon) %>%
  summarise(n=n(), n_beat=sum(beats_rw,na.rm=TRUE), pct_beat=round(100*n_beat/n,1),
            n_sig05=sum(sig05,na.rm=TRUE), pct_sig05=round(100*n_sig05/n,1),
            n_cw_sig05=sum(cw_sig05,na.rm=TRUE), pct_cw_sig05=round(100*n_cw_sig05/n,1),.groups="drop") %>%
  mutate(horizon=factor(horizon, levels=horizons_order)) %>%
  arrange(horizon)
out$by_horizon <- by_h

# h1 vs h12
out$h1 <- as.list(by_h %>% filter(horizon=="h1"))
out$h12 <- as.list(by_h %>% filter(horizon=="h12"))

# share of wins occurring at h12
out$share_wins_at_h12 <- round(100*out$h12$n_beat/n_beat,1)

# DIF vs ECM
by_spec <- results_all %>%
  group_by(model_type) %>%
  summarise(n=n(), n_beat=sum(beats_rw,na.rm=TRUE), pct_beat=round(100*n_beat/n,1),
            n_sig05_among_beat=sum(beats_rw&sig05,na.rm=TRUE),
            pct_sig05_among_beat=round(100*n_sig05_among_beat/n_beat,1),.groups="drop")
out$by_spec <- by_spec

# by method
by_method <- results_all %>%
  group_by(method) %>%
  summarise(n=n(), n_beat=sum(beats_rw,na.rm=TRUE), pct_beat=round(100*n_beat/n,1),
            n_sig05=sum(sig05,na.rm=TRUE), pct_sig05=round(100*n_sig05/n,1),
            pct_sig05_among_beat=round(100*sum(beats_rw&sig05,na.rm=TRUE)/n_beat,1),
            .groups="drop") %>%
  arrange(desc(pct_beat))
out$by_method <- by_method

# linear vs nonlinear
linear_methods <- c("lm","lasso","ridge","pcr")
by_class <- results_all %>%
  mutate(class = ifelse(method %in% linear_methods, "linear", "nonlinear")) %>%
  group_by(class) %>%
  summarise(n=n(), n_beat=sum(beats_rw,na.rm=TRUE), pct_beat=round(100*n_beat/n,1),.groups="drop")
out$by_class <- by_class

# by country
by_country <- results_all %>%
  group_by(country) %>%
  summarise(n=n(), n_beat=sum(beats_rw,na.rm=TRUE), pct_beat=round(100*n_beat/n,1),
            n_sig05=sum(sig05,na.rm=TRUE), pct_sig05=round(100*n_sig05/n,1),.groups="drop") %>%
  arrange(desc(pct_beat))
out$by_country <- by_country

# by economic model (with Colombo and Pelagatti labels)
by_model <- results_all %>%
  group_by(model) %>%
  summarise(n=n(), n_beat=sum(beats_rw,na.rm=TRUE), pct_beat=round(100*n_beat/n,1),
            n_sig05=sum(sig05,na.rm=TRUE), pct_sig05=round(100*n_sig05/n,1),.groups="drop") %>%
  mutate(model_label = unname(model_labels[model])) %>%
  relocate(model_label, .after = model) %>%
  arrange(desc(pct_beat))
out$by_model <- by_model

# lm specifically in ecm vs dif (post-patch check)
lm_check <- results_all %>%
  filter(method == "lm") %>%
  group_by(model_type) %>%
  summarise(n=n(), n_beat=sum(beats_rw,na.rm=TRUE), pct_beat=round(100*n_beat/n,1),.groups="drop")
out$lm_check <- lm_check

# best method overall (lowest mean relative rmse)
best_method_rmse <- results_all %>%
  group_by(method) %>%
  summarise(mean_rel_rmse = round(mean(relative_rmse, na.rm=TRUE),3), .groups="drop") %>%
  arrange(mean_rel_rmse)
out$best_method_rmse <- best_method_rmse

# ---- Common evaluation windows (04, sheet Janela_Avaliacao) ----
windows <- purrr::map_dfr(countries, function(cty) {
  f <- paste0(OUT, "/Resultados_Todos_", cty, ".xlsx")
  if (!file.exists(f) || !"Janela_Avaliacao" %in% excel_sheets(f)) return(NULL)
  read_excel(f, sheet = "Janela_Avaliacao")
})
if (nrow(windows) > 0) {
  windows <- windows %>%
    mutate(horizon = factor(horizon, levels = horizons_order),
           common_start = as.Date(common_start), common_end = as.Date(common_end)) %>%
    arrange(country, horizon)
}
out$eval_windows <- windows

# ---- MCS: random walk p-value and survival from 06 (mcs_rw_pvalue_<cty>.csv) ----
mcs_rw <- purrr::map_dfr(countries, function(cty) {
  f <- paste0(OUT, "/tables/mcs_rw_pvalue_", cty, ".csv")
  if (!file.exists(f)) { warning("Missing: ", f); return(NULL) }
  read.csv(f, stringsAsFactors = FALSE)
})
mcs_rw <- mcs_rw %>%
  transmute(country, horizon = factor(horizon, levels = horizons_order),
            rw_mcs_p_value, rw_in_mcs = as.logical(rw_survived),
            eval_start, eval_end, n_obs) %>%
  arrange(country, horizon)
out$mcs_rw <- mcs_rw

mcs_summary <- mcs_rw %>%
  group_by(country) %>%
  summarise(n_horizons = n(), n_rw_survives = sum(rw_in_mcs), pct_rw_survives = round(100*n_rw_survives/n_horizons,1),
            .groups="drop")
out$mcs_summary <- mcs_summary

overall_mcs <- mcs_rw %>%
  summarise(n=n(), n_rw=sum(rw_in_mcs), pct_rw=round(100*n_rw/n,1))
out$overall_mcs <- overall_mcs

# first horizon at which the random walk is eliminated, and last horizon of an
# uninterrupted run of survival from h1 (NA if eliminated at h1)
mcs_threshold <- mcs_rw %>%
  group_by(country) %>%
  summarise(
    first_elim_horizon = if (any(!rw_in_mcs)) as.character(horizon[which(!rw_in_mcs)[1]]) else NA_character_,
    retained_through   = { k <- which(!rw_in_mcs)[1]
                           # horizon e fator: o ultimo horizonte e o ultimo da linha (ordenado acima)
                           if (is.na(k)) as.character(horizon[length(horizon)]) else if (k == 1) NA_character_ else as.character(horizon[k - 1]) },
    survives_after_first_elim = { k <- which(!rw_in_mcs)[1]
                                  if (is.na(k)) FALSE else any(rw_in_mcs[seq_along(rw_in_mcs) > k]) },
    .groups = "drop")
out$mcs_threshold <- mcs_threshold

# ---- write everything to a readable txt + json ----
sink(file.path(OUT, "aggregate_stats_for_paper.txt"))
cat("=== OVERALL ===\n")
print(out[c("n_total","n_beat","pct_beat","n_sig05_among_beat","pct_sig05_among_beat",
            "n_cw_sig05","pct_cw_sig05","n_cw_sig05_among_beat","pct_cw_sig05_among_beat")])
cat("\n=== BY HORIZON ===\n")
print(as.data.frame(by_h))
cat("\n=== BY SPECIFICATION (DIF vs ECM) ===\n")
print(as.data.frame(by_spec))
cat("\n=== BY METHOD ===\n")
print(as.data.frame(by_method))
cat("\n=== LINEAR vs NONLINEAR ===\n")
print(as.data.frame(by_class))
cat("\n=== BY COUNTRY ===\n")
print(as.data.frame(by_country))
cat("\n=== BY ECONOMIC MODEL ===\n")
print(as.data.frame(by_model))
cat("\n=== LM METHOD ONLY (post ECM patch check) ===\n")
print(as.data.frame(lm_check))
cat("\n=== MEAN RELATIVE RMSE BY METHOD (overall) ===\n")
print(as.data.frame(best_method_rmse))
cat("\n=== COMMON EVALUATION WINDOWS BY COUNTRY x HORIZON ===\n")
print(as.data.frame(windows))
cat("\n=== MCS: RANDOM WALK P-VALUE AND SURVIVAL BY COUNTRY x HORIZON (top10 pre-filter) ===\n")
print(as.data.frame(mcs_rw))
cat("\n=== MCS: RW SURVIVAL SUMMARY BY COUNTRY ===\n")
print(as.data.frame(mcs_summary))
cat("\n=== MCS: RW SURVIVAL OVERALL ===\n")
print(as.data.frame(overall_mcs))
cat("\n=== MCS: FIRST HORIZON AT WHICH THE RW IS ELIMINATED ===\n")
print(as.data.frame(mcs_threshold))
sink()

jsonlite::write_json(out, file.path(OUT, "aggregate_stats_for_paper.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")

cat("\nDone. Wrote outputs/aggregate_stats_for_paper.txt and .json\n")
