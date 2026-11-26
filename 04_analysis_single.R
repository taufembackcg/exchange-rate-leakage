# TCC forecast: analysis (UM PAIS POR PROCESSO)
#
# ─────────────────────────────────────────────────────────────────────────────
# VERSAO CORRIGIDA -- problema 0 de review_problems.md (datas das previsoes)
# ─────────────────────────────────────────────────────────────────────────────
# Correcoes desta versao (detalhes em forecast_alignment.R):
#   1) Datas: model$pred$rowIndex indexa o painel FILTRADO pelo caret
#      (na.omit). A versao anterior usava data$date[rowIndex] no painel
#      completo e atribuia datas erradas (ate 191 meses antes). Agora as linhas
#      do caret sao reconstruidas e cada previsao e conferida contra a resposta
#      do painel; se nao bater, o script para.
#   2) Random walk: calculado nas MESMAS datas de cada previsao (erro do RW =
#      -obs). A versao anterior usava rmse_rw(), que media o RW sempre desde
#      a linha 121 do painel (2006-01), qualquer que fosse a janela do modelo.
#   3) Janela comum (opcao B): em cada horizonte, todas as especificacoes do
#      pais sao avaliadas nos mesmos meses (intersecao). RMSE, RMSE relativa e
#      DM passam a ser comparaveis entre modelos, metodos e especificacoes.
#
# Mantido da versao anterior: indexacao do RW por horizonte, horizontes h1..h12,
# escolha do bestTune pelo caret (problema 2, ainda nao tratado aqui).
#
# ─────────────────────────────────────────────────────────────────────────────
# COMO EXECUTAR
# ─────────────────────────────────────────────────────────────────────────────
# Um pais por processo, via run_04.ps1, ou manualmente:
#
#     Rscript --no-save --no-restore "04_analysis_single.R" br "D:/Data Exchange Rate/outputs"
#
# Argumentos posicionais:
#   [1] pais        (obrigatorio)  ex.: br
#   [2] models_dir  (opcional)     default: D:/Data Exchange Rate/outputs
#
# Saida: outputs/Resultados_Todos_<pais>.xlsx
#   - Resultados_Todos / Resultados_<PAIS>: mesmas colunas de antes (05 nao muda)
#     + n_obs, eval_start, eval_end, cw_test_p_value (Clark-West, ponto 30)
#   - Janela_Avaliacao: janela comum por horizonte e a especificacao que a limita
#   - RW_RMSE: RMSE do random walk na janela comum, por horizonte

# ─────────────────────────────────────────────
# 0) ARGUMENTOS E AMBIENTE
# ─────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1 || !nzchar(args[1])) {
  stop("Uso: Rscript 04_analysis_single.R <pais> [models_dir]", call. = FALSE)
}

cty        <- tolower(trimws(args[1]))
models_dir <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "D:/Data Exchange Rate/outputs"
# [3] metodos excluidos da avaliacao (opcional, separados por virgula, ex.: "brnn").
# Ponto 40: o artigo (v6) reporta os mesmos doze metodos nas duas amostras.
skip_methods <- if (length(args) >= 3 && nzchar(trimws(args[3]))) trimws(strsplit(args[3], ",")[[1]]) else character(0)

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

if (is.null(getOption("repos")) || getOption("repos")[["CRAN"]] %in% c("", "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

t0 <- Sys.time()

cat("\n===============================================================\n")
cat(" 04_analysis_single  |  pais: ", toupper(cty), "\n", sep = "")
cat(" models_dir         : ", models_dir, "\n", sep = "")
cat(" wd                 : ", getwd(), "\n", sep = "")
cat("===============================================================\n\n")

# ─────────────────────────────────────────────
# 1) PACOTES E FUNCOES AUXILIARES
# ─────────────────────────────────────────────
source("01_functions.R")
source("forecast_alignment.R")
# Variante da rodada (M1 de reviews_v7.md): painel da variante e resultados em
# outputs/ (main) ou outputs_<variante>/.
OUT <- variant_out_dir()

RMSE_manual <- function(pred, obs) sqrt(mean((pred - obs)^2))
min_common_obs <- 24   # aviso se a janela comum ficar menor que isso

# Teste de Clark e West (2007) -- ponto 30 de review_problems.md.
# O random walk (previsao zero) e um caso aninhado dos modelos (coeficientes
# nulos); sob H0 o DM com valores criticos normais fica subdimensionado. CW
# corrige o EQM do modelo pelo ruido de estimar parametros que sao zero sob H0:
#     f_t = e_rw^2 - [ e_mod^2 - (yhat_rw - yhat_mod)^2 ],  yhat_rw = 0
#         = obs^2 - (pred - obs)^2 + pred^2
# H0: E[f_t] = 0 contra H1: E[f_t] > 0 (o modelo supera o RW), estatistica t
# com variancia de longo prazo Newey-West (Bartlett, h - 1 defasagens, as
# mesmas do DM) e valor critico normal unilateral.
cw_test <- function(pred, obs, h) {
  f <- obs^2 - (pred - obs)^2 + pred^2
  n <- length(f)
  if (n < 3 || !all(is.finite(f))) return(NA_real_)
  u  <- f - mean(f)
  lr <- sum(u^2) / n
  L  <- min(h - 1, n - 1)
  if (L > 0) for (j in seq_len(L)) {
    lr <- lr + 2 * (1 - j / (L + 1)) * sum(u[(j + 1):n] * u[1:(n - j)]) / n
  }
  if (!is.finite(lr) || lr <= 0) return(NA_real_)
  stats::pnorm(mean(f) / sqrt(lr / n), lower.tail = FALSE)
}

# ─────────────────────────────────────────────
# 2) CAMINHOS E VALIDACAO
# ─────────────────────────────────────────────
models_path <- file.path(models_dir, paste0("output_models_", cty, ".rds"))
data_path   <- panel_path(cty)

if (!file.exists(models_path)) stop("Modelos nao encontrados: ", models_path, call. = FALSE)
if (!file.exists(data_path))   stop("Panel data nao encontrado: ", data_path, call. = FALSE)

cat("Tamanho do .rds: ", round(file.info(models_path)$size / 1024^2, 1), " MB\n", sep = "")
cat("Lendo ", models_path, " ...\n", sep = "")
models <- readRDS(models_path)

data <- read_csv2(
  data_path,
  col_types = cols(date = col_date(format = ""), .default = col_double())
) %>%
  arrange(date)

# ─────────────────────────────────────────────
# 3) GRADE
# ─────────────────────────────────────────────
model_list      <- c("taylor", "monetary_standard", "monetary_sticky", "monetary_sticky_uip", "monetary_flexible_uip")
model_list      <- eval_models(cty, model_list)   # M4: India sem o modelo de Taylor (01_functions.R)
wdw             <- "rolling"
model_type_list <- c("dif", "ecm")
horizon_list    <- paste0("h", 1:12)

results_list <- list()
window_list  <- list()
rw_list      <- list()

# ─────────────────────────────────────────────
# 4) LOOP POR HORIZONTE
# ─────────────────────────────────────────────
for (hzn in horizon_list) {
  h_num <- as.numeric(gsub("h", "", hzn))

  fc_all <- fa_collect(models, data, cty, hzn, model_list, model_type_list, wdw, skip_methods = skip_methods)
  if (is.null(fc_all)) {
    cat(sprintf("  %s: nenhuma previsao encontrada, ignorado.\n", hzn))
    next
  }
  fc <- fa_common_window(fc_all, data, cty, hzn)
  win <- fa_window_report(fc_all, fc, cty, hzn)
  window_list[[hzn]] <- win

  cat(sprintf("  %-3s | %3d especificacoes | janela comum %s a %s (n = %d) | limitada por %s\n",
              hzn, win$n_specs, win$common_start, win$common_end, win$n_common, win$binding_spec))
  if (win$n_common < min_common_obs) {
    cat(sprintf("  AVISO: janela comum com apenas %d meses em %s.\n", win$n_common, hzn))
  }

  # Random walk na janela comum: previsao de variacao zero, erro = -obs.
  rw_obs  <- fc %>% distinct(date, obs) %>% arrange(date)
  rmse_rw <- sqrt(mean(rw_obs$obs^2))
  rw_list[[hzn]] <- data.frame(country = cty, horizon = hzn, rmse_rw = rmse_rw,
                               n_obs = nrow(rw_obs))

  for (sp in unique(fc$spec)) {
    d_sp <- fc[fc$spec == sp, ]
    d_sp <- d_sp[order(d_sp$date), ]

    e_model <- d_sp$pred - d_sp$obs
    e_rw    <- -d_sp$obs                       # mesmas datas, mesma ordem
    rmse_val <- RMSE_manual(d_sp$pred, d_sp$obs)

    dm_pval <- tryCatch(
      forecast::dm.test(
        e1 = e_rw,
        e2 = e_model,
        h  = h_num,
        power = 2,
        alternative = "greater",
        varestimator = "bartlett"
      )[["p.value"]],
      error = function(e) NA_real_
    )
    cw_pval <- cw_test(d_sp$pred, d_sp$obs, h_num)

    results_list[[length(results_list) + 1]] <- data.frame(
      country         = cty,
      model_type      = d_sp$model_type[1],
      horizon         = hzn,
      window          = wdw,
      model           = d_sp$model[1],
      method          = d_sp$method[1],
      rmse            = rmse_val,
      relative_rmse   = rmse_val / rmse_rw,
      dm_test_p_value = dm_pval,
      cw_test_p_value = cw_pval,
      n_obs           = nrow(d_sp),
      eval_start      = min(d_sp$date),
      eval_end        = max(d_sp$date),
      stringsAsFactors = FALSE
    )
  }
}

if (length(results_list) == 0) {
  stop("Nenhum resultado calculado para ", cty,
       ". O .rds contem as chaves esperadas?", call. = FALSE)
}

# ─────────────────────────────────────────────
# 5) COMBINA E EXPORTA
# ─────────────────────────────────────────────
results_all <- bind_rows(results_list)
window_all  <- bind_rows(window_list)
rw_all      <- bind_rows(rw_list)

sheets <- list(results_all, results_all, window_all, rw_all)
names(sheets) <- c("Resultados_Todos", paste0("Resultados_", toupper(cty)),
                   "Janela_Avaliacao", "RW_RMSE")

if (!dir.exists(OUT)) dir.create(OUT)
namefile <- paste0(OUT, "/Resultados_Todos_", cty, ".xlsx")
writexl::write_xlsx(sheets, path = namefile)

rm(models)
gc()

tf <- Sys.time()

cat("\n---------------------------------------------------------------\n")
cat(" pais       : ", toupper(cty), "\n", sep = "")
cat(" linhas     : ", nrow(results_all), "\n", sep = "")
cat(" arquivo    : ", namefile, "\n", sep = "")
cat(" RW RMSE na janela comum:\n")
print(rw_all[, c("horizon", "rmse_rw", "n_obs")], row.names = FALSE)
cat(" tempo      : ", format(round(difftime(tf, t0), 2)), "\n", sep = "")
cat("---------------------------------------------------------------\n")
cat("OK\n")
