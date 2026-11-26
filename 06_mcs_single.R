# ___________________________________________________________________________________
# _____________________ 06_mcs_single.R (UM PAIS POR PROCESSO) ______________________
# ___________________________________________________________________________________
#
# Substitui 06_mcs_analysis_top10.R. Corrige dois bugs confirmados na extracao
# dos resultados do MCS (verificados com um exemplo sintetico do pacote MCS,
# NAO com os dados do projeto -- ver Codes/inspect_mcs_object.R):
#
# BUG 1 -- colunas em branco (Rank / MCS p-value)
#   O script original tentava casar os nomes de coluna de mcs_out@show contra
#   uma lista de candidatos adivinhados (ex.: "p.Value.for.H_0.M_k."). O nome
#   real, apos make.names(), e "p.Value.for.H_.0.M_k." (ha um ponto extra antes
#   do "0"). Nenhum candidato batia, entao as colunas saiam vazias.
#   CORRIGIDO: mcs_out@show e sempre uma matriz 3 colunas, SEMPRE na ordem
#   Avg.Loss | p-Value for H_{0,M_k} | MCS p-Value (confirmado no pacote MCS).
#   Agora as colunas sao lidas por POSICAO, nao por nome, o que elimina essa
#   classe de bug por completo.
#
# BUG 2 -- "Random Walk" listado nao significa que ele sobreviveu
#   mcs_out@show contem TODOS os modelos originalmente testados (sobreviventes
#   E eliminados), em ordem de eliminacao -- exatamente a tabela impressa por
#   print(mcs_out). O script original tratava a mera PRESENCA do nome "Random
#   Walk" nessa tabela como prova de sobrevivencia, o que esta errado: um
#   modelo eliminado no primeiro passo ainda aparece como a primeira linha.
#   CORRIGIDO: a sobrevivencia real vem de mcs_out@Info$included (nomes dos
#   modelos no conjunto final) vs mcs_out@Info$excluded (eliminados).
#
# Tambem corrigido: ordenacao de horizonte como fator numerico (h1..h12) em
# vez de texto (h1, h10, h11, h12, h2, h3, ...).
#
# BUG 3 -- datas das previsoes (problema 0 de review_problems.md)
#   preds$date <- data_cty$date[preds$rowIndex] usava o painel completo,
#   mas rowIndex indexa o painel filtrado pelo caret (na.omit). As perdas de
#   modelos diferentes e do random walk eram alinhadas em meses errados.
#   CORRIGIDO: previsoes extraidas por forecast_alignment.R (o mesmo codigo do
#   04): datas reconstruidas e conferidas contra o painel, janela comum a todas
#   as especificacoes do pais naquele horizonte (opcao B), perda do random walk
#   = obs^2 nas mesmas datas. O Top-10 e escolhido pela perda media NESSA
#   janela comum, e nao mais na janela propria de cada modelo.
#
# ─────────────────────────────────────────────────────────────────────────────
# COMO EXECUTAR
# ─────────────────────────────────────────────────────────────────────────────
# Um pais por processo (mesmo motivo do 03/04: os .rds sao grandes demais para
# carregar em sequencia numa unica sessao R). Use o driver run_06.ps1, ou
# manualmente:
#
#     Rscript --no-save --no-restore "06_mcs_single.R" br "D:/Data Exchange Rate/outputs"
#
# Argumentos posicionais:
#   [1] pais        (obrigatorio)  ex.: br
#   [2] models_dir  (opcional)     default: D:/Data Exchange Rate/outputs
#
# ─────────────────────────────────────────────────────────────────────────────
# SAIDA
# ─────────────────────────────────────────────────────────────────────────────
# outputs/tables/table_mcs_full_<pais>.tex
#   Tabela longa (Horizon | Model | Avg. Loss | p-value (elim.) | MCS p-value),
#   uma linha por modelo testado em cada horizonte, na ordem de eliminacao
#   (pior/eliminado primeiro), exatamente as tres colunas numericas impressas
#   por print(mcs_out). Destinada ao Material Suplementar (Secao S2).
#
# outputs/tables/mcs_rw_pvalue_<pais>.csv
#   Uma linha por horizonte: o MCS p-value do Random Walk e se ele sobreviveu
#   (TRUE/FALSE, via Info$included). Usado para montar a tabela do artigo
#   principal (Secao 4.4) com o p-value real do RW por pais/horizonte.

# ─────────────────────────────────────────────
# 0) ARGUMENTOS E AMBIENTE
# ─────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1 || !nzchar(args[1])) {
  stop("Uso: Rscript 06_mcs_single.R <pais> [models_dir]", call. = FALSE)
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
cat(" 06_mcs_single  |  pais: ", toupper(cty), "\n", sep = "")
cat(" models_dir    : ", models_dir, "\n", sep = "")
cat(" wd            : ", getwd(), "\n", sep = "")
cat("===============================================================\n\n")

# ─────────────────────────────────────────────
# 1) PACOTES E FUNCOES AUXILIARES
# ─────────────────────────────────────────────
library(MCS)
library(dplyr)
library(purrr)
library(readr)
library(kableExtra)
library(stringr)

source("01_functions.R")
source("forecast_alignment.R")
# Variante da rodada (M1 de reviews_v7.md): painel da variante e resultados em
# outputs/ (main) ou outputs_<variante>/.
OUT <- variant_out_dir()

dir.create(file.path(OUT, "tables"), showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────────────
# 2) PARAMETROS (mesma configuracao do 06_mcs_analysis_top10.R original)
# ─────────────────────────────────────────────
mcs_alpha      <- 0.10
mcs_B          <- 5000
mcs_statistic  <- "Tmax"
mcs_k_min      <- 3   # ponto 27: comprimento do bloco = max(3, h)
# CORRECAO (ponto 27 de review_problems.md): as perdas de previsoes h passos a
# frente sao autocorrelacionadas ate a defasagem h - 1 (MA(h-1)). Com bloco
# fixo de 3, o bootstrap subestimava essa dependencia para h >= 4 e os
# p-valores do MCS ficavam baixos demais (o RW era eliminado com facilidade).
mcs_block_len  <- function(h) max(mcs_k_min, h)
# Semente do bootstrap (2026-11-25): sem semente, os p-valores do MCS mudavam
# de uma rodada para outra e as combinacoes perto de 0,10 podiam trocar de
# lado (Japao, h3, window: 0,096 e depois 0,104). A semente depende so do pais
# e do horizonte, entao o resultado nao depende da ordem dos jobs, e os cinco
# desenhos usam os mesmos numeros aleatorios para o mesmo pais e horizonte.
mcs_seed_base  <- 20261125L
mcs_seed       <- function(cty, h) {
  mcs_seed_base + 100L * sum(utf8ToInt(tolower(cty)) * c(1L, 7L)) + as.integer(h)
}
n_top_models   <- 10
min_common_obs <- 10

model_list      <- c("taylor", "monetary_standard", "monetary_sticky",
                     "monetary_sticky_uip", "monetary_flexible_uip")
model_list      <- eval_models(cty, model_list)   # M4: India sem o modelo de Taylor (01_functions.R)
model_type_list <- c("dif", "ecm")
horizon_list    <- paste0("h", 1:12)

# Nomes dos modelos (problema 13): nomes de Colombo & Pelagatti (2020, Secao
# 4.1), na ordem da Tabela 1; o modelo de precos flexiveis com juros
# (Frenkel-Bilson, ausente em Colombo) recebe nome no mesmo estilo.
model_labels <- c(
  monetary_standard      = "Monetary model",
  monetary_sticky        = "Monetary model, sticky prices",
  monetary_sticky_uip    = "Monetary model, sticky prices, UIP deviations",
  monetary_flexible_uip  = "Monetary model, UIP deviations",
  taylor                 = "Taylor rule model"
)

# ─────────────────────────────────────────────
# 3) CAMINHOS E VALIDACAO
# ─────────────────────────────────────────────
models_path <- file.path(models_dir, paste0("output_models_", cty, ".rds"))
data_path   <- panel_path(cty)

if (!file.exists(models_path)) stop("Modelos nao encontrados: ", models_path, call. = FALSE)
if (!file.exists(data_path))   stop("Panel data nao encontrado: ", data_path, call. = FALSE)

cat("Lendo ", models_path, " (",
    round(file.info(models_path)$size / 1024^2, 1), " MB) ...\n", sep = "")
models_cty <- readRDS(models_path)

data_cty <- read_csv2(
  data_path,
  col_types = cols(date = col_date(format = ""), .default = col_double())
) %>% arrange(date)

# ─────────────────────────────────────────────
# 4) HELPERS (identicos ao script original -- nao alterados)
# ─────────────────────────────────────────────
sanitize_loss_matrix <- function(loss_mat) {
  if (is.null(loss_mat) || ncol(loss_mat) < 2) return(NULL)

  row_ok <- apply(loss_mat, 1, function(z) all(is.finite(z)) && !any(is.na(z)))
  loss_mat <- loss_mat[row_ok, , drop = FALSE]
  if (nrow(loss_mat) < min_common_obs) return(NULL)

  keep <- sapply(seq_len(ncol(loss_mat)), function(j) {
    x <- loss_mat[, j]
    length(x) >= min_common_obs && is.finite(sd(x)) && sd(x) > 0
  })
  loss_mat <- loss_mat[, keep, drop = FALSE]
  if (ncol(loss_mat) < 2) return(NULL)

  if (ncol(loss_mat) >= 2) {
    dup <- duplicated(as.data.frame(t(loss_mat)), fromLast = FALSE)
    loss_mat <- loss_mat[, !dup, drop = FALSE]
  }

  if (ncol(loss_mat) < 2) return(NULL)
  loss_mat
}

# Matriz de perdas Top-N + random walk na janela comum (opcao B).
# Retorna list(loss_mat, window) ou NULL.
build_topN_loss_matrix <- function(models_cty, data_cty, cty, hzn, wdw = "rolling") {
  fc_all <- fa_collect(models_cty, data_cty, cty, hzn, model_list, model_type_list, wdw, skip_methods = skip_methods)
  if (is.null(fc_all)) return(NULL)
  fc  <- fa_common_window(fc_all, data_cty, cty, hzn)
  win <- fa_window_report(fc_all, fc, cty, hzn)
  if (win$n_common < min_common_obs) return(NULL)

  common_dates <- sort(unique(fc$date))

  # Perda media de cada especificacao na janela comum -> Top-N
  mean_loss <- fc %>%
    group_by(spec) %>%
    summarise(loss = mean((pred - obs)^2), .groups = "drop") %>%
    arrange(loss)
  top_specs <- head(mean_loss$spec, n_top_models)

  # Random walk: previsao zero, perda = obs^2 nas mesmas datas
  rw_obs <- fc %>% distinct(date, obs) %>% arrange(date)
  loss_cols <- list()
  loss_cols[["Random Walk"]] <- rw_obs$obs[match(common_dates, rw_obs$date)]^2

  for (sp in top_specs) {
    d_sp  <- fc[fc$spec == sp, ]
    ord   <- match(common_dates, d_sp$date)
    label <- paste0(model_labels[d_sp$model[1]], " [", d_sp$method[1], "] (", d_sp$model_type[1], ")")
    loss_cols[[label]] <- (d_sp$pred[ord] - d_sp$obs[ord])^2
  }

  loss_mat <- do.call(cbind, loss_cols)
  colnames(loss_mat) <- names(loss_cols)
  list(loss_mat = sanitize_loss_matrix(loss_mat), window = win)
}

# ─────────────────────────────────────────────
# 5) EXTRACAO CORRIGIDA
# ─────────────────────────────────────────────
# mcs_out@show: matriz [modelos x 3], SEMPRE nesta ordem de colunas
#   (Avg.Loss, p-Value for H_{0,M_k}, MCS p-Value), linhas em ordem de
#   eliminacao (pior/eliminado primeiro). Lido por posicao, nao por nome.
# mcs_out@Info$included / $excluded: nomes dos modelos no conjunto final
#   (sobreviventes) e dos eliminados, respectivamente -- fonte de verdade
#   para sobrevivencia, independente do texto impresso.
extract_mcs_full <- function(mcs_out, cty, hzn) {
  show_mat <- mcs_out@show
  if (is.null(show_mat) || nrow(show_mat) == 0) return(NULL)

  included <- mcs_out@Info$included

  data.frame(
    country      = cty,
    horizon      = hzn,
    model        = rownames(show_mat),
    avg_loss     = as.numeric(show_mat[, 1]),
    p_value_elim = as.numeric(show_mat[, 2]),
    mcs_p_value  = as.numeric(show_mat[, 3]),
    survived     = rownames(show_mat) %in% included,
    stringsAsFactors = FALSE
  )
}

# ─────────────────────────────────────────────
# 6) LOOP DE HORIZONTES
# ─────────────────────────────────────────────
all_results <- list()
all_windows <- list()

for (hzn in horizon_list) {

  cat("Rodando MCS: ", toupper(cty), " / ", hzn, " ... ", sep = "")
  mcs_k <- mcs_block_len(as.integer(sub("h", "", hzn)))

  built <- build_topN_loss_matrix(models_cty, data_cty, cty, hzn)
  loss_mat <- if (is.null(built)) NULL else built$loss_mat
  if (is.null(loss_mat) || ncol(loss_mat) < 2) {
    cat("dados insuficientes, ignorado.\n")
    next
  }
  all_windows[[hzn]] <- built$window
  cat("janela ", format(built$window$common_start), " a ", format(built$window$common_end),
      " (n = ", nrow(loss_mat), ") ... ", sep = "")

  set.seed(mcs_seed(cty, as.integer(sub("h", "", hzn))))
  mcs_out <- tryCatch(
    suppressWarnings(
      MCSprocedure(
        Loss = loss_mat,
        alpha = mcs_alpha,
        B = mcs_B,
        statistic = mcs_statistic,
        k = mcs_k
      )
    ),
    error = function(e) {
      cat("ERRO: ", conditionMessage(e), "\n", sep = "")
      NULL
    }
  )
  if (is.null(mcs_out)) next

  res <- extract_mcs_full(mcs_out, cty, hzn)
  if (is.null(res)) {
    cat("sem resultados extraiveis.\n")
    next
  }

  all_results[[hzn]] <- res
  n_elim <- sum(!res$survived)
  cat(nrow(res), " modelos, ", n_elim, " eliminados.\n", sep = "")
}

if (length(all_results) == 0) {
  stop("Nenhum resultado de MCS produzido para ", cty, call. = FALSE)
}

results_df <- bind_rows(all_results)

# ─────────────────────────────────────────────
# 7) TABELA COMPLETA (Material Suplementar) -- ordem numerica de horizonte,
#    ordem de eliminacao preservada dentro de cada horizonte
# ─────────────────────────────────────────────
results_df <- results_df %>%
  mutate(horizon_f = factor(horizon, levels = horizon_list)) %>%
  arrange(horizon_f)

tab <- results_df %>%
  transmute(
    Horizon = horizon,
    Model = paste0("\\raggedright ", model),   # quebra entre palavras, sem hifenizar
    `Avg. Loss` = sprintf("%.4f", avg_loss),
    `p-value (elim.)` = sprintf("%.3f", p_value_elim),
    `MCS p-value` = sprintf("%.3f", mcs_p_value)
  )

# CORRIGIDO: o linesep padrao do kableExtra insere \addlinespace a cada 5
# linhas, sem relacao com o agrupamento real dos dados (cada horizonte tem
# ~11 linhas, uma por modelo candidato). Aqui o espaco extra e inserido
# apenas quando o horizonte muda de fato.
n_tab <- nrow(tab)
linesep_vec <- ifelse(as.character(tab$Horizon[-n_tab]) != as.character(tab$Horizon[-1]),
                       "\\addlinespace", "")

tex <- kbl(
  tab, format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE,
  caption = paste0("Model Confidence Set -- full elimination sequence, Top-",
                   n_top_models, " pre-filter -- ", toupper(cty)),
  label = paste0("mcs_full_", cty),
  # Problema 13: nomes de Colombo (mais longos) -> coluna Model com largura
  # fixa, p{5.5cm}, quebrando em duas linhas, e espacamento entre colunas de
  # 3pt (abaixo). p{} simples nao exige o pacote array.
  align = c("c", "p{5.5cm}", "c", "c", "c"),
  linesep = linesep_vec
) %>%
  kable_styling(latex_options = c("repeat_header"), font_size = 8) %>%
  footnote(
    general = paste0(
      "For each horizon, models are listed in elimination order (worst/first-eliminated at top). ",
      "Avg. Loss is the mean squared forecast loss; p-value (elim.) is the p-value at the step the ",
      "model was tested for elimination; MCS p-value is the running maximum, i.e. the model's final ",
      "MCS p-value. A model is retained in the final Model Confidence Set when its MCS p-value is at ",
      "least alpha = ", mcs_alpha, ". Statistic = ", mcs_statistic, ", block length k = max(", mcs_k_min,
      ", h), B = ", mcs_B, "."
    ),
    general_title = "Note:", footnote_as_chunk = TRUE, threeparttable = TRUE
  )

tex <- sub("\\selectfont", "\\selectfont\\setlength{\\tabcolsep}{3pt}", as.character(tex), fixed = TRUE)
fname_tex <- paste0(OUT, "/tables/table_mcs_full_", cty, ".tex")
writeLines(as.character(tex), fname_tex)

# ─────────────────────────────────────────────
# 8) RESUMO: p-value do Random Walk por horizonte (artigo principal, Secao 4.4)
# ─────────────────────────────────────────────
windows_df <- bind_rows(all_windows) %>%
  transmute(horizon, eval_start = common_start, eval_end = common_end, n_obs = n_common)

rw_summary <- results_df %>%
  filter(model == "Random Walk") %>%
  left_join(windows_df, by = "horizon") %>%
  transmute(
    country = country,
    horizon = factor(horizon, levels = horizon_list),
    rw_mcs_p_value = mcs_p_value,
    rw_survived = survived,
    eval_start, eval_end, n_obs
  ) %>%
  arrange(horizon)

fname_csv <- paste0(OUT, "/tables/mcs_rw_pvalue_", cty, ".csv")
write_csv(rw_summary, fname_csv)

# ─────────────────────────────────────────────
# 9) LOG FINAL
# ─────────────────────────────────────────────
tf <- Sys.time()

cat("\n---------------------------------------------------------------\n")
cat(" pais              : ", toupper(cty), "\n", sep = "")
cat(" horizontes         : ", length(all_results), "/", length(horizon_list), "\n", sep = "")
cat(" tabela completa    : ", fname_tex, "\n", sep = "")
cat(" p-values do RW     : ", fname_csv, "\n", sep = "")
print(rw_summary)
cat(" tempo              : ", format(round(difftime(tf, t0), 2)), "\n", sep = "")
cat("---------------------------------------------------------------\n")
cat("OK\n")
