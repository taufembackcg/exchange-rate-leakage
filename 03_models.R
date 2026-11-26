# TCC forecast: train models (UM PAIS POR PROCESSO)
# Wagner Scheid Mallmann
#
# Substitui, num unico script:
#   - 03_models_v3 - testes.r            (ajuste de todos os metodos)
#   - 03_models_LM_single.R              (LM-ECM pela Eq. 3 de Colombo & Pelagatti)
#   - 03_models_v3-testes_brnn_na_fix.r  (retreino do brnn com NA)
#
# ─────────────────────────────────────────────────────────────────────────────
# CORRECAO: problema 1 de review_problems.md (janela de treino com o futuro)
# ─────────────────────────────────────────────────────────────────────────────
# O painel do 02 guarda, na linha t, os preditores X(t) e a variavel
# dependente y(t+h) = e(t+h) - e(t) (previsao direta). Com
#     trainControl(method = "timeslice", initialWindow = 120, horizon = 1)
# a previsao da linha t era feita com um modelo treinado nas linhas t-120..t-1.
# As respostas das linhas t-h+1..t-1 so sao conhecidas em t+1..t+h-1, depois
# da origem da previsao: para h > 1, o modelo era estimado com h-1 meses de
# cambio futuro.
#
# Agora cada especificacao/horizonte recebe fatias explicitas:
#     treino  = linhas i .. i+119
#     previsao = linha i+119+h
# ou seja, a janela termina em t-h para a previsao da linha t, e todas as
# respostas usadas no treino (datadas t-119..t) ja sao observadas na origem t.
# Para h = 1 nada muda. Cada especificacao perde as h-1 primeiras previsoes.
#
# As fatias sao construidas sobre as linhas COMPLETAS da especificacao
# (resposta + preditores), que e exatamente o conjunto de linhas que o caret
# usa depois do na.omit; rowIndex continua indexando essas linhas, como o
# forecast_alignment.R (04/06) espera.
#
# ─────────────────────────────────────────────────────────────────────────────
# CORRECAO: problema 2 de review_problems.md (tuning escolhido com o futuro)
# ─────────────────────────────────────────────────────────────────────────────
# Antes, o caret avaliava cada valor da grade em TODAS as fatias e o bestTune
# era o valor com o menor erro medio no periodo de avaliacao inteiro: um unico
# parametro escolhido ex post com os proprios erros fora da amostra.
#
# Agora segue Colombo & Pelagatti (2020, Secao 4.2): "the tuning parameters
# were selected using 10-fold cross-validation applied to the first estimate
# of the 10-year rolling window. The parameters selected on the initial window
# were used on all subsequent windows."
#   1) Tuning: validacao cruzada 10-fold nas 120 primeiras linhas completas da
#      especificacao (a primeira janela de estimacao), com a grade original.
#   2) Previsao: as fatias rolantes (problema 1) sao estimadas com o parametro
#      escolhido fixo; cada metodo e ajustado uma vez por janela.
# A grade e o erro de CV da etapa 1 ficam guardados no objeto, em
# $tuning_first_window. O LM nao tem parametro de tuning.
#
# ─────────────────────────────────────────────────────────────────────────────
# REVISAO v8 (reviews_v7.md, M2, M3, M6 e M1; substitui o tuning na 1a janela)
# ─────────────────────────────────────────────────────────────────────────────
# M2 -- re-tuning anual: as previsoes de cada especificacao sao agrupadas em
#   blocos de 12 origens consecutivas. No inicio de cada bloco os parametros de
#   tuning sao escolhidos de novo, so com as 120 linhas da janela de estimacao
#   daquela origem (todas observadas na origem), e ficam fixos nas 12 previsoes
#   do bloco. tune_yearly_roll() faz isso; a grade de cada bloco fica em
#   $tuning_by_year.
# M3 -- validacao cruzada em blocos com intervalo: a janela de 120 linhas e
#   dividida em 10 blocos consecutivos de 12; cada bloco e deixado de fora uma
#   vez e as h-1 linhas de cada lado dele saem do treino (blocked_folds()). As
#   respostas de linhas a menos de h meses de distancia compartilham meses de
#   cambio; dobras aleatorias favoreciam ajustes com sobreajuste.
# M6 -- quebras do M1 (EUA 05/2020, China 01/2024, Indonesia 01/2011): dummies
#   em tempo real, como em Medeiros et al. (2021, JBES), a partir da origem em
#   que o mes da quebra e publicado (M1 com um mes de defasagem). ECM: degrau
#   (1 a partir da publicacao). DIF: 1 nos h meses em que a variacao de h meses
#   atravessa a quebra. Entram so nos modelos com m_hat (nao no de Taylor). Em
#   janelas onde a dummy e constante ela sai (preProcess "zv" do caret).
# M1 -- variantes com violacoes de tempo real (RUN_VARIANT, 01_functions.R):
#   window  janela termina em t-1 (fatias do 03 antigo, problema 1);
#   scaling padronizacao (center/scale) e nos/base das splines na amostra inteira;
#   nolags e filters mudam so o painel do 02.
# O brnn saiu da lista de metodos (ponto 40: o artigo usa doze metodos).
#
# Problema 3 (hiato do produto recursivo) e corrigido no 01_functions.R
# (output_gap_m, usado pelo 02); aqui apenas o modelo de Taylor da Indonesia e
# excluido (ver excluded_models).
#
# Problema 4 (PIB mensal em tempo real) e corrigido no 01/02.
# Problema 19 (nos das B-splines de SP/SPI) e corrigido aqui: ver secao 2,
# spline_glmnet() (spline_glmnet.R), base B-spline construida dentro de cada janela.
#
# ─────────────────────────────────────────────────────────────────────────────
# COMO EXECUTAR
# ─────────────────────────────────────────────────────────────────────────────
#     Rscript --no-save --no-restore "03_models.R" br "D:/Data Exchange Rate/outputs"
# Argumentos posicionais:
#   [1] pais        (obrigatorio)  ex.: br
#   [2] models_dir  (opcional)     default: D:/Data Exchange Rate/outputs
#   [3] n_cores     (opcional)     default: 15
# Saida: <models_dir>/output_models_<pais>.rds (mesma estrutura de antes).

# ─────────────────────────────────────────────
# 0) ARGUMENTOS E AMBIENTE
# ─────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || !nzchar(args[1])) {
  stop("Uso: Rscript 03_models.R <pais> [models_dir] [n_cores]", call. = FALSE)
}
cty        <- tolower(trimws(args[1]))
models_dir <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "D:/Data Exchange Rate/outputs"
n_cores    <- if (length(args) >= 3 && nzchar(args[3])) as.integer(args[3]) else 15L
# [4] metodos a pular (opcional, separados por virgula, ex.: "brnn"). Ponto 37:
# no conjunto de Colombo e Pelagatti (amostra desde 1973, ~500 janelas) o brnn
# com ate 1500 re-estimacoes por especificacao tornava o 03 inviavel.
skip_methods <- if (length(args) >= 4 && nzchar(trimws(args[4]))) trimws(strsplit(args[4], ",")[[1]]) else character(0)

# Progresso (ponto 38): o .bat manda toda a saida do R para o log, entao a janela
# ficava muda. Cada linha de progresso vai para o log e tambem para a janela do
# console (dispositivo "CON" do Windows, que nao passa pelo redirecionamento).
# Fora de um console (ex.: RStudio) a linha vai so para o log.
progress_con <- tryCatch(file("CON", "w"), error = function(e) NULL)
progress <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
  cat(msg, "\n", sep = "")
  if (!is.null(progress_con)) {
    try({ cat(msg, "\n", sep = "", file = progress_con); flush(progress_con) }, silent = TRUE)
  }
}

# Como list(), mas sem avaliar os metodos em skip_methods (nao sao estimados).
# Avalia um metodo por vez, com linhas de progresso (inicio e tempo gasto).
list_skip <- function(...) {
  ex  <- as.list(substitute(list(...)))[-1]
  ex  <- ex[!names(ex) %in% skip_methods]
  env <- parent.frame()
  out <- list()
  for (nm in names(ex)) {
    t1 <- Sys.time()
    progress("    ", nm, " ...")
    out[nm] <- list(eval(ex[[nm]], envir = env))     # list() preserva um NULL
    progress("    ", nm, " ok (", format(round(difftime(Sys.time(), t1), 1)), ")")
  }
  out
}

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

t0 <- Sys.time()
cat("\n===============================================================\n")
cat(" 03_models  |  pais: ", toupper(cty), "\n", sep = "")
cat(" models_dir: ", models_dir, "   n_cores: ", n_cores, "\n", sep = "")
cat(" metodos pulados: ", if (length(skip_methods)) paste(skip_methods, collapse = ", ") else "nenhum", "\n", sep = "")
cat("===============================================================\n\n")

source("01_functions.R")

set.seed(123)
options(timeout = 60 * 60 * 24)
cl <- parallel::makePSOCKcluster(n_cores, timeout = 60 * 60 * 24)
doParallel::registerDoParallel(cl)
on.exit({
  try(parallel::stopCluster(cl), silent = TRUE)
  try(doParallel::stopImplicitCluster(), silent = TRUE)
}, add = TRUE)

ntuneLength <- 15
window_len  <- 120
retune_every <- 12   # M2: re-tuning a cada 12 origens (anual)
cv_folds     <- 10   # M3: blocos da validacao cruzada dentro da janela

variant <- run_variant()
cat(" variante (RUN_VARIANT): ", variant, "\n", sep = "")

# Teste de tempo (timing_test.bat; vazio na rodada real): MAX_FORECASTS limita
# cada especificacao as primeiras N previsoes e ONLY_HORIZONS aos horizontes
# listados (ex.: "h1,h12"). O custo cresce em proporcao ao numero de previsoes
# (janelas de 120 linhas, um tuning a cada 12 previsoes), entao o tempo por
# previsao medido aqui extrapola para a rodada inteira (timing_estimate.R).
# ONLY_HORIZONS tambem e usado na rodada real das variantes de vazamento (M1):
# seus .bat fixam "h1,h3,h6,h12"; o main roda h1 a h12.
max_fc <- suppressWarnings(as.integer(Sys.getenv("MAX_FORECASTS", "0")))
if (is.na(max_fc)) max_fc <- 0L
only_h <- trimws(strsplit(Sys.getenv("ONLY_HORIZONS", ""), ",")[[1]])
if (max_fc > 0 || length(only_h)) cat(" TESTE DE TEMPO: max_forecasts =", max_fc, "| horizontes =", only_h, "\n")

# ─────────────────────────────────────────────
# 1) FATIAS DE PREVISAO DIRETA (correcao do problema 1)
# ─────────────────────────────────────────────
# n: numero de linhas completas da especificacao; h: horizonte.
# createTimeSlices(horizon = h) devolve treino i..i+119 e teste i+120..i+119+h;
# guardamos apenas a ultima linha do teste (i+119+h).
# Variante "window" (violacao): treino i..i+119 e previsao da linha i+120, ou
# seja, a janela termina em t-1 e contem h-1 respostas ainda nao observadas.
direct_slices <- function(n, h, window = window_len) {
  if (variant == "window") {
    sl <- caret::createTimeSlices(seq_len(n), initialWindow = window, horizon = 1, fixedWindow = TRUE)
    out <- list(index = sl$train, indexOut = sl$test)
  } else {
    sl <- caret::createTimeSlices(seq_len(n), initialWindow = window,
                                  horizon = h, fixedWindow = TRUE)
    out <- list(index    = sl$train,
                indexOut = lapply(sl$test, function(x) x[length(x)]))
  }
  if (max_fc > 0 && length(out$index) > max_fc) {        # so no teste de tempo
    out$index <- out$index[seq_len(max_fc)]; out$indexOut <- out$indexOut[seq_len(max_fc)]
  }
  out
}

# ─────────────────────────────────────────────
# 1b) VALIDACAO CRUZADA EM BLOCOS COM INTERVALO (M3)
# ─────────────────────────────────────────────
# n linhas (a janela), k blocos consecutivos; o bloco b fica de fora e as h-1
# linhas de cada lado dele saem do treino.
blocked_folds <- function(n, h, k = cv_folds) {
  blocks <- split(seq_len(n), cut(seq_len(n), k, labels = FALSE))
  gap    <- h - 1L
  out    <- lapply(blocks, function(b) b)
  trn    <- lapply(blocks, function(b) setdiff(seq_len(n), max(1L, min(b) - gap):min(n, max(b) + gap)))
  names(out) <- names(trn) <- sprintf("Fold%02d", seq_along(blocks))
  list(index = trn, indexOut = out)
}

# ─────────────────────────────────────────────
# 1c) RE-TUNING ANUAL + PREVISAO ROLANTE (M2)
# ─────────────────────────────────────────────
# f: formula; d: painel; h: horizonte.
# grid_args: list(tuneGrid = ...) ou list(tuneLength = ...); NULL = sem tuning (LM).
# ...: argumentos fixos do metodo (preProcess, ntree), usados no tuning e na previsao.
# Para cada bloco de 12 origens: tuning na janela da primeira origem do bloco
# (blocked_folds) e previsao das 12 origens com o parametro escolhido. Cada
# chamada do caret recebe so as linhas do bloco (janelas + linhas previstas);
# rowIndex e devolvido na numeracao das linhas completas da especificacao,
# como o forecast_alignment.R (04/06) espera.
tune_yearly_roll <- function(method, f, d, h, grid_args = NULL, ...) {
  dcc <- d[stats::complete.cases(d[, all.vars(f)]), ]
  sl  <- direct_slices(nrow(dcc), h)
  roll_ctrl <- function(idx, out) trainControl(method = "cv", index = idx, indexOut = out,
                                               allowParallel = TRUE, savePredictions = TRUE,
                                               returnResamp = "none")
  groups <- if (is.null(grid_args)) list(seq_along(sl$index))
            else split(seq_along(sl$index), ceiling(seq_along(sl$index) / retune_every))
  preds <- list(); tunes <- list(); last <- NULL
  for (g in seq_along(groups)) {
    ii   <- groups[[g]]
    rows <- min(unlist(sl$index[ii])):max(unlist(sl$indexOut[ii]))
    off  <- min(rows) - 1L
    sub  <- dcc[rows, , drop = FALSE]
    idx  <- lapply(sl$index[ii],    function(x) x - off)
    out  <- lapply(sl$indexOut[ii], function(x) x - off)
    tg   <- NULL
    if (!is.null(grid_args)) {
      tw <- sl$index[[ii[1]]]                      # janela da primeira origem do bloco
      bf <- blocked_folds(length(tw), h)
      set.seed(123)
      tuned <- suppressWarnings(do.call(train, c(
        list(f, data = dcc[tw, , drop = FALSE], method = method, metric = "RMSE",
             trControl = trainControl(method = "cv", index = bf$index, indexOut = bf$indexOut,
                                      allowParallel = TRUE)),
        grid_args, list(...))))
      tg <- tuned$bestTune
      tunes[[g]] <- data.frame(block = g, first_rowIndex = sl$indexOut[[ii[1]]], tg, check.names = FALSE)
    }
    fit_args <- c(list(f, data = sub, method = method, metric = "RMSE", trControl = roll_ctrl(idx, out)),
                  if (!is.null(tg)) list(tuneGrid = tg), list(...))
    last <- suppressWarnings(do.call(train, fit_args))
    p <- last$pred[, c("rowIndex", "pred", "obs")]
    p$rowIndex <- p$rowIndex + off
    preds[[g]] <- p
  }
  last$pred <- do.call(rbind, preds)
  last$bestTune <- NULL            # tuning varia por bloco: o 04/06 le $pred inteiro
  last$tuning_by_year <- if (length(tunes)) do.call(rbind, tunes) else NULL
  last
}

# ─────────────────────────────────────────────
# 2) B-SPLINES COM NOS POR JANELA: spline_glmnet() em spline_glmnet.R
#    (problema 19; truncamento nos nos de fronteira, ponto 34)
# ─────────────────────────────────────────────
source("spline_glmnet.R")

# ─────────────────────────────────────────────
# 2b) GAMBOOST COM DUMMIES LINEARES (M6)
# ─────────────────────────────────────────────
# O gamboost do caret usa base B-spline (bbs) para todo preditor; numa dummy de
# quebra (so 0 e 1, quase sempre constante) o ajuste fica singular. Aqui cada
# dummy entra com base linear (bols) e os demais preditores com bbs, como no
# gamboost padrao (mesmos graus de liberdade). Sem dummies, o modelo e o mesmo.
gamboost_lin <- local({
  m <- caret::getModelInfo("gamboost", regex = FALSE)$gamboost
  m$library <- c("mboost", "plyr", "import")
  m$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
    # gamboost procura bbs/bols pelo nome: o mboost precisa estar anexado (o caret
    # original faz import::from(mboost, bbs)); vale tambem nos processos paralelos.
    if (!"package:mboost" %in% search()) suppressPackageStartupMessages(library(mboost))
    theDots <- list(...)
    ctl <- mboost::boost_control(mstop = param$mstop)
    theDots$family <- mboost::GaussReg()
    if (!is.null(wts)) theDots$weights <- wts
    dat <- as.data.frame(x, stringsAsFactors = TRUE)
    is_dummy <- vapply(dat, function(v) all(v %in% c(0, 1)), logical(1))
    terms <- ifelse(is_dummy, paste0("bols(`", names(dat), "`)"), paste0("bbs(`", names(dat), "`)"))
    dat$.outcome <- y
    f <- stats::as.formula(paste(".outcome ~", paste(terms, collapse = " + ")))
    out <- do.call(mboost::gamboost, c(list(formula = f, data = dat, control = ctl), theDots))
    if (param$prune == "yes") {
      iters <- mboost::mstop(AIC(out))
      if (iters < out$mstop()) out <- out[iters]
    }
    out$.org.mstop <- out$mstop()
    out
  }
  m
})

# ─────────────────────────────────────────────
# 3) DADOS E GRADE
# ─────────────────────────────────────────────
data <- readr::read_csv2(
  panel_path(cty),
  col_types = readr::cols(date = readr::col_date(format = ""), .default = readr::col_double())
) |> dplyr::arrange(date)

vars <- list(
  taylor                = c("ygap_hat", "r_hat", "i_hat"),
  monetary_standard     = c("m_hat", "y_hat"),
  monetary_sticky       = c("m_hat", "y_hat", "r_hat"),
  monetary_sticky_uip   = c("m_hat", "y_hat", "r_hat", "i_hat"),
  monetary_flexible_uip = c("m_hat", "y_hat", "i_hat")
)
wdw           <- "rolling"
modeltypelist <- c("dif", "ecm")

# Problema 3 (hiato recursivo): com o hiato de Hamilton estimado em tempo real
# (min. 36 obs.), o modelo de Taylor da Indonesia nao tem previsoes a partir de
# h ~ 6 (PIB mensal comeca em 2008). Decisao: excluir o modelo de Taylor para a
# Indonesia em todos os horizontes; 04/06 ignoram chaves ausentes.
# Ponto 36 (conjunto de Colombo e Pelagatti): mesmo motivo para a Franca. O PIB
# comeca em 1980 e o franco termina em 12/1998; o hiato recursivo so existe a
# partir de ~1986, o que deixaria o modelo de Taylor com ~20-30 previsoes e
# encolheria a janela comum de todos os modelos.
excluded_models <- list(id = c("taylor"), fr = c("taylor"))
min_forecasts   <- 12   # salvaguarda: especificacoes com menos previsoes sao puladas
horizonlist   <- paste0("h", 1:12)
if (length(only_h)) horizonlist <- intersect(horizonlist, only_h)   # teste de tempo e variantes M1

# Ponto 33 de review_problems.md: margem epsilon fixa do SVR linear (svmlw2)
svr_eps_svmlw2 <- 0.001

en_tune_grid <- expand.grid(alpha = seq(0, 1, length.out = 11),
                            lambda = 10^seq(-4, 1, length.out = 50))

# ─────────────────────────────────────────────
# 3b) DUMMIES DAS QUEBRAS DO M1 (M6)
# ─────────────────────────────────────────────
# Mes da quebra (referencia). A quebra dos EUA entra no m_hat de todos os paises.
# A dummy comeca no mes em que a quebra e publicada: mes + defasagem do M1 (1;
# 0 na variante nolags, em que o M1 entra no mes de referencia), a mesma data
# em que a quebra aparece no painel do 02.
m1_breaks   <- list(us = "2020-05-01", cn = "2024-01-01", id = "2011-01-01")
m1_brk_lag  <- if (variant == "nolags") 0L else 1L
brk_names   <- intersect(names(m1_breaks), c("us", cty))
brk_ecm     <- character(0)
for (b in brk_names) {
  rel <- lubridate::`%m+%`(as.Date(m1_breaks[[b]]), months(m1_brk_lag))
  nm  <- paste0(cty, "_dbrk_", b)
  data[[nm]] <- as.numeric(data$date >= rel)                                   # ECM: degrau
  brk_ecm <- c(brk_ecm, nm)
  for (h in 1:12) {                                                            # DIF: h meses
    data[[paste0(nm, "_dif_h", h)]] <- as.numeric(data$date >= rel &
      data$date <= lubridate::`%m+%`(rel, months(h - 1)))
  }
}
brk_vars <- function(mdl_type, hzn) {
  if (mdl_type == "ecm") brk_ecm else if (length(brk_ecm)) paste0(brk_ecm, "_dif_", hzn) else character(0)
}

# ─────────────────────────────────────────────
# 3c) VARIANTE "scaling" (M1): padronizacao na amostra inteira
# ─────────────────────────────────────────────
# Violacao: cada preditor (exceto as dummies) e padronizado com a media e o
# desvio de toda a amostra antes do 03, e o center/scale do caret (por janela)
# sai. Os nos e a base das splines (SP/SPI) tambem vem da amostra inteira.
# (glmnet e kernlab ainda reescalam internamente com os dados de cada janela.)
pp <- function(...) c("zv", ...)                   # "zv": tira colunas constantes na janela (dummies)
if (variant == "scaling") {
  pred_cols <- grep(paste0("^", cty, "_(m_hat|y_hat|r_hat|i_hat|ygap_hat)"), names(data), value = TRUE)
  pred_cols <- c(pred_cols, paste0(cty, "_et"))
  for (v in pred_cols) {
    s <- stats::sd(data[[v]], na.rm = TRUE)
    if (is.finite(s) && s > 0) data[[v]] <- (data[[v]] - mean(data[[v]], na.rm = TRUE)) / s
  }
  pp <- function(...) "zv"
}

outputmodels <- list()

# ─────────────────────────────────────────────
# 4) LOOP: especificacao x horizonte x modelo
# ─────────────────────────────────────────────
n_total <- length(modeltypelist) * length(horizonlist) * length(vars)
i_spec  <- 0L
progress(toupper(cty), ": ", n_total, " especificacoes; variante ", variant, "; metodos pulados: ",
         if (length(skip_methods)) paste(skip_methods, collapse = ", ") else "nenhum")

for (mdl_type in modeltypelist) {
  for (hzn in horizonlist) {
    h_num <- as.integer(sub("h", "", hzn))

    for (mdl in names(vars)) {
      i_spec <- i_spec + 1L
      response_var <- paste0(cty, "_et_", mdl_type, "_", hzn)
      pred_vars <- if (mdl_type == "ecm") {
        c(paste0(cty, "_", vars[[mdl]]), paste0(cty, "_et"))            # X(t) e e(t)
      } else {
        paste0(cty, "_", vars[[mdl]], "_", mdl_type, "_", hzn)          # X(t) - X(t-h)
      }
      missing_cols <- setdiff(c(response_var, pred_vars), names(data))
      if (length(missing_cols) > 0) {
        warning("Colunas ausentes (", paste(missing_cols, collapse = ", "), "): ",
                mdl_type, " ", hzn, " ", mdl, " -- ignorado.")
        next
      }
      # M6: dummies de quebra nos modelos com m_hat (sem NA: nao mudam as linhas completas)
      if ("m_hat" %in% vars[[mdl]]) pred_vars <- c(pred_vars, brk_vars(mdl_type, hzn))

      key <- paste0(cty, "_", wdw, "_", mdl_type, "_", hzn, "_form_", mdl)
      if (mdl %in% excluded_models[[cty]]) {
        progress(key, " (", i_spec, "/", n_total, ") -- excluido (decisao do problema 3)")
        next
      }
      n_cc <- sum(stats::complete.cases(data[, c(response_var, pred_vars)]))
      n_fc <- length(direct_slices(max(n_cc, window_len + h_num), h_num)$index)
      if (n_cc < window_len + h_num || n_fc < min_forecasts) {
        n_fc <- if (n_cc < window_len + h_num) 0 else n_fc
        warning(key, ": apenas ", n_fc, " previsoes possiveis -- ignorado.")
        progress(key, " (", i_spec, "/", n_total, ") -- ignorado: ", n_fc, " previsoes")
        next
      }
      progress(key, " (", i_spec, "/", n_total, ")")

      # LM-ECM pela Eq. (3) de Colombo & Pelagatti: MQO irrestrito em (X(t), e(t)),
      # o mesmo conjunto de regressores dos demais metodos.
      train_formula <- reformulate(pred_vars, response = response_var)

      # Re-tuning anual (M2) com blocos (M3) + previsao rolante
      fit <- function(method, grid_args = NULL, ..., f = train_formula, d = data) {
        tune_yearly_roll(method, f, d, h_num, grid_args, ...)
      }
      grid  <- function(g) list(tuneGrid = g)
      tlen  <- list(tuneLength = ntuneLength)
      # Variante scaling: nos e base das splines na amostra inteira da especificacao
      sp_global <- if (variant == "scaling") {
        data[stats::complete.cases(data[, c(response_var, pred_vars)]), pred_vars, drop = FALSE]
      } else NULL

      outputmodels[[key]] <- list_skip(
        lm       = fit("lm", preProcess = "zv"),
        lasso    = fit("glmnet", grid(expand.grid(alpha = 1, lambda = 10^seq(-4, 1, length.out = 50))),
                       preProcess = pp("center", "scale")),
        ridge    = fit("glmnet", grid(expand.grid(alpha = 0, lambda = 10^seq(-4, 1, length.out = 50))),
                       preProcess = pp("center", "scale")),
        svmr     = fit("svmRadial", tlen, preProcess = pp("center", "scale")),
        # Ponto 33: margem epsilon do SVR = 0.001 (padrao do LiblineaR, 0.1, maior que as
        # variacoes do cambio da China: ajuste plano, previsao zero = random walk)
        svmlw2   = fit("svmLinear3", grid(expand.grid(cost = 10^seq(-3, 2, length.out = 10), Loss = c("L1", "L2"))),
                       preProcess = pp("center", "scale"), svr_eps = svr_eps_svmlw2),
        forest   = fit("rf", tlen, ntree = 1000, preProcess = "zv"),
        gam      = fit("gamSpline", tlen, preProcess = pp("center", "scale")),
        sp       = fit(spline_glmnet(interactions = FALSE, global_x = sp_global), grid(en_tune_grid),
                       preProcess = "zv"),                                    # base por janela
        spi      = fit(spline_glmnet(interactions = TRUE,  global_x = sp_global), grid(en_tune_grid),
                       preProcess = "zv"),                                    # base por janela
        earth    = fit("earth", grid(expand.grid(nprune = c(2, 5, 10, 15, 20), degree = c(1, 2))),
                       preProcess = "zv"),
        gamboost = fit(gamboost_lin, grid(expand.grid(mstop = c(50, 100, 150, 200), prune = c(TRUE, FALSE))),
                       preProcess = "zv"),
        pcr      = fit("pcr", tlen, preProcess = pp("center", "scale"))
      )
    }
  }
}

# ─────────────────────────────────────────────
# 5) SALVA
# ─────────────────────────────────────────────
dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(models_dir, paste0("output_models_", cty, ".rds"))
progress("salvando ", out_path)
saveRDS(outputmodels, out_path)
progress("fim")

tf <- Sys.time()
cat("\n---------------------------------------------------------------\n")
cat(" pais   : ", toupper(cty), "\n", sep = "")
cat(" chaves : ", length(outputmodels), "\n", sep = "")
cat(" arquivo: ", out_path, "\n", sep = "")
cat(" tempo  : ", format(round(difftime(tf, t0), 2)), "\n", sep = "")
cat("---------------------------------------------------------------\n")
cat("OK\n")
