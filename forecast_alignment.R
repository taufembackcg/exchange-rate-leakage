# ___________________________________________________________________________________
# __________________ forecast_alignment.R (usado por 04 e 06) _______________________
# ___________________________________________________________________________________
#
# Funcoes compartilhadas por 04_analysis_single.R e 06_mcs_single.R para extrair
# as previsoes fora da amostra dos objetos caret e alinha-las corretamente no
# tempo. Mantidas num unico arquivo para que 04 (RMSE relativa, DM) e 06 (MCS)
# usem exatamente as mesmas previsoes e as mesmas datas.
#
# ─────────────────────────────────────────────────────────────────────────────
# PROBLEMA CORRIGIDO (problema 0 de review_problems.md)
# ─────────────────────────────────────────────────────────────────────────────
# O 03 chama train() com formula e na.action = na.omit. O caret remove as
# linhas incompletas ANTES de criar as fatias, entao model$pred$rowIndex indexa
# as linhas do painel FILTRADO. Os scripts antigos faziam
#     preds$date <- data$date[preds$rowIndex]
# no painel NAO filtrado, que comeca em 1996 para todos os paises. Sempre que
# havia linhas iniciais incompletas, as previsoes recebiam datas anteriores as
# verdadeiras (ate 191 meses para a Indonesia). Consequencias: o RMSE do random
# walk era calculado sobre outro periodo (desde 2006), o teste DM pareava erros
# de meses diferentes e o MCS alinhava perdas de meses diferentes.
#
# CORRECAO
#   1) As linhas usadas pelo caret sao reconstruidas com complete.cases() sobre
#      as mesmas variaveis do 03 (resposta + preditores daquela especificacao).
#      A data de cada previsao e a data dessa linha filtrada.
#   2) Verificacao obrigatoria: o valor observado guardado pelo caret
#      (preds$obs) tem de ser igual a resposta do painel na data atribuida.
#      Se nao for, o script para (nao produz numeros desalinhados).
#   3) Random walk sem drift: a previsao da variacao e zero, logo o erro do RW
#      na data t e -(e_{t+h} - e_t), que e exatamente -obs. O RW e portanto
#      calculado nas MESMAS datas das previsoes dos modelos.
#   4) Janela comum (opcao B): para cada pais e horizonte, todas as
#      especificacoes (5 modelos x 2 especificacoes x 13 metodos) sao avaliadas
#      nos mesmos meses, a intersecao das datas com previsao valida.

suppressPackageStartupMessages({
  library(dplyr)
})

# Variaveis exatamente como o 03 constroi cada especificacao ----------------------
fa_vars <- list(
  taylor                = c("ygap_hat", "r_hat", "i_hat"),
  monetary_standard     = c("m_hat", "y_hat"),
  monetary_sticky       = c("m_hat", "y_hat", "r_hat"),
  monetary_sticky_uip   = c("m_hat", "y_hat", "r_hat", "i_hat"),
  monetary_flexible_uip = c("m_hat", "y_hat", "i_hat")
)

fa_spec_vars <- function(cty, mdl_type, mdl, hzn) {
  response <- paste0(cty, "_et_", mdl_type, "_", hzn)
  if (mdl_type == "ecm") {
    preds <- c(paste0(cty, "_", fa_vars[[mdl]]), paste0(cty, "_et"))
  } else {
    preds <- paste0(cty, "_", fa_vars[[mdl]], "_", mdl_type, "_", hzn)
  }
  list(response = response, preds = preds)
}

# Previsoes do bestTune (mesma regra do 04 original, incluindo o fallback
# quando o bestTune tem NA). Com o 03_models.R corrigido (problema 2), o
# tuning e fixado na primeira janela e $pred ja contem um unico valor de
# tuning; o filtro abaixo fica trivial, mas continua valido para os .rds antigos.
# Revisao v8 (re-tuning anual, M2 de reviews_v7.md): o 03 remove $bestTune e
# $pred traz uma previsao por linha, cada uma com o tuning do seu bloco; sem
# bestTune nenhum filtro e aplicado.
fa_best_preds <- function(model_obj) {
  preds <- model_obj$pred
  if (is.null(preds) || nrow(preds) == 0) return(NULL)

  if ("bestTune" %in% names(model_obj)) {
    best <- model_obj$bestTune
    for (col in names(best)) preds <- preds[preds[[col]] == best[[col]], ]
  }
  if ("bestTune" %in% names(model_obj) && any(is.na(preds$pred))) {
    preds      <- model_obj$pred
    col        <- names(model_obj$bestTune)
    candidates <- unique(preds[[col]])
    RMSEvector <- rep(Inf, length(candidates))
    for (cand in seq_along(candidates)) {
      preds_tmp <- preds[preds[[col]] == candidates[cand], ]
      if (any(is.na(preds_tmp$pred))) next
      RMSEvector[cand] <- sqrt(mean((preds_tmp$obs - as.numeric(preds_tmp$pred))^2))
    }
    best_idx <- which.min(RMSEvector)
    preds    <- preds[preds[[col]] == candidates[best_idx], ]
  }
  if (!is.numeric(preds$pred)) preds$pred <- as.numeric(preds$pred)
  preds[order(preds$rowIndex), c("rowIndex", "pred", "obs")]
}

# Atribui a data correta a cada previsao e verifica o alinhamento --------------
# Tenta (1) as variaveis da especificacao do 03; (2) as variaveis do proprio
# objeto (model_obj$terms), que cobrem o LM-ECM antigo com termo cointegrante
# caso algum pais nao tenha sido corrigido pelo 03_models_LM_single.R.
fa_assign_dates <- function(preds, data, model_obj, sv, tol = 1e-8) {
  cand_sets <- list(c(sv$response, sv$preds))
  if (!is.null(model_obj$terms)) {
    tv <- all.vars(model_obj$terms)
    if (length(tv) > 0 && all(tv %in% names(data))) cand_sets[[2]] <- tv
  }
  for (vars in cand_sets) {
    if (!all(vars %in% names(data))) next
    ok    <- stats::complete.cases(data[, vars])
    dates <- data$date[ok]
    y     <- data[[sv$response]][ok]
    if (max(preds$rowIndex) > length(dates)) next
    if (isTRUE(all(abs(y[preds$rowIndex] - preds$obs) < tol))) {
      preds$date <- dates[preds$rowIndex]
      return(preds)
    }
  }
  NULL
}

# Coleta todas as previsoes de um pais/horizonte -------------------------------
# Retorna um data.frame longo: spec, model_type, model, method, date, pred, obs.
# skip_methods (ponto 40): metodos ignorados na avaliacao (ex.: "brnn"), de modo
# que a janela comum e o Top-10 do MCS sao formados sem eles.
fa_collect <- function(models, data, cty, hzn, model_list, model_type_list,
                       wdw = "rolling", skip_methods = character(0)) {
  out <- list()
  for (mdl_type in model_type_list) {
    for (mdl in model_list) {
      key <- paste0(cty, "_", wdw, "_", mdl_type, "_", hzn, "_form_", mdl)
      model_group <- models[[key]]
      if (is.null(model_group)) next
      sv <- fa_spec_vars(cty, mdl_type, mdl, hzn)

      for (method_name in setdiff(names(model_group), skip_methods)) {
        model_obj <- model_group[[method_name]]
        if (is.null(model_obj)) next
        preds <- fa_best_preds(model_obj)
        if (is.null(preds)) next

        preds <- fa_assign_dates(preds, data, model_obj, sv)
        if (is.null(preds)) {
          stop(sprintf(paste0("Desalinhamento: %s %s %s %s [%s] -- preds$obs nao bate com ",
                              "a resposta do painel em nenhuma definicao de linhas. ",
                              "Nenhum resultado foi gravado."),
                       cty, mdl_type, hzn, mdl, method_name), call. = FALSE)
        }
        preds <- preds[!is.na(preds$pred), ]
        if (nrow(preds) == 0) next

        out[[length(out) + 1]] <- data.frame(
          spec       = paste(mdl_type, mdl, method_name, sep = "|"),
          model_type = mdl_type,
          model      = mdl,
          method     = method_name,
          date       = as.Date(preds$date),
          pred       = preds$pred,
          obs        = preds$obs,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(out) == 0) return(NULL)
  bind_rows(out)
}

# Janela comum (opcao B) -------------------------------------------------------
# Mantem apenas os meses em que TODAS as especificacoes tem previsao valida e
# verifica que o valor realizado e o mesmo em todas elas (DIF e ECM tem a mesma
# variavel dependente e_{t+h} - e_t) e igual ao do painel.
fa_common_window <- function(fc, data, cty, hzn, tol = 1e-8) {
  n_specs <- length(unique(fc$spec))
  common  <- fc %>%
    count(date) %>%
    filter(n == n_specs) %>%
    pull(date)

  fc <- fc %>% filter(date %in% common) %>% arrange(spec, date)

  y_panel <- data[[paste0(cty, "_et_dif_", hzn)]][match(fc$date, data$date)]
  if (any(is.na(y_panel)) || max(abs(y_panel - fc$obs)) > tol) {
    stop(sprintf("Valores realizados inconsistentes entre especificacoes: %s %s", cty, hzn),
         call. = FALSE)
  }
  fc
}

# Diagnostico da janela: n por especificacao antes da intersecao ---------------
fa_window_report <- function(fc_all, fc_common, cty, hzn) {
  per_spec <- fc_all %>%
    group_by(spec) %>%
    summarise(n_own = n(), first_own = min(date), last_own = max(date), .groups = "drop")
  data.frame(
    country        = cty,
    horizon        = hzn,
    n_specs        = nrow(per_spec),
    common_start   = min(fc_common$date),
    common_end     = max(fc_common$date),
    n_common       = length(unique(fc_common$date)),
    latest_start   = max(per_spec$first_own),
    binding_spec   = per_spec$spec[which.max(per_spec$first_own)],
    min_n_own      = min(per_spec$n_own),
    stringsAsFactors = FALSE
  )
}
