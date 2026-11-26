# TCC forecast: helper functions
# Wagner Scheid Mallmann
# 2025/02


# _____________________________________________________________________________
# __________________________________ PACKAGES __________________________________
# _____________________________________________________________________________

{
  pacotes_necessarios <- list(
    "lubridate", "tidyverse", "httr", "tidyr", "stringr", "zoo", "writexl",
    "jsonlite", "tsibble", "fable", "feasts", "sidrar", "forecast",
    "ipeadatar", "neverhpfilter", "caret", "doParallel", "ggplot2",
    "fredr", "tempdisagg", "rsdmx", "dplyr", "purrr", "ecb", "rlang", "urca", "elasticnet", "glmnet",
    "LiblineaR","earth","mboost","import","brnn","pls"
  )
  
  for (i in pacotes_necessarios) {
    eval(bquote(
      if (!require(.(i))) {
        install.packages(.(i))
        library(.(i))
      }
    ))
  }
  
  rm(pacotes_necessarios, i)
}

# _______________________________________________________________________________
# __________________________________ FUNCTIONS __________________________________
# _______________________________________________________________________________

#FUNÇÃO 001 : Hiato do produto via Hamilton Filter -- versao RECURSIVA (tempo real)
#
# CORRECAO: problema 3 de review_problems.md. A versao anterior aplicava
# neverhpfilter::yth_filter() uma unica vez na amostra inteira: os coeficientes
# da regressao de Hamilton usavam dados posteriores a cada data de previsao.
#
# Agora, para cada mes t, a regressao de Hamilton (2018)
#     y_s = b0 + b1 y_{s-24} + ... + b12 y_{s-35} + v_s
# e estimada por MQO apenas com s <= t (janela expansiva). A tendencia em t e o
# valor ajustado em t e o hiato e gap_t = y_t / tendencia_t (mesma definicao,
# razao, da versao anterior). Cada valor do hiato usa somente informacao ate t.
# A primeira estimativa exige min_obs observacoes da regressao (default 36,
# cerca de 3x o numero de coeficientes); antes disso o hiato e NA.
# Com min_obs = Inf nao ha recursao: o resultado coincide com yth_filter().
output_gap_m <- function(data, colname = "br_gap", h = 24, p = 12, min_obs = 36) {
  data <- data[order(data[[1]]), ]
  y    <- data[[2]]
  n    <- length(y)
  X    <- sapply(0:(p - 1), function(j) dplyr::lag(y, h + j))
  ok   <- which(stats::complete.cases(cbind(y, X)))
  gap  <- rep(NA_real_, n)

  fit_trend <- function(rows, t) {
    b <- qr.coef(qr(cbind(1, X[rows, , drop = FALSE])), y[rows])
    b[is.na(b)] <- 0
    sum(c(1, X[t, ]) * b)
  }

  if (is.infinite(min_obs)) {
    for (t in ok) gap[t] <- y[t] / fit_trend(ok, t)
  } else {
    for (k in seq_along(ok)) {
      if (k < min_obs) next
      t <- ok[k]
      gap[t] <- y[t] / fit_trend(ok[seq_len(k)], t)
    }
  }

  out <- data.frame(date = lubridate::ymd(data[[1]]), gap = gap)
  colnames(out) <- c("date", colname)
  out
}


# FUNÇÃO 002 : Erro Quadrático Médio (RMSE) do Random Walk com base no painel
# rmse_rw <- function(cty, start_window = 121) {
#   
#   data <- read_csv2(
#     paste0("outputs/", cty, "_panel_data.csv"),
#     col_types = cols(date = col_date(format = ""), .default = col_double())
#   ) %>%
#     arrange(date)
#   
#   # Usa as diferenças já calculadas no painel
#   diff_h1  <- data[[paste0(cty, "_et_dif_h1")]]
#   diff_h12 <- data[[paste0(cty, "_et_dif_h12")]]
#   
#   # RW alinhado ao período fora da amostra (a partir de start_window)
#   df_rw <- data %>%
#     select(date) %>%
#     mutate(
#       error_rw_h1  = diff_h1,
#       error_rw_h12 = diff_h12
#     ) %>%
#     slice(start_window:n()) %>%   # desloca início do RW para 121º mês
#     drop_na()
#   
#   # Calcula RMSE
#   rmse_rw_h1  <- sqrt(mean(df_rw$error_rw_h1^2,  na.rm = TRUE))
#   rmse_rw_h12 <- sqrt(mean(df_rw$error_rw_h12^2, na.rm = TRUE))
#   
#   return(list(
#     rmse_rw      = c(rmse_rw_h1, rmse_rw_h12),
#     data         = df_rw$date,
#     error_rw_h1  = df_rw$error_rw_h1,
#     error_rw_h12 = df_rw$error_rw_h12
#   ))
# }

rmse_rw <- function(cty, start_window = 121, horizons = 1:12) {
  
  data <- read_csv2(
    paste0("outputs/", cty, "_panel_data.csv"),
    col_types = cols(date = col_date(format = ""), .default = col_double())
  ) %>%
    arrange(date)
  
  # Nomes das colunas de erro para cada horizonte
  error_cols <- paste0("error_rw_h", horizons)
  diff_cols  <- paste0(cty, "_et_dif_h", horizons)
  
  # Monta o data frame com todas as diferenças já calculadas no painel
  df_rw <- data %>%
    select(date, all_of(diff_cols)) %>%
    rename_with(~ error_cols, all_of(diff_cols)) %>%
    slice(start_window:n()) %>%   # desloca início do RW para 121º mês
    drop_na()
  
  # Calcula RMSE para cada horizonte
  rmse_rw_vec <- sapply(error_cols, function(col) {
    sqrt(mean(df_rw[[col]]^2, na.rm = TRUE))
  })
  names(rmse_rw_vec) <- paste0("h", horizons)
  
  # Monta a lista de retorno, incluindo os erros de cada horizonte
  result <- list(
    rmse_rw = rmse_rw_vec,
    data    = df_rw$date
  )
  
  for (h in horizons) {
    result[[paste0("error_rw_h", h)]] <- df_rw[[paste0("error_rw_h", h)]]
  }
  
  return(result)
}



#FUNÇÃO 003 : Coleta de dados do FMI 
fetch_imf_series <- function(flowref, pattern_suffix, countries, startPeriod = "2000") {
  get_one <- function(cc) {
    key <- if (identical(cc, "G163") && identical(flowref, "IMF.STA,CPI")) {
      paste(cc, "HICP._T.IX.M", sep = ".")} else {paste(cc, pattern_suffix, sep = ".")}
    
    obj <- readSDMX(
      providerId = "IMF_DATA",
      resource   = "data",
      flowRef    = flowref,
      key        = key,
      start      = startPeriod
    )
    as.data.frame(obj)
  }
  get_safe <- purrr::possibly(get_one, otherwise = NULL)
  raw_list <- countries %>% set_names() %>% map(get_safe)
  
  raw_list %>%
    compact() %>%
    bind_rows() %>%
    transmute(
      country = COUNTRY,
      date    = TIME_PERIOD,                           # "YYYY-MM" ou "YYYY-Qn"
      value   = suppressWarnings(as.numeric(OBS_VALUE))
    ) %>%
    arrange(country, date)
}


#FUNÇÃO 004 : Split tables - dados FMI 
split_indicator_dfs <- function(df_long, indicator_suffix, code_map) {
  dfs <- imap(code_map, function(cc2, cc3) {
    df <- df_long %>%
      filter(country == cc3) %>%
      select(date, value) %>%
      arrange(date)
    
    if (nrow(df) == 0) return(NULL)
    
    # se TODOS os períodos estiverem no formato mensal "YYYY-Mmm", converte p/ Date
    is_all_monthly <- all(grepl("^\\d{4}-M\\d{2}$", df$date))
    
    if (is_all_monthly) {
      df <- df %>%
        mutate(date = as.Date(as.yearmon(sub("-M", "-", date), format = "%Y-%m")))
    } else {
      # mantém como string (ex.: trimestral "YYYY-Qn")
      df$date <- as.character(df$date)
    }
    
    df  # retorna com apenas 'date' (Date ou chr) e 'value' (num)
  }) %>% compact()
  
  # nomes dos objetos no ambiente: br_exchange, us_inflation, cn_gdp, etc.
  names(dfs) <- paste0(unname(code_map[names(dfs)]), "_", indicator_suffix)
  
  list2env(dfs, envir = .GlobalEnv)
  invisible(dfs)
}

#FUNÇÃO 005 : Desagrega PIB trimestral (gdp_q) em mensal usando Prod. Industrial mensal (ip_m)
# helper: "YYYY-Qn" -> Date (primeiro dia do trimestre)
as_date_quarter <- function(x) as.Date(zoo::as.yearqtr(x, format = "%Y-Q%q"))
# helper: time(ts) mensal -> Date
date_from_ts_time_month <- function(tt) as.Date(zoo::as.yearmon(tt))

disagg_gdp_monthly <- function(gdp_q, ip_m,
                               method = "chow-lin-maxlog",
                               conversion = "sum") {
  stopifnot("date" %in% names(gdp_q), "date" %in% names(ip_m))
  
  # ordenar
  gdp_q <- gdp_q %>% arrange(date)
  ip_m  <- ip_m  %>% arrange(date)
  
  # identificar coluna de valor (qualquer nome que não seja 'date')
  gdp_val_col <- setdiff(names(gdp_q), "date")[1]
  ip_val_col  <- setdiff(names(ip_m),  "date")[1]
  if (is.na(gdp_val_col) || is.na(ip_val_col)) {
    stop("Cada tabela precisa ter 'date' e uma coluna de valores.")
  }
  
  # garantir tipos/parse
  # gdp trimestral: date string "YYYY-Qn" -> Date do 1º dia do trimestre + numeric value
  gdp_q <- gdp_q %>%
    mutate(
      date_q = as_date_quarter(date),
      value  = as.double(.data[[gdp_val_col]])
    )
  
  # ip mensal: date em Date (se vier string, tenta converter) + numeric
  if (!inherits(ip_m$date, "Date")) ip_m <- ip_m %>% mutate(date = as.Date(date))
  ip_m <- ip_m %>%
    mutate(
      date_m = date,
      value  = as.double(.data[[ip_val_col]])
    )
  
  if (nrow(gdp_q) == 0) stop("PIB trimestral vazio após parsing.")
  if (nrow(ip_m)  == 0) stop("IP mensal vazio após parsing.")
  
  # ---------------------------
  # AUTO-DETECÇÃO DE STARTS
  # ---------------------------
  # gdp: extrai ano e trimestre do menor período
  gdp_first_str <- gdp_q$date[1]  # presume "YYYY-Qn"
  # fallback robusto caso 'date' já não esteja nesse formato:
  if (!grepl("^\\d{4}-Q[1-4]$", gdp_first_str)) {
    stop("A coluna 'date' do PIB precisa estar no formato 'YYYY-Qn'. Ex.: '2020-Q1'.")
  }
  gdp_year  <- as.integer(substr(gdp_first_str, 1, 4))
  gdp_quart <- as.integer(sub(".*Q", "", gdp_first_str))
  start_gdp <- c(gdp_year, gdp_quart)
  
  # ip: menor data (Date) -> ano e mês
  ip_first_date <- min(ip_m$date_m, na.rm = TRUE)
  start_ip <- c(as.integer(format(ip_first_date, "%Y")),
                as.integer(format(ip_first_date, "%m")))
  
  # ---------------------------
  # ts objects
  # ---------------------------
  gdp_ts <- ts(gdp_q$value, start = start_gdp, frequency = 4)
  ip_ts  <- ts(ip_m$value,  start = start_ip,  frequency = 12) %>%
    forecast::na.interp()
  
  # ---------------------------
  # temporal disaggregation
  # ---------------------------
  gdp_m_fit <- tempdisagg::td(gdp_ts ~ ip_ts, to = "month",
                              method = method, conversion = conversion) %>%
    predict()
  
  # montar saída mensal com datas corretas
  gdp_monthly <- tibble(
    date  = date_from_ts_time_month(time(gdp_m_fit)),
    value = as.numeric(gdp_m_fit)
  ) %>%
    # cortar para o intervalo observado do PIB trimestral
    filter(date >= min(gdp_q$date_q), date <= max(gdp_q$date_q))
  
  # nome da coluna de saída: igual ao nome do PIB trimestral, removendo sufixo "_q" se existir
  out_colname <- sub("_q$", "", gdp_val_col)
  names(gdp_monthly)[names(gdp_monthly) == "value"] <- out_colname
  
  gdp_monthly
}

#FUNÇÃO 006 : Interpolação Linear para dados faltantes - Um Dataframe
interp_linear_df_any <- function(df, value_cols = NULL, maxgap = Inf, extrapolate = FALSE) {
  stopifnot("date" %in% names(df))
  
  df <- df %>% arrange(date) %>%
    mutate(date = if (!inherits(date, "Date")) as.Date(date) else date)
  
  if (is.null(value_cols)) value_cols <- setdiff(names(df), "date")
  if (length(value_cols) == 0) return(df)
  
  # ── CORREÇÃO: expandir para sequência mensal completa ──────────────────────
  full_dates <- seq(min(df$date), max(df$date), by = "month")
  df <- full_join(data.frame(date = full_dates), df, by = "date")
  # ───────────────────────────────────────────────────────────────────────────
  
  df %>%
    mutate(across(
      all_of(value_cols),
      ~ zoo::na.approx(
        object = suppressWarnings(as.numeric(.x)),
        x      = as.numeric(date),
        na.rm  = extrapolate,
        maxgap = maxgap
      ),
      .names = "{.col}"
    ))
}
# interp_linear_df_any <- function(df, value_cols = NULL, maxgap = Inf, extrapolate = FALSE) {
#   stopifnot("date" %in% names(df))
#   df <- df %>% arrange(date)
#   
#   # detectar colunas de valor, se não informadas
#   if (is.null(value_cols)) value_cols <- setdiff(names(df), "date")
#   if (length(value_cols) == 0) return(df)
#   
#   df %>%
#     mutate(date = if (!inherits(date, "Date")) as.Date(date) else date) %>%
#     mutate(across(
#       all_of(value_cols),
#       ~ zoo::na.approx(
#         # >>> CORREÇÃO: o primeiro argumento é o vetor a interpolar (object)
#         object = suppressWarnings(as.numeric(.x)),
#         x       = as.numeric(date),
#         na.rm   = extrapolate,   # FALSE = sem extrapolar bordas
#         maxgap  = maxgap
#       ),
#       .names = "{.col}"
#     ))
# }

#FUNÇÃO 007 : Interpolação Linear para dados faltantes - Vários Dataframes
# (NAO MAIS USADA pelo 02 -- substituida por fill_locf_many(), problema 4)
interp_linear_many_any <- function(lst, maxgap = Inf, extrapolate = FALSE) {
  imap(lst, ~ interp_linear_df_any(.x, value_cols = setdiff(names(.x), "date"),
                                   maxgap = maxgap, extrapolate = extrapolate))
}

#FUNÇÃO 007b : Preenchimento de faltantes em TEMPO REAL (ultima observacao)
#
# CORRECAO: problema 4 de review_problems.md. zoo::na.approx() interpolava
# lacunas internas em linha reta ate o PROXIMO valor observado (informacao
# futura). Agora cada lacuna interna recebe o ultimo valor observado antes
# dela. A serie e expandida para o calendario mensal completo entre a primeira
# e a ultima observacao; nao ha extrapolacao antes do inicio.
fill_locf_df <- function(df) {
  stopifnot("date" %in% names(df))
  df <- df %>% arrange(date) %>%
    mutate(date = if (!inherits(date, "Date")) as.Date(date) else date)
  value_cols <- setdiff(names(df), "date")
  if (length(value_cols) == 0 || nrow(df) == 0) return(df)
  full_dates <- seq(min(df$date), max(df$date), by = "month")
  df <- full_join(data.frame(date = full_dates), df, by = "date") %>% arrange(date)
  df %>% mutate(across(all_of(value_cols),
                       ~ zoo::na.locf(suppressWarnings(as.numeric(.x)), na.rm = FALSE)))
}
fill_locf_many <- function(lst) lapply(lst, fill_locf_df)

#FUNÇÃO 007c : PIB mensal em TEMPO REAL via Chow-Lin com producao industrial
#
# CORRECAO: problema 4 de review_problems.md.
#  (a) Bug: o 02 converte a data trimestral do FMI para Date, mas
#      disagg_gdp_monthly() exige o texto "AAAA-Qn"; a funcao parava e o
#      tryCatch do 02 caia sempre em disagg_gdp_fallback() (interpolacao
#      linear). Chow-Lin nunca era usado, e o PIB do mes m interpolava em
#      direcao ao trimestre seguinte (5-7 meses de informacao futura).
#  (b) Tempo real: para cada mes t (origem), usa somente
#        - trimestres de PIB ja publicados: fim do trimestre + gdp_lag <= t
#        - producao industrial ja publicada: mes + ip_lag <= t
#      roda Chow-Lin (tempdisagg::td, conversao "sum") nesse vintage e guarda a
#      estimativa mensal mais recente do vintage (o ultimo mes com PI). Esse e
#      o valor do PIB em t: a melhor estimativa disponivel na origem t.
#  (c) Quando o vintage tem menos de min_quarters trimestres cobertos pela PI
#      (ou a PI termina antes do ultimo trimestre publicado, ou o Chow-Lin
#      falha), o valor em t e o ultimo trimestre publicado / 3 (mesma escala
#      mensal da conversao "sum"). A fonte de cada valor fica na coluna source.
#
# gdp_q: data.frame (date = Date do 1o mes do trimestre, valor na 2a coluna)
# ip_m : data.frame (date = Date mensal, valor na 2a coluna) ou NULL
realtime_gdp_monthly <- function(gdp_q, ip_m, gdp_lag = 3, ip_lag = 2,
                                 min_quarters = 12, method = "chow-lin-maxlog",
                                 end_date = NULL) {
  # soma k meses a datas do 1o dia do mes
  add_m <- function(d, k) {
    m <- as.integer(format(d, "%m")) - 1L + as.integer(k)
    as.Date(sprintf("%04d-%02d-01", as.integer(format(d, "%Y")) + m %/% 12L, m %% 12L + 1L))
  }
  gq <- data.frame(date = as.Date(gdp_q[[1]]), value = as.numeric(gdp_q[[2]])) %>%
    filter(!is.na(value)) %>% arrange(date)
  gq$q_end   <- add_m(gq$date, 2)
  gq$pub     <- add_m(gq$q_end, gdp_lag)

  have_ip <- !is.null(ip_m) && nrow(ip_m) > 0
  if (have_ip) {
    im <- data.frame(date = as.Date(ip_m[[1]]), value = as.numeric(ip_m[[2]])) %>%
      filter(!is.na(value)) %>% arrange(date)
    im$pub <- add_m(im$date, ip_lag)
  }

  if (is.null(end_date)) {
    end_date <- max(gq$pub, if (have_ip) max(im$pub))
  }
  origins <- seq(min(gq$pub), as.Date(end_date), by = "month")
  value   <- rep(NA_real_, length(origins))
  source  <- rep(NA_character_, length(origins))

  for (k in seq_along(origins)) {
    t  <- origins[k]
    qa <- gq[gq$pub <= t, ]
    if (nrow(qa) == 0) next
    step_value <- qa$value[nrow(qa)] / 3

    est <- NULL
    if (have_ip) {
      ia <- im[im$pub <= t, ]
      if (nrow(ia) > 0) {
        # trimestres inteiramente cobertos pela PI disponivel (sem lacunas)
        ip_start <- min(ia$date); ip_end <- max(ia$date)
        qc <- qa[qa$date >= ip_start & qa$q_end <= ip_end, ]
        if (nrow(qc) >= min_quarters && max(qc$q_end) == max(qa$q_end)) {
          months_needed <- seq(min(qc$date), ip_end, by = "month")
          iv <- ia$value[match(months_needed, ia$date)]
          if (!anyNA(iv)) {
            est <- tryCatch({
              g_ts <- stats::ts(qc$value,
                                start = c(as.integer(format(qc$date[1], "%Y")),
                                          (as.integer(format(qc$date[1], "%m")) - 1) %/% 3 + 1),
                                frequency = 4)
              i_ts <- stats::ts(iv,
                                start = c(as.integer(format(months_needed[1], "%Y")),
                                          as.integer(format(months_needed[1], "%m"))),
                                frequency = 12)
              fit <- stats::predict(tempdisagg::td(g_ts ~ i_ts, to = "month",
                                                   method = method, conversion = "sum"))
              as.numeric(fit)[length(fit)]
            }, error = function(e) NULL)
          }
        }
      }
    }
    if (!is.null(est) && is.finite(est) && est > 0) {
      value[k] <- est; source[k] <- "chow-lin"
    } else {
      value[k] <- step_value; source[k] <- "last quarter / 3"
    }
  }
  data.frame(date = origins, value = value, source = source)
}

#FUNÇÃO 007d : PIB mensal por Chow-Lin estimado UMA VEZ na amostra inteira
# Variante "filters" (M1 de reviews_v7.md): reproduz o vazamento de informacao
# da desagregacao estimada com todos os trimestres. Para isolar esse canal, o
# valor em cada origem t e a estimativa do mes t - ip_lag (o mesmo mes que o
# realtime_gdp_monthly() usa em t) e as datas e o fallback sao os do tempo real:
# so os meses em que o tempo real usa Chow-Lin recebem o valor de amostra inteira.
fullsample_gdp_monthly <- function(rt, gdp_q, ip_m, ip_lag = 2, method = "chow-lin-maxlog") {
  if (is.null(ip_m) || nrow(ip_m) == 0) return(rt)
  add_m <- function(d, k) {
    m <- as.integer(format(d, "%m")) - 1L + as.integer(k)
    as.Date(sprintf("%04d-%02d-01", as.integer(format(d, "%Y")) + m %/% 12L, m %% 12L + 1L))
  }
  gq <- data.frame(date = as.Date(gdp_q[[1]]), value = as.numeric(gdp_q[[2]])) %>%
    filter(!is.na(value)) %>% arrange(date)
  im <- data.frame(date = as.Date(ip_m[[1]]), value = as.numeric(ip_m[[2]])) %>%
    filter(!is.na(value)) %>% arrange(date)
  ip_start <- min(im$date); ip_end <- max(im$date)
  qc <- gq[gq$date >= ip_start & add_m(gq$date, 2) <= ip_end, ]
  if (nrow(qc) < 12) return(rt)
  months_needed <- seq(min(qc$date), ip_end, by = "month")
  iv <- im$value[match(months_needed, im$date)]
  if (anyNA(iv)) iv <- zoo::na.locf(iv, na.rm = FALSE)
  est <- tryCatch({
    g_ts <- stats::ts(qc$value, start = c(as.integer(format(qc$date[1], "%Y")),
                                         (as.integer(format(qc$date[1], "%m")) - 1) %/% 3 + 1),
                      frequency = 4)
    i_ts <- stats::ts(iv, start = c(as.integer(format(months_needed[1], "%Y")),
                                    as.integer(format(months_needed[1], "%m"))), frequency = 12)
    as.numeric(stats::predict(tempdisagg::td(g_ts ~ i_ts, to = "month", method = method,
                                             conversion = "sum")))
  }, error = function(e) NULL)
  if (is.null(est)) return(rt)
  full <- data.frame(date = add_m(months_needed[seq_along(est)], ip_lag), full = est)
  out  <- rt
  idx  <- match(out$date, full$date)
  use  <- out$source == "chow-lin" & !is.na(idx)
  out$value[use]  <- full$full[idx[use]]
  out$source[use] <- "chow-lin (full sample)"
  out
}

#FUNÇÃO 007e : Variantes da rodada (M1 de reviews_v7.md)
# RUN_VARIANT (variavel de ambiente, definida em cada .bat):
#   main    desenho sem vazamento (resultados principais)
#   window  janela de estimacao termina em t-1 em vez de t-h
#   nolags  sem defasagens de publicacao (CPI, M1, PI e PIB no mes de referencia)
#   filters Chow-Lin e filtro de Hamilton estimados uma vez na amostra inteira
#   scaling padronizacao e nos das splines calculados na amostra inteira
# Cada variante viola uma unica condicao de tempo real; o resto e igual ao main.
run_variants <- c("main", "window", "nolags", "filters", "scaling")
run_variant <- function() {
  v <- tolower(trimws(Sys.getenv("RUN_VARIANT", "main")))
  if (!nzchar(v)) v <- "main"
  if (!v %in% run_variants) stop("RUN_VARIANT invalido: ", v, " (use ", paste(run_variants, collapse = ", "), ")")
  v
}
# Paineis: o 02 grava outputs/<pais>_panel_data<sufixo>.csv; window e scaling
# usam o painel do main (a violacao delas esta no 03).
variant_panel_suffix <- function(v = run_variant()) switch(v, nolags = "_nolags", filters = "_filters", "")
panel_path <- function(cty, v = run_variant()) {
  file.path("outputs", paste0(cty, "_panel_data", variant_panel_suffix(v), ".csv"))
}
# Resultados (04, 05, 06, aggregate_stats): outputs/ no main, outputs_<variante>/ nas demais
variant_out_dir <- function(v = run_variant()) if (v == "main") "outputs" else paste0("outputs_", v)

# M4 de reviews_v7.md: a India e avaliada (04, 06) sem o modelo de Taylor. O hiato
# do produto so existe a partir de ~2010 e limitava a janela comum da India a
# 2020-2026; sem ele, a janela comum e a dos quatro modelos monetarios. O 03
# continua estimando o modelo de Taylor (fica no .rds, fora da avaliacao).
eval_excluded_models <- list(`in` = "taylor")
eval_models <- function(cty, model_list) setdiff(model_list, eval_excluded_models[[cty]])

# FUNÇÃO 008 : (removida -- ponto 29 de review_problems.md) add_coint_term()
# estimava um termo cointegrante (Johansen) padronizado com a amostra inteira;
# nao era usado por nenhum script.
