##---------------
# Clears all old variables
rm(list=ls())

##---------------
# Sets the directory automatically (point 22: works with Rscript/.bat and in RStudio)
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

# _____________________________________________________________________________
# _______________________________ HELPER FUNCTIONS ____________________________
# _____________________________________________________________________________

source("01_functions.R")
library(tidyverse)
library(lubridate)
fredr_set_key(read_lines("fred_key.txt"))
set.seed(123)

# _____________________________________________________________________________
# _______________________ CONJUNTO DE PAISES (ponto 36) _______________________
# _____________________________________________________________________________
# Argumento opcional: "brics" (padrao, o artigo) ou "cp" (paises de Colombo e
# Pelagatti, 2020, que tem todas as series das mesmas fontes: Canada, Japao,
# Reino Unido e Franca). No conjunto "cp" os dados comecam em 1973, o inicio da
# amostra de Colombo e Pelagatti; cada serie termina quando a fonte termina
# (o franco frances termina em 12/1998).
country_set <- {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) >= 1 && nzchar(a[1])) tolower(trimws(a[1])) else "brics"
}
# Ponto 39: "cp" passa a incluir Dinamarca e Suecia (seis paises); "cp_dkse" baixa
# so os dois (mais os EUA), para acrescenta-los sem refazer os paineis ja estimados.
if (!country_set %in% c("brics", "cp", "cp_dkse")) stop("Conjunto invalido: ", country_set, " (use brics, cp ou cp_dkse)")
start_date <- if (country_set %in% c("cp", "cp_dkse")) "1973-01-01" else "1996-01-01"
cat("02 | conjunto de paises:", country_set, "| inicio:", start_date, "\n")

# _____________________________________________________________________________
# _________________ DEFASAGENS DE PUBLICACAO (problema 4, tempo real) _________
# _____________________________________________________________________________
# Na origem t (fim do mes t) o previsor conhece:
#   - cambio e juros do proprio mes t (sem defasagem);
#   - CPI e M1 ate t-1 (divulgados cerca de um mes depois);
#   - producao industrial ate t-2;
#   - PIB dos trimestres encerrados ate t-3 (ex.: o 1o trimestre, encerrado em
#     marco, e usado a partir de junho).
cpi_m1_pub_lag <- 1   # meses
ip_pub_lag     <- 2   # meses
gdp_pub_lag    <- 3   # meses apos o fim do trimestre

# Variantes (M1 de reviews_v7.md; RUN_VARIANT definido no .bat, ver 01_functions.R).
# O 02 so muda nas variantes "nolags" e "filters"; as demais usam o painel do main.
variant <- run_variant()
if (variant == "nolags") {
  # Violacao: cada serie entra no mes de referencia (PIB no ultimo mes do trimestre)
  cpi_m1_pub_lag <- 0; ip_pub_lag <- 0; gdp_pub_lag <- 0
}
cat("02 | variante:", variant, "| defasagens CPI/M1, PI, PIB:", cpi_m1_pub_lag, ip_pub_lag, gdp_pub_lag, "\n")

# Dados brutos guardados pelo main (outputs/raw_<conjunto>.rds): as variantes
# reusam exatamente os mesmos dados, sem baixar de novo (as fontes podem revisar
# as series entre uma rodada e outra).
raw_cache  <- file.path("outputs", paste0("raw_", country_set, ".rds"))
use_cache  <- variant != "main"
if (use_cache && !file.exists(raw_cache)) stop("Rode antes o 02 do main para este conjunto (falta ", raw_cache, ")")
# TCC forecast: get data
# Wagner Scheid Mallmann
# 2025/02

# TCC forecast: get data
# Wagner Scheid Mallmann
# 2025/02


# Brasil, China, Indonésia, Índia, Rússia e África do Sul (além dos EUA como base).

# Países 100% COMPLETOS (Aprovados): BR, CN, ID, IN, RU, US, ZA
# 
# Países INCOMPLETOS (Excluir): AE, EG, ET, IR, SA

# _____________________________________________________________________________
# ________________________________ GET DATA IMF________________________________
# _____________________________________________________________________________

# ---------------------- PARÂMETROS E DOWNLOAD ----------------------
# Apenas países aprovados no diagnóstico de dados + EUA (benchmark)
if (country_set == "cp") {
  countries_imf <- c("CAN", "JPN", "GBR", "FRA", "DNK", "SWE", "USA")
  code_map <- c(CAN="ca", JPN="jp", GBR="gb", FRA="fr", DNK="dk", SWE="se", USA="us")
} else if (country_set == "cp_dkse") {
  countries_imf <- c("DNK", "SWE", "USA")
  code_map <- c(DNK="dk", SWE="se", USA="us")
} else {
  countries_imf <- c("BRA", "CHN", "IDN", "IND", "RUS", "ZAF", "USA")
  code_map <- c(BRA="br", CHN="cn", IDN="id", IND="in", RUS="ru", ZAF="za", USA="us")
}

ccodes <- unname(code_map)

if (use_cache) {
  raw <- readRDS(raw_cache)
  list2env(raw, envir = .GlobalEnv)
  cat("02 | dados brutos lidos de", raw_cache, "\n")
} else {

flow_er   <- "IMF.STA,ER"   ; patt_er   <- "XDC_USD.EOP_RT.M"
flow_cpi  <- "IMF.STA,CPI"  ; patt_cpi  <- "CPI._T.IX.M"          
flow_gdp  <- "IMF.STA,QNEA" ; patt_gdp  <- "B1GQ.V.NSA.XDC.Q"     

er_long  <- fetch_imf_series(flow_er,  patt_er,  countries_imf, startPeriod = start_date)
cpi_long <- fetch_imf_series(flow_cpi, patt_cpi, countries_imf, startPeriod = start_date)
gdp_long <- fetch_imf_series(flow_gdp, patt_gdp, countries_imf, startPeriod = start_date)

split_indicator_dfs(er_long,  "exchange",  code_map)   
split_indicator_dfs(cpi_long, "inflation", code_map)   
split_indicator_dfs(gdp_long, "gdp_q",     code_map)   

# --- CORREÇÃO DE DATAS DO FMI ---
parse_imf_date <- function(dates) {
  if (inherits(dates, "Date")) return(dates)
  dates <- as.character(dates)
  
  if (any(grepl("Q", dates, ignore.case = TRUE))) {
    y <- substr(dates, 1, 4)
    q <- as.numeric(sub(".*Q", "", dates))
    m <- (q - 1) * 3 + 1
    return(as.Date(sprintf("%s-%02d-01", y, m)))
  }
  if (all(nchar(dates[!is.na(dates)]) == 7)) {
    return(as.Date(paste0(dates, "-01")))
  }
  return(as.Date(dates))
}

ccodes <- unname(code_map)

for(cc in ccodes) {
  for(var in c("exchange", "inflation", "gdp_q")) {
    df_name <- paste0(cc, "_", var)
    if(exists(df_name)) {
      df <- get(df_name)
      if(nrow(df) > 0) {
        colnames(df)[1] <- "date" 
        df$date <- parse_imf_date(df$date)
        assign(df_name, df, envir = .GlobalEnv)
      }
    }
  }
}

# _____________________________________________________________________________
# ________________________________ GET DATA FRED ______________________________
# _____________________________________________________________________________

safe_get_fred <- function(possible_ids, obs_start = start_date) {
  for (id in possible_ids) {
    if (is.na(id)) next
    res <- tryCatch({
      df <- fredr(series_id = id, observation_start = as.Date(obs_start))
      if (nrow(df) > 0) return(df %>% select(date, value))
      NULL
    }, error = function(e) NULL)
    if (!is.null(res)) return(res)
  }
  return(data.frame(date = as.Date(character()), value = numeric()))
}

fred_interest <- list(
  br = "IRSTCI01BRM156N", ru = "IRSTCI01RUM156N", `in` = "IRSTCI01INM156N", 
  cn = "IRSTCI01CNM156N", za = "IRSTCI01ZAM156N", id = "IRSTCI01IDM156N",
  us = "IRSTCI01USM156N",
  # conjunto "cp" (ponto 36)
  ca = "IRSTCI01CAM156N", jp = "IRSTCI01JPM156N", gb = "IRSTCI01GBM156N", fr = "IRSTCI01FRM156N"
)

fred_m1 <- list(
  br = "MANMM101BRM189N", ru = "MANMM101RUM189N", `in` = "MANMM101INM189N",
  cn = "MANMM101CNM189N", za = "MANMM101ZAM189N", id = c("MANMM101IDM189N", "MYAGM1IDNM189N"),
  us = "MANMM101USM189N"
)

fred_ip <- list(
  br = "BRAPRMNTO01IXOBM", ru = "RUSPRMNTO01IXOBM", `in` = "INDPRMNTO01IXOBM",
  cn = "CHNPRINTO01IXPYM", za = "ZAFPRMNTO01IXOBM", id = c("IDNPRMNTO01IXOBM", "IDNPROINDISAAA"),
  us = "USAPRMNTO01IXOBM",
  # conjunto "cp" (ponto 36)
  ca = "CANPRMNTO01IXOBM", jp = "JPNPRMNTO01IXOBM", gb = "GBRPRMNTO01IXOBM", fr = "FRAPRMNTO01IXOBM",
  dk = "DNKPRMNTO01IXOBM", se = "SWEPRMNTO01IXOBM"
)

for (cc in ccodes) {
  df_int <- safe_get_fred(fred_interest[[cc]])
  if(nrow(df_int) > 0) colnames(df_int)[2] <- paste0(cc, "_interest")
  assign(paste0(cc, "_interest"), df_int, envir = .GlobalEnv)
  
  # df_m1 <- safe_get_fred(fred_m1[[cc]])
  # if(nrow(df_m1) > 0) colnames(df_m1)[2] <- paste0(cc, "_m1")
  # assign(paste0(cc, "_m1"), df_m1, envir = .GlobalEnv)

  df_ip <- safe_get_fred(fred_ip[[cc]])
  if(nrow(df_ip) > 0) colnames(df_ip)[2] <- paste0(cc, "_ip")
  assign(paste0(cc, "_ip"), df_ip, envir = .GlobalEnv)
}

# Juros de curto prazo da Dinamarca e da Suecia (ponto 39): a FRED nao tem a
# taxa call (IRSTCI01) desses paises. Usa a taxa interbancaria overnight
# (day-to-day) da Eurostat, irt_st_m, via DBnomics: o mesmo tipo de taxa, outra
# fonte. Dinamarca desde 01/1997, Suecia desde 01/1985 (medias mensais).
dbnomics_interest <- list(dk = "Eurostat/irt_st_m/M.IRT_DTD.DK", se = "Eurostat/irt_st_m/M.IRT_DTD.SE")
get_dbnomics <- function(id) {
  j <- jsonlite::fromJSON(paste0("https://api.db.nomics.world/v22/series/", id, "?observations=1&format=json"))
  s <- j$series$docs
  d <- data.frame(date  = as.Date(paste0(unlist(s$period), "-01")),
                  value = suppressWarnings(as.numeric(unlist(s$value))))
  d[!is.na(d$value) & d$date >= as.Date(start_date), ]
}
for (cc in intersect(names(dbnomics_interest), ccodes)) {
  if (nrow(get(paste0(cc, "_interest"))) == 0) {
    df_int <- get_dbnomics(dbnomics_interest[[cc]])
    if (nrow(df_int) > 0) colnames(df_int)[2] <- paste0(cc, "_interest")
    assign(paste0(cc, "_interest"), df_int, envir = .GlobalEnv)
  }
}

# _____________________________________________________________________________
# ________________________________ GET DATA OECD M1 ___________________________
# _____________________________________________________________________________

# Mapeamento: código OECD → código interno do script
oecd_country_map <- c(
  "Brazil"        = "br",
  "China"         = "cn",
  "Indonesia"     = "id",
  "India"         = "in",
  "Russia"        = "ru",
  "South Africa"  = "za",
  "United States" = "us",
  # conjunto "cp" (ponto 36); a Franca nao tem M1 nacional no arquivo da OCDE
  "Canada"         = "ca",
  "Japan"          = "jp",
  "United Kingdom" = "gb",
  "Denmark"        = "dk",
  "Sweden"         = "se"
)

# Leitura do CSV da OCDE
oecd_raw <- read_csv(
  "OECD.SDD.STES,DSD_STES@DF_MONAGG,4.0+all.csv",
  show_col_types = FALSE
)

# Filtra: M1, frequência mensal, países de interesse, sem transformação
oecd_m1_raw <- oecd_raw %>%
  filter(
    Measure                    == "M1",
    `Frequency of observation` == "Monthly",
    `Reference area`           %in% names(oecd_country_map),
    `Decimals`                 == 'Zero',
    `Adjustment`               == 'Neither seasonally adjusted nor calendar adjusted'
  ) %>%
  mutate(
    date  = if_else(
      grepl("^\\d{4}-\\d{2}$", TIME_PERIOD),
      as.Date(paste0(TIME_PERIOD, "-01")),
      as.Date(NA)
    ),
    value = suppressWarnings(as.numeric(OBS_VALUE)),
    cc    = oecd_country_map[`Reference area`]
  ) %>%
  filter(!is.na(date), !is.na(value), date >= as.Date(start_date)) %>%
  select(cc, date, value)

# Cria objetos cc_m1 no GlobalEnv
for (cc in ccodes) {
  df_m1 <- oecd_m1_raw %>%
    filter(cc == !!cc) %>%
    select(date, value)

  if (nrow(df_m1) > 0) {
    colnames(df_m1)[2] <- paste0(cc, "_m1")
  }

  assign(paste0(cc, "_m1"), df_m1, envir = .GlobalEnv)
}

# M1 da Franca (ponto 36): o arquivo da OCDE so tem o agregado da area do euro.
# Usa a serie M1 do IFS (FMI) arquivada na FRED, em francos, 12/1977 a 12/1998.
fred_m1_fallback <- list(fr = "MYAGM1FRM189N")
for (cc in intersect(names(fred_m1_fallback), ccodes)) {
  if (nrow(get(paste0(cc, "_m1"))) == 0) {
    df_m1 <- safe_get_fred(fred_m1_fallback[[cc]])
    if (nrow(df_m1) > 0) colnames(df_m1)[2] <- paste0(cc, "_m1")
    assign(paste0(cc, "_m1"), df_m1, envir = .GlobalEnv)
  }
}

# Guarda os dados brutos para as variantes
if (!dir.exists("outputs")) dir.create("outputs")
raw_names <- as.vector(outer(ccodes, c("exchange", "inflation", "gdp_q", "interest", "ip", "m1"),
                             paste, sep = "_"))
raw_names <- raw_names[sapply(raw_names, exists, envir = .GlobalEnv)]
saveRDS(mget(raw_names, envir = .GlobalEnv), raw_cache)
cat("02 | dados brutos guardados em", raw_cache, "\n")
}   # fim do download (main)

# _____________________________________________________________________________
# ______________________ FALTANTES: ULTIMA OBSERVACAO (problema 4) ____________
# _____________________________________________________________________________
# Antes: zoo::na.approx() (interpolacao linear ate o proximo valor observado,
# informacao futura). Agora: cada lacuna interna recebe o ultimo valor observado.

lista_interpolacao <- list()
for(cc in ccodes) {
  for(var in c("exchange", "ip", "inflation", "interest", "m1")) {
    df_name <- paste0(cc, "_", var)
    if(exists(df_name)) {
      df <- get(df_name)
      if(nrow(df) > 0) lista_interpolacao[[df_name]] <- df
    }
  }
}

lst_interp <- fill_locf_many(lista_interpolacao)

for(nome in names(lst_interp)) {
  assign(nome, lst_interp[[nome]], envir = .GlobalEnv)
}

# _____________________________________________________________________________
# ____________________ CPI E M1: DEFASAGEM DE PUBLICACAO (problema 4) _________
# _____________________________________________________________________________
# O valor do mes m so e conhecido em m + cpi_m1_pub_lag: a observacao e
# re-datada para o mes em que ja esta disponivel.
for (cc in ccodes) {
  for (var in c("inflation", "m1")) {
    df_name <- paste0(cc, "_", var)
    if (exists(df_name)) {
      df <- get(df_name)
      if (nrow(df) > 0) {
        df$date <- df$date %m+% months(cpi_m1_pub_lag)
        assign(df_name, df, envir = .GlobalEnv)
      }
    }
  }
}

# _____________________________________________________________________________
# ________________ PIB MENSAL EM TEMPO REAL: CHOW-LIN COM PI (problema 4) ______
# _____________________________________________________________________________
# Antes: disagg_gdp_monthly() nunca rodava (bug de formato de data) e todos os
# paises caiam em disagg_gdp_fallback(), uma interpolacao linear entre
# trimestres que usava o trimestre seguinte. Agora: realtime_gdp_monthly()
# (01_functions.R) roda Chow-Lin, para cada mes t, so com o PIB e a PI ja
# publicados em t, e guarda a estimativa mais recente desse vintage.

gdp_source_log <- list()
for (cc in ccodes) {
  gdp_q_name <- paste0(cc, "_gdp_q")
  ip_name    <- paste0(cc, "_ip")
  gdp_q <- if (exists(gdp_q_name, envir = .GlobalEnv)) get(gdp_q_name, envir = .GlobalEnv) else NULL
  ip_m  <- if (exists(ip_name, envir = .GlobalEnv)) get(ip_name, envir = .GlobalEnv) else NULL
  if (is.null(gdp_q) || nrow(gdp_q) == 0) next

  rt <- realtime_gdp_monthly(gdp_q, ip_m, gdp_lag = gdp_pub_lag, ip_lag = ip_pub_lag)
  # Variante "filters": Chow-Lin estimado uma vez na amostra inteira (mesmas datas)
  if (variant == "filters") rt <- fullsample_gdp_monthly(rt, gdp_q, ip_m, ip_lag = ip_pub_lag)
  gdp_source_log[[cc]] <- rt %>% mutate(country = cc)
  assign(paste0(cc, "_gdp"), rt %>% select(date, value), envir = .GlobalEnv)
}

if (!dir.exists("outputs")) dir.create("outputs")
gdp_source_all <- bind_rows(gdp_source_log)
# (conjunto cp_dkse: arquivo proprio, para nao apagar o diagnostico do conjunto cp)
write_csv2(gdp_source_all, paste0(if (country_set == "cp_dkse") "outputs/gdp_realtime_source_dkse"
                                  else paste0("outputs/gdp_realtime_source_", country_set),
                                  variant_panel_suffix(variant), ".csv"))
print(gdp_source_all %>% count(country, source) %>% tidyr::pivot_wider(names_from = source, values_from = n))

# _____________________________________________________________________________
# ______________________  OUTPUT GAP VIA HAMILTON FILTER ______________________
# _____________________________________________________________________________

for (cc in ccodes) {
  gdp_m <- get(paste0(cc, "_gdp"))
  if(!is.null(gdp_m) && nrow(gdp_m) > 0) {
    gap_df <- tryCatch(output_gap_m(data = gdp_m, colname = paste0(cc, "_gap")), error = function(e) NULL)
    # Variante "filters": regressao de Hamilton estimada uma vez na amostra inteira
    # (min_obs = Inf). Para isolar o vazamento, o hiato fica NA nos mesmos meses
    # em que o hiato recursivo e NA (mesma amostra do main).
    if (variant == "filters" && !is.null(gap_df)) {
      full_df <- tryCatch(output_gap_m(data = gdp_m, colname = paste0(cc, "_gap"), min_obs = Inf),
                          error = function(e) NULL)
      if (!is.null(full_df)) {
        full_df[[2]][is.na(gap_df[[2]])] <- NA
        gap_df <- full_df
      }
    }
    if(!is.null(gap_df)) assign(paste0(cc, "_gap"), gap_df, envir = .GlobalEnv)
  }
}

# _____________________________________________________________________________
# ______________________  RENOMEAR COLUNAS E LOGS _____________________________
# _____________________________________________________________________________

apply_log <- function(df_name, new_col_name) {
  if (exists(df_name, envir = .GlobalEnv)) {
    df <- get(df_name, envir = .GlobalEnv)
    if(nrow(df) > 0) {
      base_col <- names(df)[2]
      colnames(df)[2] <- sub("_log", "", new_col_name) 
      df[[new_col_name]] <- log(df[[2]])
      assign(df_name, df, envir = .GlobalEnv)
    }
  }
}

for(cc in ccodes) {
  apply_log(paste0(cc, "_exchange"), paste0(cc, "_et"))
  apply_log(paste0(cc, "_m1"), paste0(cc, "_m1_log"))
  apply_log(paste0(cc, "_gdp"), paste0(cc, "_gdp_log"))
  apply_log(paste0(cc, "_inflation"), paste0(cc, "_inflation_log"))
}

# _____________________________________________________________________________
# ______________________________  EXCHANGE RATE _______________________________
# _____________________________________________________________________________

make_exchange_exprs <- function(cc, et_col) {
  exprs_list <- vector("list", 4L * 12L)
  idx <- 1L
  for (h in 1:12) {
    exprs_list[[idx]]       <- eval(parse(text = sprintf("quote(dplyr::lead(%s, %d) - %s)", et_col, h, et_col)))
    names(exprs_list)[idx]  <- sprintf("%s_et_dif_h%d", cc, h)
    idx <- idx + 1L
    exprs_list[[idx]]       <- eval(parse(text = sprintf("quote(dplyr::lead(%s, %d) - %s)", et_col, h, et_col)))
    names(exprs_list)[idx]  <- sprintf("%s_et_ecm_h%d", cc, h)
    idx <- idx + 1L
  }
  exprs_list
}

# Exclui EUA da iteração do painel (benchmark)
exchange_ccodes <- setdiff(ccodes, "us")

for (cc in exchange_ccodes) {
  df_name <- paste0(cc, "_exchange")
  if (exists(df_name, envir = .GlobalEnv) && nrow(get(df_name)) > 0) {
    et_col  <- paste0(cc, "_et")
    df      <- get(df_name)
    exprs   <- make_exchange_exprs(cc, et_col)
    df      <- df %>% mutate(!!!exprs)
    assign(df_name, df)
  }
}

# _____________________________________________________________________________
# ______________________ DEFINIÇÃO DA FUNÇÃO make_panel_exprs _________________
# _____________________________________________________________________________

make_panel_exprs <- function(cc) {
  dif_exprs <- list()
  # Define as variáveis de diferença para os 5 fundamentos
  fundamentals <- c("m_hat", "y_hat", "r_hat", "i_hat", "ygap_hat")
  
  for (fund in fundamentals) {
    # Cria a string de referência base (ex: br_m_hat)
    base_name <- paste0(cc, "_", fund)
    
    for (h in 1:12) {
      col_name <- paste0(base_name, "_dif_h", h)
      # Cria a expressão de diferença: variável(t) - variável(t-h)
      # Usamos .data para garantir acesso seguro dentro do mutate
      expr_str <- substitute( .data[[base]] - dplyr::lag(.data[[base]], h), 
                              list(base = base_name, h = h))
      dif_exprs[[col_name]] <- expr_str
    }
  }
  return(dif_exprs)
}

# _____________________________________________________________________________
# ______________________________ BUILD PANEL (VERSÃO ESTÁVEL) _________________
# _____________________________________________________________________________

for (cc in exchange_ccodes) {
  
  # Certifica que os dados base existem
  ref_df <- get(paste0(cc, "_exchange"))
  
  # ... (mantenha o seu código de ensure_df aqui) ...
  
  df <- reduce(
    list(
      get(paste0(cc, "_exchange")), get(paste0(cc, "_gap")), get(paste0(cc, "_gdp")),
      get(paste0(cc, "_inflation")), get(paste0(cc, "_interest")), get(paste0(cc, "_m1")),
      us_exchange, us_gap, us_gdp, us_inflation, us_interest, us_m1
    ),
    full_join, by = "date"
  )
  
  # Mapeamento seguro das variáveis: acessa os objetos pelo nome no GlobalEnv
  # sem precisar de parse() ou eval()
  df <- df %>%
    arrange(date) %>%
    mutate(
      !!paste0(cc, "_m_hat")    := .data[[paste0(cc, "_m1_log")]] - .data[["us_m1_log"]],
      !!paste0(cc, "_y_hat")    := .data[[paste0(cc, "_gdp_log")]] - .data[["us_gdp_log"]],
      !!paste0(cc, "_r_hat")    := .data[[paste0(cc, "_inflation_log")]] - .data[["us_inflation_log"]],
      !!paste0(cc, "_i_hat")    := .data[[paste0(cc, "_interest")]] - .data[["us_interest"]],
      !!paste0(cc, "_ygap_hat") := .data[[paste0(cc, "_gap")]] - .data[["us_gap"]]
    )
  
  # ---- DIF: deltas em t-k, para k = 1:12 ----
  dif_exprs <- make_panel_exprs(cc)
  df <- df %>% mutate(!!!dif_exprs)
  
  assign(paste0(cc, "_panel_data"), df)
}

# CORRECAO (ponto 29 de review_problems.md): o termo cointegrante
# (add_coint_term: Johansen nos 120 primeiros meses, padronizado com media e
# desvio da amostra inteira) foi removido. Nenhum script o usava: o ECM do 03
# usa X(t) e e(t) em nivel (Eq. 3 de Colombo e Pelagatti).

# _____________________________________________________________________________
# ___________________________________ EXPORT __________________________________
# _____________________________________________________________________________

if(!dir.exists("outputs")) dir.create("outputs")

for (cc in exchange_ccodes) {
  df_name <- paste0(cc, "_panel_data")
  if(exists(df_name, envir = .GlobalEnv)) {
    write_csv2(
      get(df_name),
      panel_path(cc, variant)    # outputs/<pais>_panel_data<sufixo da variante>.csv
    )
  }
}
